#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 16:47:13 2018

@author: mellis
"""
# Own modules
from load import load_coeff
from load import load_ener
from load import load_ham
from load import load_QM
from load import load_pos
from load import load_inp
from load import load_xyz
from load import load_tintf
from load import load_frc
from load import load_sigma
from load import load_K

from Plot import plot_utils
from Plot import plot_coeff
from Plot import plot_norm
from Plot import plot_QM
from Plot import plot_ham
from Plot import plot_ener
from Plot import plot_frc
from Plot import plot_tintf
from Plot import plot_pos
from Plot import plot_K

# External Modules
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import datetime
from matplotlib.widgets import MultiCursor
from collections import OrderedDict
import time
import re
import pandas as pd
import os


class Params(object):
    """
    Will store all the parameters needed to plot the graphs. Will also
    calculate parameters such as alpha etc...
    """
    def __init__(self, folder, reps, plot_params):
        self.folder = folder
        self.reps = reps
        self.plot_params = plot_params
        self.non_qlk_params = [i for i in self.plot_params if 'qm_r' not in i]
        self._set_coeff_params()
        self._correct_plot_params()
        self._get_inp_data()

        self._set_title()
        self.title = r""
#        self.title += "   (**tanh_width**)"

        self.colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
                       '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
                       'r', 'g', 'b']
        self.colors = [i for j in range(50) for i in self.colors]
        self._use_control = False
#        if self.num_reps == 1: self._use_control = False
        self.max_time = 'all'  # (in fs)
        self.min_time = 0  # (in fs)  NOT WORKING CAN ONLY USE 0
        self.quick_stride = 0  # (in fs)
        self.slow_stride = 0  # (in fs)
        self.dt = self.run_inp_params['NUCLEAR_TIMESTEP']
        self.atoms_to_plot = 'all'

        self.worst_reps = {}
        self.best_reps = {}
        self._fix_load_timings()

    def _get_inp_data(self):
        """
        Will recursively retrieve all the data from the run.inp file.
        """
        with open(self.folder+'run.inp', 'r') as f:
            inp_file = f.read().split('\n')
        self.nested_inp_params = {}
        self.run_inp_params = {}
        load_inp.parse_inp_file(self.nested_inp_params,
                                self.run_inp_params,
                                inp_file)
        self.atoms_per_site = self.run_inp_params['NUMBER_ATOMS_PER_SITE']
        self.num_states = self.run_inp_params['NUMBER_DIABATIC_STATES']
        self.num_active_atoms = self.atoms_per_site * self.num_states

    def _fix_load_timings(self):
        """
        Will convert the times to load into step numbers to pass into the
        load_* functions.
        """
        dt = self.run_inp_params['NUCLEAR_TIMESTEP']

        if type(self.max_time) == str:
            if self.max_time == 'all':
                self.max_step = 'all'
            else:
                raise SystemExit("""Sorry I don't recognise the setting for
max_time, you can use all, or specify a maximum time in fs.""")
        else:
            self.max_step = int(self.max_time/dt)
        self.min_time = int(self.min_time/dt)
        self.quick_stride = int(self.quick_stride/dt)
        self.slow_stride = int(self.slow_stride/dt)

        if type(self.max_step) != str and self.max_step < 1:
            self.max_step = 1
        if self.min_time < 0:
            self.min_time = 0
        if self.quick_stride < 1:
            self.quick_stride = 1
        if self.slow_stride < 1:
            self.slow_stride = 1

    def _set_title(self):
        """
        Will set the title of the plot according to the parameters given.
        """
        params_convert = {'|u|^2': 'Diab Coeffs',
                          '|c|^2': 'Adiab Coeffs',
                          'qm_t': "Quantum Momentum",
                          "adiab_state": "Adiabatic States",
                          "norm": "Norm",
                          'site_ener': 'site energy differences',
                          "qm_r": "Quantum Momentum",
                          "fl_fk": "history force state difference",
                          "ylk/sum(ylk)": "Rlk denominator",
                          'blank': "",
                          "energy_cons": "Energy Conservation",
                          "coup": r"H$_{12}$",
                          "force": "Nuc. Frc", }
        name_plot_params = self.plot_params[:]
        if 'qm_r' in name_plot_params and 'qm_t' in name_plot_params:
            name_plot_params.remove("qm_r")
        if len(name_plot_params) == 1:
            params_joined = str(params_convert.get(name_plot_params[0]))
        else:
            tmpStr = ', '.join([str(params_convert.get(i))
                                for i in name_plot_params[:-1]])
            params_joined = tmpStr + " and " \
                + str(params_convert.get(name_plot_params[-1]))
        params_joined = params_joined.replace("None", "???")

        if self.plot_params[0] != "qm_r":
            self.title = r"Evolution of %s under " % (params_joined) + \
                         "**CT/Eh** for **irep** replica [H$_{12}$ =" + \
                         "**coup** meV]"
        else:
            self.title = r"Spatial distribution of QM -data " + \
                        "from **irep** replica [H$_{12}$ = **coup** meV]"

    # Will get the number of replicas and the transparency of the lines.
    def _get_transparency(self):
        """
        Will change the alpha value depending on how many reps were used in the
        simulation.
        """
        alphas = {1: 1, 2: 0.8, 3: 0.6,
                  10: 0.4, 50: 0.3, 100: 0.2,
                  200: 0.15, 500: 0.1, 10000: 0.05}
        self.all_alphas = {}
        keys = sorted(alphas.keys())
        for i in range(len(alphas)-1):
            curr_key = keys[i]
            next_key = keys[i+1]
            fit = np.polyfit([curr_key, next_key],
                             [alphas[curr_key],
                              alphas[next_key]],
                             1)
            for i in range(curr_key, next_key):
                self.all_alphas[i] = np.polyval(fit, i)
        self.alpha = self.all_alphas[self.num_reps]

    def _correct_plot_params(self):
        """
        Will hopefully correct typos in the plotting parameters.
        """
        self.plot_params = pd.Series([i.strip().lower()
                                      for i in self.plot_params])
        self.plot_params = list(self.plot_params.unique())

        for short_name in ['adiab_state', 'coup']:
            self.plot_params = [short_name if short_name in i else i
                                for i in self.plot_params]

        self.plot_paramsC = [i for i in self.plot_params
                             if any(i == j for j in ('|u|^2', '|c|^2'))]

    def _set_coeff_params(self):
        self.check_control_coeff = {}
        self.all_coeff_lines = {}
        self.avg_coeff_lines = {}
        self.coeff_widg_axes = {}
        self.coeff_plot_axes = {}


class LoadData(Params):
    """
    Will load all data in the folder. Takes 3 inputs:
        * folder      = folder in which to look for data
        * reps        = which replicas to plot
                        (can be 'all' or list/range of integers)
        * plot_params = which parameters to plot (can be a list of strings)
                        valid parameters:
                            - 'norm'
                            - '|C|^2'
                            - '|u|^2'
                            - 'adiab_state'
                            - 'qm'
                            - 'site_ener'
                            - 'fk_fl' or 'fl_fk'
        * avg_on      = Whether to average the data or not (default = True)

    NOTE: Should be in plot_utils
    """
    # Decide which things to load when the params are inputted
    #     (avoids loading multiple times of some dependencies)
    dependencies = {'qm_r':        ['pos', 'qm'],
                    'pos':         ['pos'],
                    '|c|^2':       ['|c|^2'],
                    'qm_t':        ['qm'],
                    'rlk':        ['rlk'],
                    'site_ener':   ['ham'],
                    'coup':        ['ham'],
                    '|u|^2':       ['|u|^2'],
                    'energy_cons': ['tot_ener'],
                    'adiab_state': ['ad_ener'],
                    'ylk/sum(ylk)': ['fl_fk', '|c|^2'],
                    'fl_fk':       ['fl_fk'],
                    'norm':        ['|u|^2'],
                    'forces':      ['force'],
                    'alpha':       ['qm', 'rlk'],
                    'pos_sigma':   ['sigma', 'pos'],
                    'sum(ylk)':    ['fl_fk', '|c|^2'],
                    'k':           ['k'],
                    }

    def __init__(self, folder, reps, plot_params='all', avg_on=True):
        if folder[-1] != '/':
            folder += '/'
        self.folder = folder
        self.reps = reps
        self.plot_params = plot_params
        self.avg_on = avg_on
        Params.__init__(self, folder, reps, plot_params)

        self.decideDependencies()
        self.load_timings = OrderedDict()

        self.load_all_ham_data()
        self.load_all_di_coeffs()
        self.load_all_ad_coeffs()
        self.load_ad_ener()
        self.load_rlk()
        self.load_qm()
        self.load_pos()
        self.load_hist_f()
        self.load_tot_ener()
        self.load_force()
        self.load_sigmas()
        self.load_k()

        self.calc_alpha()
        self._get_transparency()

        if self.avg_on:
            self._average_data()

        self.print_timing_info(self.load_timings,
                               "Timing Data for Reading Data")

    def decideDependencies(self):
        """
        Will decide which parameters to load from the inputs, i.e. if the |u|^2
        , and norm are requested then it makes no sense to load |u|^2 twice...
        """
        self.load_params = []
        for i in self.plot_params:
            for j in self.dependencies[i.lower()]:
                self.load_params.append(j)
        self.load_params = list(set(self.load_params))

        # Print some feedback
        for i in self.load_params:
            print("Loading %s" % i)

    def __load_ham(self, stride, max_step):
        """
        Will load the hamiltonian (should always use load_all_ham_data instead
        of this).
        """
        self.load_timings['H'] = time.time()
        self.all_ham_data = load_ham.load_all_ham_in_folder(
                                                       self.folder,
                                                       reps=self.reps,
                                                       max_step=max_step,
                                                       min_step=self.min_time,
                                                       stride=stride
                                                           )
        allLens = [len(self.all_ham_data[i][0]) for i in self.all_ham_data]
        if all(allLens):
            # Average data
            self.avg_ham_data = plot_utils.avg_H_data_dict(self.all_ham_data)
            tmp = plot_utils.get_coup_data(self.avg_ham_data, 'avg_ham')
            self.avg_site_ener, self.avg_couplings = tmp[0], tmp[1]
            self.avg_avg_couplings, self.Stimesteps = tmp[2], tmp[3]
            self.all_site_ener = [plot_utils.get_coup_data(self.all_ham_data,
                                                           ham_key)
                                  for ham_key in self.all_ham_data]
            self.all_site_ener = [[i[0], i[3]] for i in self.all_site_ener]

            self.load_timings['H'] = time.time() - self.load_timings['H']

            # Find metadata
            self.num_ham_steps = len(self.avg_ham_data['avg_ham'][0])
            self.num_states = len(self.all_ham_data[
                                  list(self.all_ham_data.keys())[0]][0][0])
            self.num_reps = len(self.all_ham_data)
            return True
        return False

    def load_all_ham_data(self):
        """
        Will load all the hamiltonian data that can be foun_set_Qlk_controld in
        the folder specified (dependent on which reps are requested)

        N.B. This is always called as it is normally outputted and it is used
        to find metadata such as how many steps and reps have been loaded. This
        can be much better though. I will improve it when I have some time!
        """
        print_step = self.nested_inp_params['FORCE_EVAL']['MIXED']['ADIABATIC']['PRINT']['HAMILTONIAN']['EACH']['MD'][0]
        if type(self.max_step) == str:
            max_step = self.max_step
        else:
            max_step = int(self.max_step/print_step)

        if 'ham' in self.load_params:
            exitCode = self.__load_ham(self.quick_stride, max_step)
        else:
            if type(max_step) == str or \
                                   (type(max_step) == int and max_step <= 100):
                stride = 1
            else:
                ham_file = [i for i in os.listdir(self.folder)
                            if 'hamil' in i][0]
                with open(self.folder + ham_file, 'r') as f:
                    ltxt = f.read().split('\n')
                tmp = load_xyz.get_xyz_step_metadata(ltxt, ham_file)
                _, _, lines_in_step, num_title_lines = tmp
                num_steps = len(ltxt) / (lines_in_step)
                stride = int(num_steps / 100)
                if stride < 1:
                    stride = 1
                max_step = num_steps

            exitCode = self.__load_ham(stride, max_step)

        self.coupling = "?"
        if exitCode:
            self._get_coupling()
            self.coupling = "%.2g" % self.coupling

    def load_all_di_coeffs(self):
        """
        Loads all the diabatic coefficients, no input.
        Saves diabatic coeffs as self.all_Dcoeff_data
        """
        if '|u|^2' in self.load_params:
            self.load_timings['di coeff'] = time.time()
            print_step = self.nested_inp_params['FORCE_EVAL']['MIXED']['ADIABATIC']['PRINT']['COEFFICIENTS']['EACH']['MD'][0]

            if type(self.max_step) == str:
                max_step = self.max_step
            else:
                max_step = int(self.max_step/print_step)

            tmpList = ['xyz', 'coeff']
            self.all_Dcoeff_data = load_coeff.load_all_coeff_in_folder(
                                             self.folder,
                                             filename_must_contain=tmpList,
                                             filename_must_not_contain=['ad'],
                                             reps=self.reps,
                                             max_step=max_step,
                                             min_step=self.min_time,
                                             stride=self.quick_stride
                                                                      )
            if not self.all_Dcoeff_data:
                raise SystemExit("Sorry I can't find any coeff data in " +
                                 "folder:\n\n\t$s" % self.folder)
            self.load_timings['di coeff'] = time.time() - \
                self.load_timings['di coeff']

            # Find metadata
            Keys = list(self.all_Dcoeff_data.keys())
            self.num_reps = len(self.all_Dcoeff_data)
            self.num_di_coeff_steps = len(self.all_Dcoeff_data[Keys[0]][0])
            self.num_states = len(self.all_Dcoeff_data[Keys[0]][3][0])

    def load_all_ad_coeffs(self):
        """
        Loads all the adiabatic coefficients, no input.
        Saves adiabatic coeffs as self.all_Acoeff_data.
        """
        if '|c|^2' in self.load_params:
            try:
                print_step = self.nested_inp_params['MOTION']['CTMQC']['PRINT']['AD_COEFF']['EACH']['MD'][0]
            except KeyError:
                print_step = self.nested_inp_params['FORCE_EVAL']['MIXED']['ADIABATIC']['PRINT']['COEFFICIENTS']['EACH']['MD'][0]

            if type(self.max_step) == str:
                max_step = self.max_step
            else:
                max_step = int(self.max_step/print_step)
            self.load_timings['ad coeff'] = time.time()
            self.all_Acoeff_data = plot_utils.load_Acoeff_data(
                                                       self.folder,
                                                       self.reps,
                                                       self.all_ham_data,
                                                       max_step=max_step,
                                                       min_step=self.min_time,
                                                       stride=self.quick_stride
                                                              )
            self.load_timings['ad coeff'] = time.time() - \
                self.load_timings['ad coeff']

            # Find metadata
            Keys = list(self.all_Acoeff_data.keys())
            self.num_reps = len(self.all_Acoeff_data)
            self.num_ad_coeff_steps = len(self.all_Acoeff_data[Keys[0]][0])
            self.num_states = len(self.all_Acoeff_data[Keys[0]][3][0])

    def load_ad_ener(self):
        """
        Loads the adiabatic energy
        """
        if 'ad_ener' in self.load_params:
            self.load_timings['adiab ener'] = time.time()
            if type(self.max_step) == str:
                max_step = self.max_step
            else:
                max_step = self.max_step*self.dt

            self.all_ad_ener_data = load_ener.load_all_ener_ad(
                                               self.folder,
                                               reps=self.reps,
                                               max_time=max_step,
                                               min_time=self.min_time * self.dt
                                                              )
            if not self.all_ad_ener_data:
                raise IOError("Can't find any data, please check folder.")
            self.load_timings['adiab ener'] = time.time() - \
                self.load_timings['adiab ener']

            # Find metadata
            Keys = list(self.all_ad_ener_data.keys())
            cols = self.all_ad_ener_data[Keys[0]].columns
            self.num_reps = len(Keys)
            self.num_ad_ener_steps = len(self.all_ad_ener_data[Keys[0]])
            self.num_states = len([i for i in cols if 'State' in i])

    def load_tot_ener(self):
        """
        Loads the total energy files. These files contain the kinetic energy,
        temperature, potential energy, total energy and CPU time taken.
        """
        if 'tot_ener' in self.load_params:
            self.load_timings['tot ener'] = time.time()

            self.all_tot_ener = load_ener.load_all_ener_dat(
                                                        folder=self.folder,
                                                        reps=self.reps,
                                                        max_time=self.max_time
                                                           )

            self.load_timings['tot ener'] = time.time() - \
                self.load_timings['tot ener']

            # Find metadata
            Keys = list(self.all_tot_ener.keys())
            self.num_reps = len(Keys)
            self.num_tot_ener_steps = len(self.all_tot_ener[Keys[0]])

    def load_sigmas(self):
        """
        Will load the sigma data from a simulation
        """
        if 'sigma' in self.load_params:
            self.load_timings['sigma'] = time.time()
            if type(self.max_step) == str:
                max_step = self.max_step
            else:
                max_step = self.max_step*self.dt

            self.all_sigma = load_sigma.load_all_sigma_in_folder(
                                                       folder=self.folder,
                                                       reps=self.reps,
                                                       max_step=max_step,
                                                       min_step=self.min_time,
                                                       stride=self.quick_stride
                                                           )
            self.load_timings['sigma'] = time.time() - \
                self.load_timings['sigma']

            # Find metadata
            Keys = list(self.all_sigma.keys())
            self.num_reps = len(Keys)

    def load_k(self):
        """
        Will load the K data from a simulation (the k data is the C_ctmqc)
        """
        if 'k' in self.load_params:
            self.load_timings['K'] = time.time()
            if type(self.max_step) == str:
                max_step = self.max_step
            else:
                max_step = self.max_step*self.dt

            self.all_K = load_K.load_all_K_in_folder(
                                                       folder=self.folder,
                                                       reps=self.reps,
                                                       max_step=max_step,
                                                       min_step=self.min_time,
                                                       stride=self.quick_stride
                                                           )
            self.load_timings['K'] = time.time() - \
                self.load_timings['K']

            # Find metadata
            Keys = list(self.all_K.keys())
            self.num_reps = len(Keys)

    def calc_alpha(self):
        """
        Will calculate the alpha variable (qlk = alpha - rlk)
        """
        if any('alpha' in param for param in self.plot_params):
            self.load_timings['alpha'] = time.time()

            self.all_alpha = {}
            for i, qm_key in enumerate(self.all_Qlk_data):
                qlk = self.all_Qlk_data[qm_key]
                alpha = qlk[0][0] + self.Rlk_data[0][0]
                alpha = ((alpha, qlk[0][1]), qlk[1])
                self.all_alpha[i] = alpha

            self.load_timings['alpha'] = time.time() - \
                self.load_timings['alpha']

    def load_qm(self):
        """
        Will load the quantum momentum file into the format in load_QM.
        """
        if 'qm' in self.load_params:
            self.load_timings['QM'] = time.time()
            print_step = self.nested_inp_params['MOTION']['CTMQC']['PRINT']['QM']['EACH']['MD'][0]
            if type(self.max_step) == str:
                max_step = self.max_step
            else:
                max_step = int(self.max_step/print_step)
            self.all_Qlk_data = load_QM.load_all_Qlk_in_folder(
                                                        self.folder,
                                                        reps=self.reps,
                                                        max_step=max_step,
                                                        min_step=self.min_time,
                                                        stride=self.slow_stride
                                                               )
            self.load_timings['QM'] = time.time() - self.load_timings['QM']

            # Find metadata
            Keys = list(self.all_Qlk_data.keys())
            self.num_reps = len(Keys)
            cols = self.all_Qlk_data[Keys[0]][0][1]
            self.num_qm_steps = len(cols)
            self.num_active_atoms = max(cols[0, :, 0])
            self.num_states = max(cols[0, :, 2])

    def load_rlk(self):
        """
        Will load the Rlk file into the format in load_QM.
        """
        if 'rlk' in self.load_params:
            self.load_timings['Rlk'] = time.time()
            print_step = self.nested_inp_params['MOTION']['CTMQC']['PRINT']['RLK']['EACH']['MD'][0]
            if type(self.max_step) == str:
                max_step = self.max_step
            else:
                max_step = int(self.max_step/print_step)
            self.Rlk_data = load_QM.load_all_Rlk_in_folder(
                                                        self.folder,
                                                        reps=self.reps,
                                                        max_step=max_step,
                                                        min_step=self.min_time,
                                                        stride=self.slow_stride
                                                               )
            self.load_timings['Rlk'] = time.time() - self.load_timings['Rlk']

            # Find metadata
            cols = self.Rlk_data[0][1]
            self.num_rlk_steps = len(cols)
            self.num_active_atoms = max(cols[0, :, 0])
            self.num_states = max(cols[0, :, 2])

    def load_pos(self):
        """
        Will load the positions
        """
        if 'pos' in self.load_params:
            self.load_timings['positions'] = time.time()
            print_step = self.nested_inp_params['MOTION']['PRINT']['TRAJECTORY']['EACH']['MD'][0]
            if type(self.max_step) == str:
                max_step = self.max_step
            else:
                max_step = int(self.max_step/print_step)

            self.all_pos_data = load_pos.load_all_pos_in_folder(
                                                       self.folder,
                                                       reps=self.reps,
                                                       max_step=max_step,
                                                       min_step=self.min_time,
                                                       stride=self.slow_stride
                                                               )
            self.load_timings['positions'] = time.time() - \
                self.load_timings['positions']

            # Find metadata
            Keys = list(self.all_pos_data.keys())
            self.num_reps = len(Keys)
            cols = self.all_pos_data[Keys[0]][0][1]
            self.num_pos_steps = len(cols)
            self.num_atoms = len(cols[0])

    def load_hist_f(self):
        """
        Will load all the time-integrated history forces in a folder.
        """
        if 'fl_fk' in self.load_params:
            self.load_timings['history forces'] = time.time()
            print_step = self.nested_inp_params['MOTION']['CTMQC']['PRINT']['T_INT_FORCE']['EACH']['MD'][0]

            if type(self.max_step) == str:
                max_step = self.max_step
            else:
                max_step = int(self.max_step/print_step)

            self.all_tintf_data = load_tintf.load_all_tintf_in_folder(
                                                        self.folder,
                                                        reps=self.reps,
                                                        max_step=max_step,
                                                        min_step=self.min_time,
                                                        stride=self.slow_stride
                                                                     )
            self.load_timings['history forces'] = time.time() - \
                self.load_timings['history forces']

            # Find metadata
            Keys = list(self.all_tintf_data.keys())
            self.num_reps = len(Keys)
            cols = self.all_tintf_data[Keys[0]][0][1]
            self.num_histf_steps = len(cols)
            self.num_active_atoms = sum(cols[0, :, 1] == '1')
            self.num_states = len(set(cols[0, :, 1]))

    def load_force(self):
        """
        Will load all the nuclear forces.
        """
        if 'force' in self.load_params:
            self.load_timings['forces'] = time.time()
            print_step = self.nested_inp_params['MOTION']['PRINT']['FORCES']['EACH']['MD'][0]
            if type(self.max_step) == str:
                max_step = self.max_step
            else:
                max_step = int(self.max_step/print_step)

            self.all_frc_data = load_frc.load_all_frc_in_folder(
                                                        self.folder,
                                                        reps=self.reps,
                                                        max_step=max_step,
                                                        min_step=self.min_time,
                                                        stride=self.slow_stride
                                                               )
            self.load_timings['forces'] = time.time() - \
                self.load_timings['forces']

            # Find metadata
            Keys = list(self.all_frc_data.keys())
            self.num_reps = len(Keys)
            self.num_frc_steps = len(self.all_frc_data[Keys[0]][0])
            self.num_atoms = len(self.all_frc_data[Keys[0]][0][0][0])

    def _get_coupling(self):
        """
        Will get the average coupling and print it to the console
        """
        self.coupling = np.mean(
                               np.abs(self.avg_ham_data['avg_ham'][0][:, 0, 1])
                               )
        self.coupling *= 27e3  # convert to meV

    # Will average the coefficient data (ham is averaged by default)
    def _average_data(self):
        """
        Will average the coefficient data and save as new arrays
        """
        self.load_timings['Averaging: '] = OrderedDict()
        if "|u|^2" in self.load_params:
            self.load_timings['Averaging: ']['di coeff'] = time.time()
            self.all_Dcoeff_data_avg = plot_utils.avg_coeff_data(
                                                          self.all_Dcoeff_data
                                                                )
            self.load_timings['Averaging: ']['di coeff'] = time.time() - \
                self.load_timings['Averaging: ']['di coeff']

        if "|c|^2" in self.load_params:
            self.load_timings['Averaging: ']['ad coeff'] = time.time()
            self.all_Acoeff_data_avg = plot_utils.avg_coeff_data(
                                                           self.all_Acoeff_data
                                                                )
            self.load_timings['Averaging: ']['ad coeff'] = time.time() - \
                self.load_timings['Averaging: ']['ad coeff']

        if 'ad_ener' in self.load_params:
            self.load_timings['Averaging: ']['adiab_ener'] = time.time()
            self.all_ad_ener_data_avg = \
                plot_utils.avg_E_data_dict(self.all_ad_ener_data)
            self.load_timings['Averaging: ']['adiab_ener'] = \
                time.time() - self.load_timings['Averaging: ']['adiab_ener']

        if 'pos' in self.load_params:
            self.load_timings['Averaging: ']['pos'] = time.time()
            self.avg_pos_data = plot_utils.avg_pos_data(self.all_pos_data)
            self.load_timings['Averaging: ']['pos'] = time.time() - \
                self.load_timings['Averaging: ']['pos']

        if 'qm' in self.load_params:
            self.load_timings['Averaging: ']['qm'] = time.time()
            self.avg_Qlk_data = plot_utils.avg_Qlk_data(self.all_Qlk_data)
            self.load_timings['Averaging: ']['qm'] = time.time() - \
                self.load_timings['Averaging: ']['qm']

        if 'fl_fk' in self.load_params:
            self.load_timings['Averaging: ']['history forces'] = time.time()
            self.sum_tintf_data = plot_utils.sum_hist_f_data(
                                                            self.all_tintf_data
                                                            )
            if 'ylk/sum(ylk)' in self.plot_params:
                self.sum_tintf_CC_data = plot_utils.sum_hist_f_CC_data(
                                                           self.all_tintf_data,
                                                           self.all_Acoeff_data
                                                                      )
            if 'sum(ylk)' in self.plot_params:
                self.sum_ylk = plot_utils.sum_Ylk_data(self.all_tintf_data,
                                                       self.all_Acoeff_data
                                                       )
            self.load_timings['Averaging: ']['history forces'] = time.time() -\
                self.load_timings['Averaging: ']['history forces']

        if 'force' in self.load_params:
            self.load_timings['Averaging: ']['force'] = time.time()
            self.avg_frc_data = plot_utils.avg_pos_data(self.all_frc_data)
            # Rename the average key (using average_pos function)
            self.avg_frc_data['avg_frc'] = self.avg_frc_data['avg_pos']
            del self.avg_frc_data['avg_pos']
            self.load_timings['Averaging: ']['force'] = time.time() - \
                self.load_timings['Averaging: ']['force']

        if 'tot_ener' in self.load_params:
            self.load_timings['Averaging: ']['Tot. Energy'] = time.time()
            self.tot_ener_mean = pd.concat(
                                           self.all_tot_ener.values()
                                          ).groupby('Step').mean()
            self.load_timings['Averaging: ']['Tot. Energy'] = time.time() -  \
                self.load_timings['Averaging: ']['Tot. Energy']

    def print_timing_info(self, timing_dict, title=""):
        """
        Will print any timings info.

        Inputs:
            * timing_dict  =>  the dictionary containing timing data
            * title        =>  the title of the timing data
        """
        max_len = 50
        print(" "+"-"*max_len)
        num_spaces = (len(title)-1)*2
        print("|"+" "*int(max_len - num_spaces) + title +
              " "*int(max_len - num_spaces) + " |")
        print("|"+" "*max_len+"|")
        plot_utils.print_timings(timing_dict, max_len=max_len)
        print(" "+"-"*max_len)


class Plot(LoadData, Params, plot_norm.Plot_Norm, plot_coeff.Plot_Coeff,
           plot_ener.Adiab_States, plot_ham.Coupling, plot_QM.QM_R,
           plot_QM.QM_t, plot_ham.Site_Ener, plot_tintf.fl_fk,
           plot_tintf.fl_fk_CC, plot_ener.Energy_Cons, plot_frc.Plot_Frc,
           plot_QM.Rlk, plot_QM.Alpha, plot_pos.PlotPos, plot_pos.PlotPosSig,
           plot_tintf.sumYlk):
    """
    Will handle plotting of (hopefully) any parameters. Pass a list of string
    with the parameters that are to be plotted. E.g. Plot(['|u|^2', '|C|^2'])
    adiab_state and this class should plot them

    Inputs:
        plot_params    =>  A list containing the parameters needing plotting.
                           Possible parameters are:
                               * |u|^2        = Diabatic populations
                               * |C|^2        = Adiabatic populations
                               * norm         = The norm of the diabatic coeffs
                               * qm           = The Quantum Momentum
                               * adiab_state = The adiabatic energy levels
                               * site_ener    = The site energies vs time
                               * fl_fk        = The difference in history force
                                                states.

        folder         =>  The folder containing the data
        reps           =>  Which replica numbers to plot (can be 'all')
        plot           =>  Can be True or False or 'close' to create the graphs
                               but not display.
                             If True will plot to screen,
                                 if False will just analyse data.
                             If 'close' then will not display the graph but
                                 hold it in memory (as self.f).
    """

    def __init__(self, plot_params, folder, reps, plot):
        print("Starting on folder ", folder)
        self.plot_params = plot_params
        self._correct_plot_params()
#        Params.__init__(self, folder, reps, self.plot_params)
        LoadData.__init__(self, folder, reps, self.plot_params)
        if self.atoms_to_plot == 'all':
            self.atoms_to_plot = range(1, self.num_active_atoms+1)

        self.reps = reps
        self.xlabel = "Time (fs)"
        self.plot_info = {}
        if type(plot) == str:
            if plot.lower() == 'close':
                self.plot = True
            self._create_ax_fig_layout(close=True)
        else:
            self.plot = plot
            self._create_ax_fig_layout(close=False)

        # REMOVE WHEN ALL CLASS PLOTS ARE CREATED #
        self.plot_all_reps = True
        self.avg_on = True
        self.fill_between = True
        #############################################

        self.plot_blank()

        # Quantum Momentum
        if 'qm_r' in self.plot_params:
            self.mQMRPlot = plot_QM.QM_R.__init__(self,
                                                  self.axes['qm_r'])
        if 'qm_t' in self.plot_params:
            self.mQMTPlot = plot_QM.QM_t.__init__(self,
                                                  self.axes['qm_t'])
        if 'rlk' in self.plot_params:
            self.mRlkPlot = plot_QM.Rlk.__init__(self,
                                                 self.axes['rlk'])
        if 'alpha' in self.plot_params:
            self.mAlphaPlot = plot_QM.Alpha.__init__(self, self.axes['alpha'])

        if 'k' in self.plot_params:
            self.mKPlot = plot_K.K.__init__(self, self.axes['k'])

        # Coefficients
        if 'norm' in self.plot_params:
            self.mNormPlot = plot_norm.Plot_Norm.__init__(self,
                                                          self.axes['norm'])
        if '|u|^2' in self.plot_params:
            self.mDCoeffPlot = plot_coeff.Plot_Coeff.__init__(
                                                             self,
                                                             self.axes['|u|^2']
                                                             )
        if '|c|^2' in self.plot_params:
            self.mACoeffPlot = plot_coeff.Plot_Coeff.__init__(
                                                             self,
                                                             self.axes['|c|^2']
                                                             )

        # Energies
        if 'adiab_state' in self.plot_params:
            self.mAStatesPlot = plot_ener.Adiab_States.__init__(
                                                       self,
                                                       self.axes['adiab_state']
                                                               )
        if 'coup' in self.plot_params:
            self.mCoupPlot = plot_ham.Coupling.__init__(self,
                                                        self.axes['coup'])
        if 'site_ener' in self.plot_params:
            self.mSiteEnerPlot = plot_ham.Site_Ener.__init__(
                                                         self,
                                                         self.axes['site_ener']
                                                            )
        if 'energy_cons' in self.plot_params:
            self.mEnerConsPlot = plot_ener.Energy_Cons.__init__(
                                                       self,
                                                       self.axes['energy_cons']
                                                               )

        # Forces
        if 'ylk/sum(ylk)' in self.plot_params:
            self.mTintfCCPlot = plot_tintf.fl_fk_CC.__init__(
                                                      self,
                                                      self.axes['ylk/sum(ylk)']
                                                            )
        if 'sum(ylk)' in self.plot_params:
            self.mTintfCCPlot = plot_tintf.sumYlk.__init__(
                                                          self,
                                                          self.axes['sum(ylk)']
                                                            )
        if 'fl_fk' in self.plot_params:
            self.mTintfPlot = plot_tintf.fl_fk.__init__(self,
                                                        self.axes['fl_fk'])

        if 'force' in self.plot_params:
            self.mForcePlot = plot_frc.Plot_Frc.__init__(
                                                         self,
                                                         self.axes['force']
                                                         )
        # Positions
        if 'pos' in self.plot_params:
            self.mPosPlot = plot_pos.PlotPos.__init__(self,
                                                      self.axes['pos'])
        if 'pos_sigma' in self.plot_params:
            self.mPosSigPlot = plot_pos.PlotPosSig.__init__(
                                                        self,
                                                        self.axes['pos_sigma'],
                                                            )

        # Finish up
        if self.plot:
            self.__finalise()
        self.print_final_info()

        if type(plot) == str and plot == 'close':
            plt.close('all')

    def _clean_widget_axes(self, axis):
        """
        Will remove spines and ticks from an axis

        Inputs:
            * axis  =>  a plt.axis

        Outpus:
            axis (plt.axis)
        """
        # Remove the ticks
        axis.set_xticks([])
        axis.set_yticks([])
        for side in ['top', 'bottom', 'left', 'right']:
            axis.spines[side].set_visible(False)

        return axis

    def _N_pane_special_case(self, num_panes, ax_name, plot_params):
        """
        Will create a special case control panel for axes with 2 panes. The
        slightly odd ax_name input is to enable multiple substring matches
        in a string. It is a list of substrings that are in the axis name.

        E.g. if ax_name = ['fl', 'fk'] all axes with both fl and fk in their
             name will be found. If ax_name = ['qm_t'] then all axes with qm_t
             in their name will be found etc...
        """
        if all([any(name in j for j in plot_params) for name in ax_name]):
            ax_name = '_'.join(ax_name).strip('_')
            if ax_name not in plot_params:
                return
            widg_ax = []
            for i in range(num_panes):
                widg_ax.append(plt.subplot2grid(
                              (len(plot_params)*num_panes, 7),
                              (plot_params.index(ax_name)*num_panes+i, 0),
                              colspan=1)
                                               )
            self.axes[ax_name] = [widg_ax,
                                  plt.subplot2grid(
                                          (len(plot_params), 7),
                                          (plot_params.index(ax_name), 1),
                                          colspan=6
                                           )]
            for i in range(num_panes):
                self.axes[ax_name][0][i] = self._clean_widget_axes(
                                                       self.axes[ax_name][0][i]
                                                                  )

    # Decides what arrangement of axes to use
    def _create_ax_fig_layout(self, close=False):
        """
        Will create the layout for the plots and assign an axis to each
        plotting parameter. This will create the self.f and self.axes
        variables.

        The self.axes variable is a dictionary with the plot parameter as a key
        and the axis that has been assigned to it as the value.
        """
        if self.plot:
            self.f = plt.figure()
            if close:
                plt.close('all')
            self.axes = OrderedDict()

            plot_params = self.plot_params[:]

            # Give structure to the OrderedDict
            for i in plot_params:
                self.axes[i] = ''

            if len(plot_params) <= 4:

                for i, param in enumerate(plot_params):
                    a = []
                    if self._use_control:
                        a.append(plt.subplot2grid((len(plot_params), 7),
                                                  (i, 0),
                                                  colspan=1))
                        a.append(plt.subplot2grid((len(plot_params), 7),
                                                  (i, 1),
                                                  colspan=6))
                        a[0] = self._clean_widget_axes(a[0])
                    else:
                        a.append('')
                        a.append(plt.subplot2grid((len(plot_params), 1),
                                                  (i, 0)))
                    self.axes[param] = a
            else:
                plt.close()
                raise SystemExit("Sorry I don't have any way to handle more"
                                 + " than 3 plots at the same time yet!")
            if self._use_control:
                self._N_pane_special_case(3, ['fl', 'fk', 'cc'], plot_params)
                self._N_pane_special_case(3, ['fl', 'fk'], plot_params)
                self._N_pane_special_case(3, ['qm_t'], plot_params)
                self._N_pane_special_case(3, ['alpha'], plot_params)
        else:
            self.axes = {i: '' for i in self.plot_params}

    def plot_blank(self):
        """
        A function to quickly plot on a blank axis. This can be editted by
        anyone wanting to quickly plot something without making it
        something that is officialy plotted. Things should be plotted on the
        'blank' axis.
        """

        if 'blank' in self.plot_params:
            widgs, ax = self.axes['blank']
            hdata, _, timesteps = self.avg_ham_data['avg_ham']
            all_Es = np.array([plot_utils.calc_U_matrix(H)[0] for H in hdata])
            ax.plot(timesteps, all_Es[:, 0])
            ax.plot(timesteps, all_Es[:, 1])
            ax.set_ylabel(r"E$_{i}$ [Ha]")
            ax.set_xlabel(r"Time [fs]")

    def print_final_info(self):
        """
        Will print some information about the data just plotted. This will
        include stuff from the input file such as whether the commutator was
        used, whether the run was Ehrenfest or not, the data/time of the plot
        , how many replicas were used and how many molecules.
        """
        if self.run_inp_params['METHOD_PROPAGATION'] == 'CTMQC':
            # Build the strs list
            str_sections = {'Vitals': [], 'Ehrenfest': [], 'CTMQC': []}
            for i in self.plot_info:
                str_sections[i] = self.plot_info[i]

            # Section 1
            strs = str_sections['Vitals']
            folder = self.folder.strip('/')
            folder = folder[folder.strip('/').rfind('/'):].strip('/')
            strs.append("Folder = " + folder)
            strs.append(datetime.datetime.strftime(datetime.datetime.now(),
                                                   "Time of plot = %H:%M"))
            strs.append(datetime.datetime.strftime(datetime.datetime.now(),
                                                   "Date = %d/%m/%Y"))
            strs.append("dt  =  %.3f" % (
                                        self.run_inp_params['NUCLEAR_TIMESTEP']
                                        ))
            strs.append("Scaling Factor = %.3g" % (
                                          self.run_inp_params['SCALING_FACTOR']
                                                  ))
            strs.append("Coupling  =  %s meV" % self.coupling)

            # Section 2
            if not self.run_inp_params['USE_QM']:
                strs = str_sections['Ehrenfest']
                if not self.run_inp_params.get('FAST_EHRENFEST'):
                    strs.append("Commutator was ON")
                else:
                    strs.append("Commutator was OFF")
                strs.append("Num Replicas = %i" % self.num_reps)

            else:
                strs = str_sections['CTMQC']
                if not self.run_inp_params.get('FAST_EHRENFEST'):
                    strs.append("Commutator was ON")
                else:
                    strs.append("Commutator was OFF")
                strs.append("Num Replicas = %i" % self.num_reps)
                strs.append("Initial Width =" +
                            " %.3g" % self.run_inp_params['INITIAL_SIGMA'])
                strs.append("Tanh Width = " +
                            "%.2g" % (self.run_inp_params['TANH_WIDTH']))

            # Find max len of string and add a hash to the end of each line
            max_len_str = np.max([len(i) for j in str_sections
                                  for i in str_sections[j]]) + 4
            str_sections = {sect: [i + ' ' * (max_len_str-len(i)) + '#'
                                   for i in str_sections[sect]]
                            for sect in str_sections}

            # Do the printing
            tabN = 6
            tab = " "*tabN
            for i in str_sections:
                if not str_sections[i]:
                    continue
                strs = str_sections[i]
                print("\n"+'#'*tabN+'#'*max_len_str+'#')
                print(i + ' '*(max_len_str-len(i))+tab+'#')
                print('-'*len(i) + ' '*(max_len_str-len(i))+tab+'#')
                for line in strs:
                    print(tab+" "*max_len_str+'#')
                    print(tab+line)

                print(tab+" "*max_len_str+'#')
                print('#'*tabN+'#'*max_len_str+'#')

    def _fill_in_the_title(self):
        """
        Will replace certain words in the title with parameters used in the
        graph etc...
        """
        if self.run_inp_params['METHOD_PROPAGATION'] == 'CTMQC':
            # Rep num
            self.title = self.title.replace("**irep**", str(self.num_reps))
            if self.num_reps > 1:
                self.title = self.title.replace("replica", "replicas")

            # Using Comm
            if self.run_inp_params.get('FAST_EHRENFEST'):
                self.title = self.title.replace("**comm**", "off")
            else:
                self.title = self.title.replace("**comm**", "on")

            # Ehrenfest or CTMQC
            if self.run_inp_params['USE_QM']:
                self.title = self.title.replace("**CT/Eh**", "CTMQC")
                tmp = str(self.run_inp_params['TANH_WIDTH'])
                self.title = self.title.replace("**tanh_width**",
                                                "tanh width = " + tmp)
            else:
                self.title = self.title.replace("**CT/Eh**", "Ehrenfest")
        else:
            tmp = self.run_inp_params['METHOD_PROPAGATION'].title()
            tmp = tmp.replace("_", " ")
            self.title = "%s propagation CP2K" % (tmp)

        self.title = self.title.replace("**coup**", "%s" % self.coupling)
        self.title = re.sub(r"\*\*.*\*\*", "", self.title)
        self.title = re.sub(r"\( *\)", "", self.title)

    # Will finish off the plots
    def __finalise(self):
        """
        Will finish off the plots by adding necessary (communal) labels etc..
        e.g. will put the time (fs) label on the lowest x axis.
        """
        self._fill_in_the_title()
        if any([j in self.plot_params for j in ('|u|^2', '|c|^2')]):
            # Set legend
            if '|u|^2' in self.plot_params:
                num_states = len(self.all_Dcoeff_data_avg[3][0])
            elif '|c|^2' in self.plot_params:
                num_states = len(self.all_Acoeff_data_avg[3][0])
            elif 'adiab_state' in self.plot_params:
                num_states = len(self.state_cols_AS)
            elif "site_ener" in self.plot_params:
                num_states = len(self.avg_ham_data['avg_ham'][0][0])
            labels = ["State %i" % (i+1) for i in range(num_states)]
            patches = [mpatches.Patch(color=self.colors[i], label=lab)
                       for i, lab in enumerate(labels)]
            if num_states < 10:
                fit = np.polyfit([2, 10], [18, 10], 1)
                self.f.legend(handles=patches,
                              fontsize=np.polyval(fit, num_states),
                              labels=labels,
                              loc='upper left')
        self.f.suptitle(self.title, fontsize=20)

        # For all axes
        for ax in self.axes:
            AX = self.axes[ax][1]
            AX.spines['top'].set_visible(False)
            AX.spines['right'].set_visible(False)
            AX.spines['bottom'].set_visible(False)
            AX.spines['left'].set_visible(True)
            AX.grid('on', alpha=0.5)
            if type(self.max_time) != str and ax != 'qm_r':
                AX.set_xlim([self.min_time*0.95*self.dt, self.max_time*1.05])
#            AX.set_ylabel(AX.get_ylabel, fontsize=27)
        # For last axis
        try:
            lastNonQlkAxis = self.axes[self.non_qlk_params[-1]][1]
            lastNonQlkAxis.set_xlabel("Time (fs)", fontsize=27)
        except IndexError:
            pass
        self.axes[self.plot_params[-1]][1].spines['bottom'].set_visible(True)

        if 'qm_r' not in self.plot_params and len(self.non_qlk_params) > 1:
            self.multi = MultiCursor(self.f.canvas,
                                     [i[1] for i in self.axes.values()],
                                     color='r',
                                     lw=1)

#        if self._use_control:
#            plt.subplots_adjust(top=0.945,
#                                bottom=0.075,
#                                left=0.14,
#                                right=0.99,
#                                hspace=0.084,
#                                wspace=1.00)
        self.f.tight_layout()
        plt.show()
