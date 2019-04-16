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

# External Modules
import numpy as np
import datetime
from collections import OrderedDict
import time
import re
import pandas as pd
import os
import difflib

# Decide which things to load when the params are inputted
#  also gives a comprehensive list of all possible plotting
#  parameters
dependencies = {'qm_r':         ['pos', 'qm'],
                'pos':          ['pos'],
                '|c|^2':        ['|c|^2'],
                'qm_t':         ['qm'],
                'rlk':          ['rlk'],
                'site_ener':    ['ham'],
                'coup':         ['ham'],
                '|u|^2':        ['|u|^2'],
                'energy_cons':  ['tot_ener'],
                'energy_drift':  ['tot_ener'],
                'adiab_state':  ['ad_ener'],
                'ylk/sum(ylk)': ['fl_fk', '|c|^2'],
                'fl_fk':        ['fl_fk'],
                'norm':         ['|u|^2'],
                'tot_force':    ['force'],
                'alpha':        ['qm', 'rlk'],
                'pos_sigma':    ['sigma', 'pos'],
                'sum(ylk)':     ['fl_fk', '|c|^2'],
                'k':            ['k'],
                'pos3d':        ['pos'],
                'fl':           ['fl_fk'],
                'fk':           ['fl_fk'],
                "qm_force":     ['qm_frc'], 
                "ad_force":     ['ad_frc'],
                }

class Params(object):
    """
    Will store all the parameters needed to plot the graphs. Will also
    calculate parameters such as alpha etc...
    """
    def __init__(self, folder, reps, plot_params, maxTime='all',
                 minTime=0):
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
        self._use_control = True
#        if self.num_reps == 1: self._use_control = False
        self.max_time = maxTime  # (in fs)
        self.min_time = minTime  # (in fs)  NOT WORKING CAN ONLY USE 0
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
                          "tot_force": "Nuc. Frc", }
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

        allPossParams = np.array(list(dependencies.keys()))
        for iprm, param in enumerate(self.plot_params):
            possParams = [difflib.SequenceMatcher(None,
                                                 param,
                                                 paramReal).ratio()
                          for paramReal in allPossParams]
            possParams = np.array(possParams)
            newParamInd = np.argmax(possParams)
            if (possParams[newParamInd] > 0.7):
                self.plot_params[iprm] = allPossParams[newParamInd]
            else:
               likelyParams = allPossParams[possParams > 0.5]
               msg = "Sorry I don't understand the parameter `%s`\n" % param
               msg += "\nThe likely parameters are:\n\t* %s" % "\n\t* ".join(likelyParams)
               msg += "\n\nAll possible parameters are:\n\t* %s" % "\n\t* ".join(allPossParams)
               raise SystemExit(msg)


                                                 
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
    def __init__(self, folder, reps, plot_params='all', avg_on=True,
                 minTime=0, maxTime='all'):
        if folder[-1] != '/':
            folder += '/'
        self.folder = folder
        self.reps = reps
        self.plot_params = plot_params
        self.avg_on = avg_on
        Params.__init__(self, folder, reps, plot_params,
                        minTime=minTime, maxTime=maxTime)

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
        self.load_qm_frc()
        self.load_ad_frc()

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
            for j in dependencies[i.lower()]:
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
        if type(self.max_step) == str and self.max_step == "all":
            max_step = self.max_step
        else:
            max_step = int(self.max_step/print_step)

        if 'ham' in self.load_params:
            exitCode = self.__load_ham(self.quick_stride, max_step)
            self.all_meta_ham = {}
            for Hkey in self.all_ham_data:
                tmp = plot_utils.get_coup_data(self.all_ham_data, Hkey)
                site_ener, couplings, _, timesteps = tmp
                self.all_meta_ham[Hkey] = {'site_ener_diff': site_ener,
                                           'coup': couplings,
                                           'tsteps': timesteps}
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
                num_steps = len(ltxt) // (lines_in_step)
                stride = int(num_steps / 1000)
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

    def get_qm_type(self):
        """
        Will check if the quantum momentum type saved is the Qlk kind or QM_0.
        """
        qmFiles = [i for i in os.listdir(self.folder) if 'run-QM' in i]
        if 'QM-' in qmFiles[0]:
            self.QM_type = "Qlk"
        else:
            self.QM_type = "QM_0"

    def load_qlk(self, max_step):
        """
        Will load the Quantum Momentum in the Qlk form
        """
        self.all_Qlk_data = load_QM.load_all_Qlk_in_folder(
                                                        self.folder,
                                                        reps=self.reps,
                                                        max_step=max_step,
                                                        min_step=self.min_time,
                                                        stride=self.slow_stride
                                                                   )
        # Find metadata
        Keys = list(self.all_Qlk_data.keys())
        self.num_reps = len(Keys)
        cols = self.all_Qlk_data[Keys[0]][1]
        self.num_qm_steps = len(cols)
        self.num_active_atoms = max(cols[0, :, 0])
        self.num_states = max(cols[0, :, 2])
        
    def load_qm_0(self, max_step):
        """
        Will load the Quantum Momentum in the QM_0 form
        """
        self.all_QM_0_data = load_QM.load_all_QM_0_in_folder(
                                                        self.folder,
                                                        reps=self.reps,
                                                        max_step=max_step,
                                                        min_step=self.min_time,
                                                        stride=self.slow_stride
                                                            )
        # Find metadata
        Keys = list(self.all_QM_0_data.keys())
        self.num_reps = len(Keys)
        cols = self.all_QM_0_data[Keys[0]][0][1]
        self.num_qm_steps = len(cols)
        self.num_active_atoms = cols.shape[0]

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
            self.get_qm_type()
            if self.QM_type == "Qlk":
                self.load_qlk(max_step)
            elif self.QM_type == "QM_0":
                self.load_qm_0(max_step)

            self.load_timings['QM'] = time.time() - self.load_timings['QM']

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
            cols = self.all_pos_data[Keys[0]][1]
            self.num_pos_steps = len(cols)
            self.num_atoms = len(cols[0])

    def load_hist_f(self):
        """
        Will load all the adiab. mom. in a folder.
        """
        if 'fl_fk' in self.load_params:
            self.load_timings['adiab. mom.'] = time.time()
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
            self.load_timings['adiab. mom.'] = time.time() - \
                self.load_timings['adiab. mom.']

            # Find metadata
            Keys = list(self.all_tintf_data.keys())
            self.num_reps = len(Keys)
            tmp = self.all_tintf_data[Keys[0]][0]
            self.num_histf_steps, self.num_states, self.num_active_atoms, _ = tmp.shape

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

    def load_qm_frc(self):
      """
      Will load all the quantum momentum forces.
      """
      if 'qm_frc' in self.load_params:
         self.load_timings['qm_forces'] = time.time()
         print_step = self.nested_inp_params['MOTION']['CTMQC']
         print_step = print_step['PRINT']['QM_FORCE']['EACH']
         print_step = print_step['MD'][0]
         if type(self.max_step) == str:
            max_step = self.max_step
         else:
            max_step = int(self.max_step/print_step)

         self.all_qm_frc_data = load_frc.load_all_qm_frc_in_folder(
                                                     self.folder,
                                                     reps=self.reps,
                                                     max_step=max_step,
                                                     min_step=self.min_time,
                                                     stride=self.slow_stride
                                                            )
         Keys = list(self.all_qm_frc_data.keys())
         self.num_reps = len(Keys)
         #self.num_active_atoms = len(self.all_qm_frc_data[Keys[0]][0][1]) 

         self.load_timings['qm_forces'] = time.time() - \
            self.load_timings['qm_forces']

    def load_ad_frc(self):
      """
      Will load the adiabatic forces
      """
      if 'ad_frc' in self.load_params:
         self.load_timings['adiab_forces'] = time.time()
         print_step = self.nested_inp_params['MOTION']['CTMQC']
         print_step = print_step['PRINT']['AD_FORCES']['EACH']
         print_step = print_step['MD'][0]
         if type(self.max_step) == str:
            max_step = self.max_step
         else:
            max_step = int(self.max_step/print_step)

         self.all_ad_frc_data = load_frc.load_all_ad_frc_in_folder(
                                                     self.folder,
                                                     reps=self.reps,
                                                     max_step=max_step,
                                                     min_step=self.min_time,
                                                     stride=self.slow_stride
                                                            )
         Keys = list(self.all_ad_frc_data.keys())
         self.num_reps = len(Keys)
         tmp = self.all_ad_frc_data[Keys[0]][0].shape
         self.num_steps, self.num_states, self.num_atoms, _ = tmp

         self.load_timings['adiab_forces'] = time.time() - \
            self.load_timings['adiab_forces']


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
            self.load_timings['Averaging: ']['adiab ener'] = time.time()
            self.all_ad_ener_data_avg = \
                plot_utils.avg_E_data_dict(self.all_ad_ener_data)
            self.load_timings['Averaging: ']['adiab ener'] = \
                time.time() - self.load_timings['Averaging: ']['adiab ener']

        if 'pos' in self.load_params:
            self.load_timings['Averaging: ']['pos'] = time.time()
            self.avg_pos_data = plot_utils.avg_pos_data(self.all_pos_data)
            self.load_timings['Averaging: ']['pos'] = time.time() - \
                self.load_timings['Averaging: ']['pos']

        if 'qm' in self.load_params:
            self.load_timings['Averaging: ']['qm'] = time.time()
            if self.QM_type == "Qlk":
                self.avg_Qlk_data = plot_utils.avg_Qlk_data(self.all_Qlk_data)
            elif self.QM_type == "QM_0":
                self.avg_QM0_data = plot_utils.avg_QM0_data(self.all_QM_0_data)
            self.load_timings['Averaging: ']['qm'] = time.time() - \
                self.load_timings['Averaging: ']['qm']

        if 'fl_fk' in self.load_params:
            self.load_timings['Averaging: ']['adiab. mom.'] = time.time()
            self.avg_tintf_data = plot_utils.avg_hist_f_data(
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
            self.load_timings['Averaging: ']['adiab. mom.'] = time.time() -\
                self.load_timings['Averaging: ']['adiab. mom.']

        if 'tot_force' in self.load_params:
            self.load_timings['Averaging: ']['force'] = time.time()
            self.avg_frc_data = plot_utils.avg_pos_data(self.all_frc_data)
            self.load_timings['Averaging: ']['force'] = time.time() - \
                self.load_timings['Averaging: ']['force']

        if 'qm_frc' in self.load_params:
            self.load_timings['Averaging: ']['qm_force'] = time.time()
            self.avg_frc_data = plot_utils.avg_pos_data(self.all_qm_frc_data)
            self.load_timings['Averaging: ']['qm_force'] = time.time() - \
                self.load_timings['Averaging: ']['qm_force']

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

