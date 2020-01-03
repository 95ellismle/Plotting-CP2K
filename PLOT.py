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
from load import load_dlk

from Plot import plot_utils
from Plot import plot_coeff
from Plot import plot_norm
from Plot import plot_QM
from Plot import plot_dlk
from Plot import plot_ham
from Plot import plot_ener
from Plot import plot_frc
from Plot import plot_tintf
from Plot import plot_pos
from Plot import plot_K
from Plot import plot_rabi
from Plot import plot_sig

# External Modules
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import datetime
from matplotlib.widgets import MultiCursor
from mpl_toolkits.mplot3d import Axes3D
from collections import OrderedDict
import time
import re
import pandas as pd
import os
import difflib


# Decide which things to load when the params are inputted
#  also gives a comprehensive list of all possible plotting
#  parameters (must all be lowercase)
dependencies = {'qm_r':             ['pos', 'qm'],
                'pos':              ['pos'],
                'com':              ['pos'],
                'pos_stddev':       ['pos'],
                '|c|^2':            ['|c|^2'],
                'qm_t':             ['qm'],
                'rlk':              ['rlk'],
                'ri0':              ['ri0'],
                'site_ener':        ['ham'],
                'coup':             ['ham'],
                '|u|^2':            ['|u|^2'],
                'energy_cons':      ['tot_ener'],
                'energy_drift':     ['tot_ener'],
                'ylk/sum(ylk)':     ['ad_mom', '|c|^2'],
                'norm':             ['|u|^2$OR$|C|^2'],
                'tot_force':        ['force'],
                'alpha':            ['alpha'],
                'sum(ylk)':         ['ad_mom', '|c|^2'],
                'k':                ['k'],
                'pos3d':            ['pos'],
                'fl':               ['ad_mom'],
                'fk':               ['ad_mom'],
                'coherence':        ['|c|^2'],
                "qm_force":         ['qm_frc'], 
                "rabi":             ['ham'],
                "dlk":              ['dlk'],
                "vel":              ['vel'],
                "sigma":            ['sigma'],
                "ad_mom":           ['ad_mom'],
                "ad_frc":           ['ad_frc'],
                'ad_ener':          ['ad_ener'],
                "h":                ['ham'],
                "eqs26":            ['ad_mom', 'qm', '|c|^2'],
                "admom_diff_adfrc": ['ad_mom', 'ad_frc'],
                }

qmVars = ('ri0', 'rlk', 'qm', 'ad_mom',)

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
        self.units = 'au'

        self._set_title()
        self.title = r""
#        self.title += "   (**tanh_width**)"

        self.colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
                       '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
                       '#ff0000', '#00ff00', '#0000ff']
        self.fs_to_AUt = 41.34137457575099
        self.colors = [i for j in range(50) for i in self.colors]
        self._use_control = True
#        if self.num_reps == 1: self._use_control = False
        self.max_time = maxTime  # (in fs)
        self.min_time = minTime  # (in fs)  NOT WORKING CAN ONLY USE 0
        self.quick_stride = 0  # (in fs)
        self.slow_stride = 0  # (in fs)
        self.dt = self.run_inp_params['NUCLEAR_TIMESTEP']
        self.atoms_to_plot = [1]

        self.worst_reps = {}
        self.best_reps = {}
        self._fix_load_timings()

        print("\n\nREMEMBER TO CHECK UNITS FOR THE THING YOU'RE PLOTTING\n\n")
        self.units = self.units.lower()
        self.allowed_units = ('au', 'cp2k')
        if not any([self.units == j for j in self.allowed_units]):
            msg = "Units string = `%s'" % self.units
            msg += "\nAllowed Units = `[%s]'" % ','.join(self.allowed_units)
            msg += "\n\nI don't know how to handle the unit string."
            msg += "\n\nPlease change it in the PLOT.py file"
            msg += "\n(The unit string tells me how to handle units, options are above)"
            print(msg)
            raise SystemExit("\nERROR: INCORRECT PARAMETER\n\n")

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

        # If we're using the tully model we only have 1 atom
        t_model = self.run_inp_params.get("TULLY_MODEL")
        if t_model is not None and t_model > 0:
           self.num_active_atoms = 1

    def _get_active_atoms(self):
       """
       Will detemine which atoms are active from the AOM_COEFF.include file.
       """
       self.AOM_file = "%sAOM_COEFF.include" % self.folder
       if not os.path.isfile(self.AOM_file):
         print("Can't find the file: `%s`" % self.AOM_file)
       else:
         with open(self.AOM_file, 'r') as f:
            data = np.loadtxt(f, dtype=str, usecols=[-2, -1, 0])
            dataF = np.sum(data[:,[0,1]].astype(float), axis=1)
            dataS = data[:,-1]
            self.electronically_active_atoms = np.arange(0, len(dataF))[dataF != 0]
            self.active_atoms = np.arange(0, len(dataF))[dataS != 'XX']
            if any(j in self.load_params for j in qmVars):
               self.active_atoms = self.electronically_active_atoms

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
        self.quick_stride = int(self.quick_stride/dt)
        self.slow_stride = int(self.slow_stride/dt)

        if type(self.max_step) != str and self.max_step < 1:
            self.max_step = 1

        if self.min_time < 0:
            self.min_time = 0
        self.min_step= int(self.min_time/dt)

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
                          "ad_ener": "Adiabatic State Energy",
                          "norm": "Norm",
                          'site_ener': 'Site Energy Differences',
                          "qm_r": "Quantum Momentum",
                          "ad_mom": "history force state difference",
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

        for short_name in ['ad_ener', 'coup']:
            self.plot_params = [short_name if short_name in i else i
                                for i in self.plot_params]

        allPossParams = np.array(list(dependencies.keys()))
        for iprm, param in enumerate(self.plot_params):
            quickCheck = [i for i in allPossParams if i.lower() == param.lower()]
            if len(quickCheck):
               self.plot_params[iprm] = quickCheck[0]
               continue
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
                            - 'ad_ener'
                            - 'qm'
                            - 'site_ener'
                            - 'fk_fl' or 'ad_mom'
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
        self._get_active_atoms()
        self.load_timings = OrderedDict()

        self.load_all_ham_data()
        self.load_all_di_coeffs()
        self.load_all_ad_coeffs()
        self.load_ad_ener()
        self.load_rlk()
        self.load_RI0()
        self.load_alpha()
        self.load_qm()
        self.load_dlk()
        self.load_pos()
        self.load_vel()
        self.load_adiab_mom()
        self.load_tot_ener()
        self.load_force()
        self.load_sigmas()
        self.load_k()
        self.load_qm_frc()
        self.load_ad_frc()

        self._get_transparency()

        if self.avg_on:
            self._average_data()
        self._calculate_new_data()

        self.print_timing_info(self.load_timings,
                               "Timing Data for Reading Data")

    def decideDependencies(self):
        """
        Will decide which parameters to load from the inputs, i.e. if the |u|^2
        , and norm are requested then it makes no sense to load |u|^2 twice...
        """
        self.load_params = []
        doAtEnd = []
        for i in self.plot_params:
            for j in dependencies[i.lower()]:
                if '$OR$' in j: doAtEnd.append(j); continue
                self.load_params.append(j)

        for i in doAtEnd:
            newDep = i.split("$OR$")
            inDepAlready = [j.lower() in self.load_params for j in newDep]
            if any(inDepAlready): continue
            else: self.load_params.append(newDep[0])

        self.load_params = list(set(self.load_params))

        # Print some feedback
        for i in self.load_params:
            print("Loading %s" % i)

    def __load_ham(self, stride, max_step):
        """
        Will load the hamiltonian (should always use load_all_ham_data instead
        of this).
        """
        if 'ham' in self.load_params:
           self.load_timings['H'] = time.time()
           self.all_ham_data = load_ham.load_all_ham_in_folder(
                                                          self.folder,
                                                          reps=self.reps,
                                                          max_step=max_step,
                                                          min_step=self.min_step,
                                                          stride=stride
                                                              )
           if self.units == 'au':
              for fname in self.all_ham_data:
                  self.all_ham_data[fname][0] /= 27000
                  self.all_ham_data[fname][2] *= self.fs_to_AUt
                  
   
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

        exitCode = False
        if 'ham' in self.load_params:
            exitCode = self.__load_ham(self.quick_stride, max_step)
            self.all_meta_ham = {}
            for Hkey in self.all_ham_data:
                tmp = plot_utils.get_coup_data(self.all_ham_data, Hkey)
                site_ener, couplings, _, timesteps = tmp
                self.all_meta_ham[Hkey] = {'site_ener_diff': site_ener,
                                           'coup': couplings,
                                           'tsteps': timesteps}
        #else:
        #    if type(max_step) == str or \
        #                           (type(max_step) == int and max_step <= 100):
        #        stride = 1
        #    else:
        #        ham_file = [i for i in os.listdir(self.folder)
        #                    if 'hamil' in i][0]
        #        with open(self.folder + ham_file, 'r') as f:
        #            ltxt = f.read().split('\n')
        #        tmp = load_xyz.get_xyz_step_metadata(ltxt, ham_file)
        #        _, _, lines_in_step, num_title_lines = tmp
        #        num_steps = len(ltxt) // (lines_in_step)
        #        stride = int(num_steps / 1000)
        #        if stride < 1:
        #            stride = 1
        #        max_step = num_steps
        #    exitCode = self.__load_ham(stride, max_step)

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

            tmpList = ['xyz', 'coeff']
            self.all_Dcoeff_data = load_coeff.load_all_coeff_in_folder(
                                             self.folder,
                                             filename_must_contain=tmpList,
                                             filename_must_not_contain=['ad'],
                                             reps=self.reps,
                                             max_time=self.max_time,
                                             min_time=self.min_time,
                                             stride=self.quick_stride
                                                                      )
            if self.units == 'au':
               for f in self.all_Dcoeff_data:
                  self.all_Dcoeff_data[f][2] *= 41.34137458

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

            self.load_timings['ad coeff'] = time.time()
            if 'all_ham_data' in self.__dict__.keys():
               self.all_Acoeff_data = plot_utils.load_Acoeff_data(
                                                       self.folder,
                                                       self.reps,
                                                       self.all_ham_data,
                                                       max_time=self.max_time,
                                                       min_time=self.min_time,
                                                       stride=self.quick_stride
                                                              )
            else:
               self.all_Acoeff_data = plot_utils.load_Acoeff_data(
                                                       self.folder,
                                                       self.reps,
                                                       max_time=self.max_time,
                                                       min_time=self.min_time,
                                                       stride=self.quick_stride
                                                              )
            if self.units == 'au':
               for f in self.all_Acoeff_data:
                  self.all_Acoeff_data[f][2] *= 41.34137458
            if not self.all_Acoeff_data:
                raise SystemExit("Sorry I can't find any coeff data in " +
                                 "folder:\n\n\t$s" % self.folder)

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
            if self.units == 'au':
               for f in self.all_ad_ener_data:
                  self.all_ad_ener_data[f]['time'] /= 0.0241884254

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

            if self.units == 'au':
               for f in self.all_tot_ener:
                  self.all_tot_ener[f]['time'] /= 0.0241884254

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

            self.all_sigma = load_sigma.load_all_list_in_folder(
                                                       folder=self.folder,
                                                       reps=self.reps,
                                                       filename_must_contain=['sigma'],
                                                       max_step=max_step,
                                                       min_step=self.min_step,
                                                       stride=self.quick_stride
                                                           )
            if self.units == 'au':
               for fname in self.all_sigma:
                  self.all_sigma[fname][1] *= self.fs_to_AUt
            self.load_timings['sigma'] = time.time() - \
                self.load_timings['sigma']

            # Find metadata
            Keys = list(self.all_sigma.keys())
            self.num_reps = len(Keys)

    def load_alpha(self):
        """
        Will load the alpha data from a simulation
        """
        if 'alpha' in self.load_params:
            self.load_timings['alpha'] = time.time()
            if type(self.max_step) == str:
                max_step = self.max_step
            else:
                max_step = self.max_step*self.dt

            self.all_alpha = load_sigma.load_all_list_in_folder(
                                                       folder=self.folder,
                                                       reps=self.reps,
                                                       filename_must_contain=['alpha'],
                                                       max_step=max_step,
                                                       min_step=self.min_step,
                                                       stride=self.quick_stride
                                                           )
            if self.units == 'au':
               for fname in self.all_alpha:
                  self.all_alpha[fname][1] *= self.fs_to_AUt
            self.load_timings['alpha'] = time.time() - \
                self.load_timings['alpha']

            # Find metadata
            Keys = list(self.all_alpha.keys())
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
                                                     min_step=self.min_step,
                                                     stride=self.quick_stride
                                                    )
            self.load_timings['K'] = time.time() - \
                self.load_timings['K']

            # Find metadata
            Keys = list(self.all_K.keys())
            self.num_reps = len(Keys)


    def get_qm_type(self):
        """
        Will check if the quantum momentum type saved is the Qlk kind or QM_0.
        """
        qmFiles = [i for i in os.listdir(self.folder) if 'run-QM' in i and 'frc' not in i]
        if 'QM-' in qmFiles[0]:
            self.QM_type = "Qlk"
        else:
            self.QM_type = "QM_0"

    def load_dlk(self):
        """
        Will load the Quantum Momentum in the Qlk form
        """
        if 'dlk' in self.load_params:
            self.load_timings['dlk'] = time.time()

            self.all_dlk_data = load_dlk.load_all_dlk_in_folder(
                                                            self.folder,
                                                            reps=self.reps,
                                                            max_time=self.max_time,
                                                            min_time=self.min_time,
                                                            stride=self.slow_stride
                                                                       )
            if self.units == 'au':
               for fname in self.all_dlk_data:
                  self.all_dlk_data[fname]['time'] *= self.fs_to_AUt
            self.load_timings['dlk'] = time.time() - self.load_timings['dlk']

            # Find metadata
            Keys = list(self.all_dlk_data.keys())
            self.num_reps = len(Keys)
            data = self.all_dlk_data[Keys[0]]
            self.num_dlk_steps = max(data['time'])
            self.num_active_atoms = max(data['v'])
            self.num_pair_states = max(data['l'])
        
    def load_qlk(self):
        """
        Will load the Quantum Momentum in the Qlk form
        """
        self.load_timings['Qlk'] = time.time()
        self.all_Qlk_data = load_QM.load_all_Qlk_in_folder(
                                                        self.folder,
                                                        reps=self.reps,
                                                        max_time=self.max_time,
                                                        min_time=self.min_time,
                                                        stride=self.slow_stride
                                                                   )
        if self.units == 'au':
           for fname in self.all_Qlk_data:
              self.all_Qlk_data[fname]['time'] *= self.fs_to_AUt
        # Find metadata
        Keys = list(self.all_Qlk_data.keys())
        self.num_reps = len(Keys)
        data = self.all_Qlk_data[Keys[0]]
        self.num_qm_steps = max(data['time'])
        self.num_active_atoms = max(data['v'])
        self.num_pair_states = max(data['l'])
        self.load_timings['Qlk'] = time.time() - self.load_timings['Qlk']
        
    def load_qm_0(self, max_step):
        """
        Will load the Quantum Momentum in the QM_0 form
        """
        self.all_QM_0_data = load_QM.load_all_QM_0_in_folder(
                                                        self.folder,
                                                        reps=self.reps,
                                                        max_step=max_step,
                                                        min_step=self.min_step,
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
                self.load_qlk()
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
                                                        reps=[1],
                                                        max_time=self.max_time,
                                                        min_time=self.min_time,
                                                        stride=self.slow_stride
                                                               )
            if self.units == 'au':
                self.Rlk_data['time'] *= self.fs_to_AUt
            self.load_timings['Rlk'] = time.time() - self.load_timings['Rlk']

            # Find metadata
            if 'num_reps' not in self.__dict__: self.num_reps = 1
            self.num_rlk_steps = len(set(self.Rlk_data['time']))
            self.num_active_atoms = len(set(self.Rlk_data['v']))
            self.num_pair_states = len(set(self.Rlk_data['l']))

    def load_RI0(self):
        """
        Will load the Rlk file into the format in load_QM.
        """
        if 'ri0' in self.load_params:
            self.load_timings['RI0'] = time.time()
            self.RI0_data = load_QM.load_all_RI0_in_folder(
                                                        self.folder,
                                                        reps=self.reps,
                                                        max_time=self.max_time,
                                                        min_time=self.min_time,
                                                        stride=self.slow_stride
                                                               )
            for f in self.RI0_data:
               if self.units == 'au':
                   self.RI0_data[f][2] *= self.fs_to_AUt
            self.load_timings['RI0'] = time.time() - self.load_timings['RI0']

            # Find metadata
            self.num_reps = len(self.RI0_data)

    def load_vel(self, force_load=False):
        """
        Will load the velocities
        """
        if 'vel' in self.load_params or force_load:
            self.load_timings['velocities'] = time.time()

            self.all_vel_data = load_pos.load_all_vel_in_folder(
                                                       self.folder,
                                                       reps=self.reps,
                                                       max_time=self.max_time,
                                                       min_time=self.min_time,
                                                       stride=self.slow_stride
                                                               )
            self.load_timings['velocities'] = time.time() - \
                self.load_timings['velocities']

            # Find metadata
            Keys = list(self.all_vel_data.keys())
            self.num_reps = len(Keys)
            cols = self.all_vel_data[Keys[0]][1]
            self.num_vel_steps = len(cols)
            self.num_atoms = len(cols[0])

    def load_pos(self):
        """
        Will load the positions
        """
        if 'pos' in self.load_params:
            self.load_timings['positions'] = time.time()

            self.all_pos_data = load_pos.load_all_pos_in_folder(
                                                       self.folder,
                                                       reps=self.reps,
                                                       max_time=self.max_time,
                                                       min_time=self.min_time,
                                                       stride=self.slow_stride
                                                               )
            if self.units == 'au':
               for fname in self.all_pos_data:
                  #self.all_pos_data[fname][0] /= 0.52917720859 
                  self.all_pos_data[fname][2] *= self.fs_to_AUt

            self.load_timings['positions'] = time.time() - \
                self.load_timings['positions']

            # Find metadata
            Keys = list(self.all_pos_data.keys())
            self.num_reps = len(Keys)
            cols = self.all_pos_data[Keys[0]][1]
            self.num_pos_steps = len(cols)
            self.num_atoms = len(cols[0])

    def load_adiab_mom(self):
        """
        Will load all the adiab. mom. in a folder.
        """
        if 'ad_mom' in self.load_params:
            self.load_timings['adiab. mom.'] = time.time()
            print_step = self.nested_inp_params['MOTION']['CTMQC']['PRINT']['AD_MOM']['EACH']['MD'][0]

            if type(self.max_step) == str:
                max_step = self.max_step
            else:
                max_step = int(self.max_step/print_step)

            self.all_adMom_data = load_tintf.load_all_adMom_in_folder(
                                                        self.folder,
                                                        reps=self.reps,
                                                        max_step=max_step,
                                                        min_step=self.min_step,
                                                        stride=self.slow_stride
                                                                     )
            if self.units == 'au':
               for f in self.all_adMom_data:
                  self.all_adMom_data[f]['time'] *= self.fs_to_AUt

            self.load_timings['adiab. mom.'] = time.time() - \
                self.load_timings['adiab. mom.']

            # Find metadata
            Keys = list(self.all_adMom_data.keys())
            self.num_reps = len(Keys)
            self.num_histf_steps = len(self.all_adMom_data[Keys[0]]['time'])

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
                                                        min_step=self.min_step,
                                                        stride=self.slow_stride,
                                                               )
            if self.units == 'au':
               for f in self.all_frc_data:
                  self.all_frc_data[f][2] *= self.fs_to_AUt
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
                                                     min_step=self.min_step,
                                                     stride=self.slow_stride
                                                            )
         if self.units == 'au':
            for f in self.all_qm_frc_data:
               self.all_qm_frc_data[f][2] *= self.fs_to_AUt

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
                                                     max_time=self.max_time,
                                                     min_time=self.min_time,
                                                     stride=self.slow_stride
                                                            )
         if self.units == 'cp2K':
            for key in self.all_ad_frc_data:
               self.all_ad_frc_data[key]['time'] *= self.fs_to_AUt 
         Keys = list(self.all_ad_frc_data.keys())
         self.num_reps = len(Keys)

         self.load_timings['adiab_forces'] = time.time() - \
            self.load_timings['adiab_forces']

    def _calculate_new_data(self):
       """
       Calculate any additional properties.
       """
       if 'eqs26' in self.plot_params:
           self._calc_eqS26()


    def _calc_eqS26(self):
       """
       Calls the neccessary function to calculate equation S26.
       """
       self.load_timings['Calculating: '] = OrderedDict()
       self.load_timings['Calculating: ']['EqS26'] = time.time()

       self.eqS26 = plot_utils.calc_eqS26(self.all_adMom_data, self.all_Acoeff_data, self.all_Qlk_data)
       
       self.load_timings['Calculating: ']['EqS26'] = time.time() - self.load_timings['Calculating: ']['EqS26']
      


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

        if any([j in self.plot_params for j in ('com', 'pos_stddev')]):
            self.load_timings['Averaging: ']['com'] = time.time()
            keys = list(self.all_pos_data.keys())
            mask = self.all_pos_data[keys[0]][1] != 'Ne'
            self.all_COM = plot_utils.calc_COM(self.all_pos_data, mask)
            self.load_timings['Averaging: ']['com'] = time.time() - \
                self.load_timings['Averaging: ']['com']

        if 'dlk' in self.load_params:
            self.load_timings['Averaging: ']['dlk'] = time.time()
            self.avg_dlk_data = plot_utils.avg_Qlk_data(self.all_dlk_data)
            self.load_timings['Averaging: ']['dlk'] = time.time() - \
                self.load_timings['Averaging: ']['dlk']

        if 'qm' in self.load_params:
            self.load_timings['Averaging: ']['qm'] = time.time()
            if self.QM_type == "Qlk":
                self.avg_Qlk_data = plot_utils.avg_Qlk_data(self.all_Qlk_data)
            elif self.QM_type == "QM_0":
                self.avg_QM0_data = plot_utils.avg_QM0_data(self.all_QM_0_data)
            self.load_timings['Averaging: ']['qm'] = time.time() - \
                self.load_timings['Averaging: ']['qm']

        if 'ad_mom' in self.load_params:
            self.load_timings['Averaging: ']['adiab. mom.'] = time.time()
            self.avg_tintf_data = plot_utils.avg_Qlk_data(self.all_adMom_data)
            if 'ylk/sum(ylk)' in self.plot_params:
                self.sum_tintf_CC_data = plot_utils.sum_hist_f_CC_data(
                                                           self.all_adMom_data,
                                                           self.all_Acoeff_data
                                                                      )
            if 'sum(ylk)' in self.plot_params:
                self.sum_ylk = plot_utils.sum_Ylk_data(self.all_adMom_data,
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


class Plot(LoadData, Params, plot_norm.Plot_Norm, plot_coeff.Plot_Coeff,
           plot_ener.Ad_State, plot_ham.Coupling, plot_QM.QM_R,
           plot_QM.Qlk_t, plot_ham.Site_Ener, plot_tintf.AdMom,
           plot_ener.Energy_Cons, plot_frc.Plot_Frc,
           plot_QM.Rlk, plot_QM.Alpha, plot_pos.PlotPos, plot_pos.PosStd,
           plot_K.K, plot_QM.QM0_t, plot_pos.Pos3D,
           plot_pos.COM, plot_rabi.Rabi, plot_coeff.Plot_Deco,
           plot_frc.QM_Frc, plot_frc.Ad_Frc, plot_dlk.dlk,
           plot_pos.PlotVel, plot_sig.Sig, plot_ham.H, plot_QM.RI0): # , plot_ener.Epot):
    """
    Will handle plotting of (hopefully) any parameters. Pass a list of string
    with the parameters that are to be plotted. E.g. Plot(['|u|^2', '|C|^2'])
    ad_ener and this class should plot them

    Inputs:
        plot_params    =>  A list containing the parameters needing plotting.
                           Possible parameters are:
                               * |u|^2        = Diabatic populations
                               * |C|^2        = Adiabatic populations
                               * norm         = The norm of the diabatic coeffs
                               * qm           = The Quantum Momentum
                               * ad_ener = The adiabatic energy levels
                               * site_ener    = The site energies vs time
                               * ad_mom        = The difference in history force
                                                states.
                           more can be found in the dependencies dict at the top of this file.

        folder         =>  The folder containing the data (string)
        reps           =>  Which replica numbers to plot (can be 'all' or a list of int)
        plot           =>  Can be True or False or 'close' to create the graphs
                               but not display.
                             If True will plot to screen,
                                 if False will just analyse data.
                             If 'close' then will not display the graph but
                                 hold it in memory (as self.f).
    """

    def __init__(self, plot_params, folder, reps='all', plot=True,
                 minTime=0, maxTime='all'):
        print("Starting on folder ", folder)
        self.plot_params = plot_params
        self._correct_plot_params()
#        Params.__init__(self, folder, reps, self.plot_params)
        LoadData.__init__(self, folder, reps, self.plot_params,
                          minTime=minTime, maxTime=maxTime)

        if self.atoms_to_plot == 'all':
            self.atoms_to_plot = range(1, self.num_active_atoms+1)

        self.reps = reps
        self.unitsTime = {'cp2k': 'fs', 'au': 'AU'}
        self.unitsLength= {'cp2k': r'bohr', 'au': 'bohr'}
        self.xlabel = "Time (%s)" % self.unitsTime[self.units]
        self.plot_info = {}
        if type(plot) == str:
            if plot.lower() == 'close':
                self.plot = True
            self._create_ax_fig_layout(close=True)
        else:
            self.plot = plot
            self._create_ax_fig_layout(close=False)
        self.plot = plot

        self._plotAllParams()

    def plotEqS26(self, axes):
       """
       Will plot the equation S26 from SI of Min, 17.
       """
       widgAx, plotAx = axes
       plotAx.plot(self.eqS26[0], self.eqS26[1])

    def plotAdFrc_vs_dfdt(self, axes):
       """
       Will plot the derivative of the adiabatic momenta and the adiabatic force.
       These 2 should be exactly on top of each other (in a perfect world), in reality
       they'll be slightly out.
       """
       dt = self.dt * self.fs_to_AUt
       wax, pax = axes
       
       count = 0
       for eF, mF in zip(self.all_ad_frc_data, self.all_adMom_data):
           f = self.all_adMom_data[mF]
           F = self.all_ad_frc_data[eF]
           
           for l in f['l'].unique():
               for v in f['v'].unique():
                   if count == 1:
                     lab1, lab2 = r"$\frac{\delta \mathbf{f}_{l, \nu}^{(I)}}{\delta t}$", r"F$_{ad \ \nu, l}^{ \ (I)}$"
                   else: lab1, lab2 = "", ""

                   actAt = self.active_atoms[v-1] + 1
                   
                   Flv = F[(F['l'] == l) & (F['v'] == actAt)]
                   flv = f[(f['l'] == l) & (f['v'] == v)]
                   Flv.index = Flv['time']
                   flv.index = flv['time']
                   
                   if self.units == 'cp2k':
                      lnD, = pax.plot(flv.index, np.gradient(flv['f(x)'], flv['time']*self.fs_to_AUt), color='r',
                               alpha=0.3, label=lab1)
                   else:
                      lnD, = pax.plot(flv.index, np.gradient(flv['f(x)'], flv['time']), color='r',
                               alpha=0.3, label=lab1)
                   lnA, = pax.plot(Flv.index*self.fs_to_AUt, Flv['Fad(x)'], color='b', alpha=0.3, label=lab2)
                   
                   count += 1
       pax.legend()
       pax.set_ylabel(r"F$_{ad \ \nu, l}^{ \ (I)}$")
                   
               
    def _plotAllParams(self):
        """
        Will plot the necessary graphs by calling the plotting class if the
        correct flag is present.
        """
        self.plot_blank()

        # NACV
        if 'dlk' in self.plot_params:
            self.mDlkPlot = plot_dlk.dlk.__init__(self, self.axes['dlk'])

        # Quantum Momentum
        if 'qm_r' in self.plot_params:
            self.mQMRPlot = plot_QM.QM_R.__init__(self,
                                                  self.axes['qm_r'])
        if 'qm_t' in self.plot_params:
            if self.QM_type == "Qlk":
                self.mQlkTPlot = plot_QM.Qlk_t.__init__(self,
                                                        self.axes['qm_t'])
            elif self.QM_type == "QM_0":
                self.mQM0TPlot = plot_QM.QM0_t.__init__(self,
                                                        self.axes['qm_t'])
        if 'rlk' in self.plot_params:
            self.mRlkPlot = plot_QM.Rlk.__init__(self,
                                                 self.axes['rlk'])
        if 'ri0' in self.plot_params:
            plot_QM.RI0.__init__(self, self.axes['ri0'])

        if 'sigma' in self.plot_params:
            self.mSigPlot = plot_sig.Sig.__init__(self,
                                                  self.axes['sigma'])
        if 'alpha' in self.plot_params:
            self.mAlphaPlot = plot_QM.Alpha.__init__(self, self.axes['alpha'])

        if 'k' in self.plot_params:
            self.mKPlot = plot_K.K.__init__(self, self.axes['k'])

        if 'eqs26' in self.plot_params:
            self.plotEqS26(self.axes['eqs26'])

        # Coefficients
        if 'norm' in self.plot_params:
            self.mNormPlot = plot_norm.Plot_Norm.__init__(self,
                                                          self.axes['norm'])
        if '|u|^2' in self.plot_params:
            self.mDCoeffPlot = plot_coeff.Plot_Coeff.__init__(
                                                             self,
                                                             self.axes['|u|^2'])
        if "rabi" in self.plot_params:
            self.mRabiPlot = plot_rabi.Rabi.__init__(self, 
                                                     self.axes['rabi'])
        if '|c|^2' in self.plot_params:
            self.mACoeffPlot = plot_coeff.Plot_Coeff.__init__(
                                                             self,
                                                             self.axes['|c|^2'])
        if 'coherence' in self.plot_params:
            self.mCoherencePlot = plot_coeff.Plot_Deco.__init__(
                                                             self,
                                                             self.axes['coherence'])
        # Energies
        if 'ad_ener' in self.plot_params:
            self.mAStatesPlot = plot_ener.Ad_State.__init__(
                                                       self,
                                                       self.axes['ad_ener'])
        if 'coup' in self.plot_params:
            self.mCoupPlot = plot_ham.Coupling.__init__(self,
                                                        self.axes['coup'])
        if 'site_ener' in self.plot_params:
            self.mSiteEnerPlot = plot_ham.Site_Ener.__init__(
                                                         self,
                                                         self.axes['site_ener'])
        if 'energy_cons' in self.plot_params:
            self.mEnerConsPlot = plot_ener.Energy_Cons.__init__(
                                                       self,
                                                       self.axes['energy_cons'])
        if 'energy_drift' in self.plot_params:
            self.mEnerConsPlot = plot_ener.Energy_Cons.__init__(
                                                       self,
                                                       self.axes['energy_drift'])
        if 'h' in self.plot_params:
            self.mHPlot = plot_ham.H.__init__(self, self.axes['h'])

        #if 'dynamic_pot_e' in self.plot_params:
        #    plot_ener.Epot.__init__(self)

        # Forces
        if 'ad_mom' in self.plot_params:
            self.mTintfPlot = plot_tintf.AdMom.__init__(self,
                                                         self.axes['ad_mom'])
        if 'tot_force' in self.plot_params:
            self.mForcePlot = plot_frc.Plot_Frc.__init__(
                                                         self,
                                                         self.axes['tot_force'])
        if 'qm_force' in self.plot_params:
            self.mQMForcePlot = plot_frc.QM_Frc.__init__(
                                                         self,
                                                         self.axes['qm_force'])
        if 'ad_frc' in self.plot_params:
            self.mAdForcePlot = plot_frc.Ad_Frc.__init__(
                                                         self,
                                                         self.axes['ad_frc'])
        if 'admom_diff_adfrc' in self.plot_params:
            self.plotAdFrc_vs_dfdt(self.axes['admom_diff_adfrc'])

        # Positions
        if 'pos' in self.plot_params:
            self.mPosPlot = plot_pos.PlotPos.__init__(self,
                                                      self.axes['pos'])
        if 'com' in self.plot_params:
            self.mCOMPlot = plot_pos.COM.__init__(self,
                                                  self.axes['com'])

        if 'pos_stddev' in self.plot_params:
            self.mPosSigPlot = plot_pos.PosStd.__init__(
                                                        self,
                                                        self.axes['pos_stddev'])
        if 'pos3d' in self.plot_params:
            self.mPos3DPlot = plot_pos.Pos3D.__init__(self,
                                                      self.axes['pos3d'])
        if 'vel' in self.plot_params:
            self.mVelPlot = plot_pos.PlotVel.__init__(self,
                                                      self.axes['vel'])

        # Finish up
        if self.plot:
            self.__finalise()
        self.print_final_info()

        if type(self.plot) == str and self.plot == 'close':
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
            self.axes = OrderedDict()
            if 'pos3d' in self.plot_params:
                if len(self.plot_params) != 1:
                    msg = "Sorry I don't support plotting 3D and 2D together"
                    msg += " yet, please don't mix `pos3d` with other inputs"
                    raise SystemExit(msg)
                self.axes['pos3d'] = [self.f.add_subplot(111, projection="3d"),
                                      False]
                return 0
            if close:
                plt.close('all')

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
            str_sections = {'Vitals': [], 'Ehrenfest': [], 'CTMQC': [], 'Best/Worst Reps':[]}
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
            use_qm = self.run_inp_params['USE_QM_COEFF'] \
                * self.run_inp_params['USE_QM_FORCE']
            if not use_qm:
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
                if "TANH_WIDTH" in self.run_inp_params:
                    strs.append("Tanh Width = " +
                                "%.2g" % (self.run_inp_params['TANH_WIDTH']))

            if len(self.best_reps) > 0:
               strs = str_sections['Best/Worst Reps']
               strs.append("  Best Replicas")
               strs.append("Due to        | Rep")
               for qualType in self.best_reps:
                  numSpacesL = 14 - len(qualType)
                  spacesL = " "*numSpacesL
                  strs.append("%s%s|  %s" % (qualType, spacesL, self.best_reps[qualType]))

            if len(self.worst_reps) > 0:
               strs = str_sections['Best/Worst Reps']
               strs.append("     ------------   ")
               strs.append("  Worst Replicas")
               strs.append("Due to        | Rep")
               for qualType in self.worst_reps:
                  numSpacesL = 14 - len(qualType)
                  spacesL = " "*numSpacesL
                  strs.append("%s%s|  %s" % (qualType, spacesL, self.worst_reps[qualType]))

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
            use_qm = self.run_inp_params['USE_QM_COEFF'] \
                * self.run_inp_params['USE_QM_FORCE']
            if use_qm:
                self.title = self.title.replace("**CT/Eh**", "CTMQC")
                if 'TANH_WIDTH' in self.run_inp_params:
                    tmp = str(self.run_inp_params['TANH_WIDTH'])
                    self.title = self.title.replace("**tanh_width**",
                                                    "tanh width = " + tmp)
            else:
                self.title = self.title.replace("**CT/Eh**", "Ehrenfest")
        else:
            tmp = self.run_inp_params['METHOD_PROPAGATION'].title()
            tmp = tmp.replace("_", " ")
            self.title = "%s propagation CP2K" % (tmp)
            self.title = ""

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
            elif 'ad_ener' in self.plot_params:
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
            if AX:
                AX.spines['top'].set_visible(False)
                AX.spines['right'].set_visible(False)
                AX.spines['bottom'].set_visible(False)
                AX.spines['left'].set_visible(True)
                AX.grid('on', alpha=0.5)
                if type(self.max_time) != str and ax != 'qm_r':
                    max_time, min_time = self.max_time, self.min_time
                    if self.units == "au":
                       max_time = self.max_time * self.fs_to_AUt
                       min_time = self.min_time * self.fs_to_AUt
                    range_ = max_time - min_time
                    AX.set_xlim([min_time - (0.05*range_),
                                 max_time + (0.05*range_)])
#            AX.set_ylabel(AX.get_ylabel, fontsize=27)
        # For last axis
        try:
            lastNonQlkAxis = self.axes[self.non_qlk_params[-1]][1]
            if lastNonQlkAxis:
                lastNonQlkAxis.set_xlabel("Time (%s)" % self.unitsTime[self.units],
                                          fontsize=27)
        except IndexError:
            pass
        ax = self.axes[self.plot_params[-1]][1]
        if ax:
            ax.spines['bottom'].set_visible(True)

        if 'qm_r' not in self.plot_params and len(self.non_qlk_params) > 1:
            self.multi = MultiCursor(self.f.canvas,
                                     [i[1] for i in self.axes.values()],
                                     color='r',
                                     lw=0.7)

        #if self._use_control:
        plt.subplots_adjust(top=0.945,
                            bottom=0.09,
                            left=0,
                            right=0.99,
                            hspace=0.3,
                            wspace=0)
        self.f.tight_layout()
        # plt.ion()
        #plt.show()  # block=False)
