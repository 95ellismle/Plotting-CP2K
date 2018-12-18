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
from load import load_tintf
from load import load_frc

from Plot import plot_utils
from Plot import plot_coeff
from Plot import plot_norm
from Plot import plot_QM
from Plot import plot_ham
from Plot import plot_ener
from Plot import plot_frc
from Plot import plot_tintf

from IO import Folders as fold

# External Modules
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import datetime
from matplotlib.widgets import MultiCursor
from collections import OrderedDict
import time
import re
import os
import pandas as pd
from multiprocessing import Pool


###############
#CTMQC_low_coup_2mol
folders = ['',
#           '/scratch/mellis/surface_hop/scripts-templates-for-aom-fssh/GENERATOR_FSSH_OS/TANH_WIDTH_CONV2_2/TANH_WIDTH=0.00/', 
#           '/scratch/mellis/surface_hop/scripts-templatess-for-aom-fssh/GENERATOR_FSSH_OS/TANH_WIDTH_CONV2_1/Scal=0.0003/TANH_WIDTH=0.00/',
#           '/scratch/mellis/surface_hop/scripts-templates-for-aom-fssh/GENERATOR_FSSH_OS/TANH_WIDTH_CONV2_1/Scal=0.03/TANH_WIDTH=0.00/',
           '/scratch/mellis/flavoured-cptk/200Rep_2mol',
#           '/scratch/mellis/flavoured-cptk/200Rep_3mol',
#           '/scratch/mellis/flavoured-cptk/200Rep_3mol_reorder',
#           '/scratch/mellis/flavoured-cptk/200Rep_8mol',
#           '/scratch/mellis/flavoured-cptk/200Rep_8mol_reorder',
#           '/scratch/mellis/flavoured-cptk/200Rep_2mol_same',
#           '/homes/mellis/Documents/Code_bits_and_bobs/Data/SH_data_state_switch',
#           '/scratch/mellis/surface_hop/scripts-templates-for-aom-fssh/GENERATOR_FSSH_OS/run-fssh-0',
           '']

#root_folders = ['',
#                '/scratch/mellis/flavoured-cptk/Tanh_width_data',
#                '']
#folders = []
#for root_folder in root_folders:
#    all_files = os.walk(root_folder)
#    for Dpath, Dnames, Fnames in all_files:
#        if 'run.inp' in Fnames:
#            folders.append(Dpath)



#folders              = ['/scratch/mellis/surface_hop/scripts-templates-for-aom-fssh/GENERATOR_FSSH_OS/TANH_WIDTH_CONV2_1/Scal=0.0003/TANH_WIDTH=8.5e-0',]# '/scratch/mellis/flavoured-cptk/200Rep_2mol_Ehren']
plotting_parameters = ['|u|^2','norm']
replicas            = 'all'
plot                = True
#######################################################





class Params(object):
    """TANH_WIDTH=0.001
    Will store all the parameters needed to plot the graphs. Will also calculate
    parameters such as alpha etc...
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
#        self.title = r"Quantum Momentum and $\sum_{J}|C_1^{ J}|^2|C_2^{ J}|^2(f_1^{J} - f_2^{J})$ for each cartesian dimension for atom 1"
#        self.title = r"Adiab coeffs evolution under **CT/Eh** -**irep** reps (with renormalisation)"

        self.colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
                       '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
                       'r','g','b',]
        self.colors = [i for j in range(50) for i in self.colors]
        
        self._use_control = False
#        if self.num_reps == 1: self._use_control = False
        
        self.max_time      = 'all'     #(in fs)
        self.min_time      = 0       #(in fs)
        self.quick_stride  = 0.01    #(in fs)
        self.slow_stride   = 0.01    #(in fs)
        self.dt = self.run_inp_params['NUCLEAR_TIMESTEP']
        
        self._fix_load_timings()
        
    def _get_inp_data(self):
        """
        Will recursively retrieve all the data from the run.inp file.
        """
        with open(self.folder+'run.inp', 'r') as f:
            inp_file = f.read().split('\n')
        self.nested_inp_params = {}
        self.run_inp_params  = {}
        load_inp.parse_inp_file(self.nested_inp_params, self.run_inp_params, inp_file)
    
    def _fix_load_timings(self):
        """
        Will convert the times to load into step numbers to pass into the 
        load_* functions.
        """
        dt = self.run_inp_params['NUCLEAR_TIMESTEP']
        
        if type(self.max_time) == str:
            if self.max_time == 'all': self.max_step = 'all'
            else: raise SystemExit("""Sorry I don't recognise the setting for 
max_time, you can use all, or specify a maximum time in fs.""")
        else: self.max_step = int(self.max_time/dt)
        self.min_time = int(self.min_time/dt)
        self.quick_stride = int(self.quick_stride/dt)
        self.slow_stride = int(self.slow_stride/dt)
        
        if type(self.max_step) != str and self.max_step < 1: self.max_step = 1 
        if self.min_time < 0: self.min_time = 0
        if self.quick_stride < 1: self.quick_stride = 1 
        if self.slow_stride< 1: self.slow_stride = 1 
        
    def _set_title(self):
        """
        Will set the title of the plot according to the parameters given.
        """
        params_convert = {'|u|^2':'Diab Coeffs', 
                          '|c|^2':'Adiab Coeffs',
                          'qm_t':"Quantum Momentum", 
                          "adiab_state":"Adiabatic States", 
                          "norm":"Norm",
                          'site_ener':'site energy differences',
                          "qm_r":"Quantum Momentum",
                          "fl_fk":"history force state difference", 
                          "fl_fk_CC":"Rlk denominator",
                          'blank':"",
                          "energy_cons":"Energy Conservation",
                          "coup": r"H$_{12}$", 
                          "force": "Nuc. Frc", }
        name_plot_params = self.plot_params[:]
        if 'qm_r' in name_plot_params and 'qm_t' in name_plot_params: 
            name_plot_params.remove("qm_r")
        if len(name_plot_params) == 1:
            params_joined = str(params_convert.get(name_plot_params[0]))
        else:
            params_joined = ', '.join([str(params_convert.get(i)) for i in name_plot_params[:-1]]) + " and " \
                                + str(params_convert.get(name_plot_params[-1]))
        params_joined = params_joined.replace("None", "???")
            
        if self.plot_params[0] != "qm_r":
            self.title = r"Evolution of %s under **CT/Eh** for **irep** replica [H$_{12}$ = **coup** meV]"%(params_joined)
        else:
            self.title = r"Spatial distribution of QM -data from **irep** replica [H$_{12}$ = **coup** meV]"
    
    # Will get the number of replicas and the transparency of the lines.
    def _get_alpha(self):
        """
        Will change the alpha value depending on how many reps were used in the 
        simulation.
        """
        alphas = {1: 1, 2: 0.8, 3: 0.6, 10: 0.4, 50: 0.3, 100: 0.2, 200: 0.15, 500: 0.1, 10000:0.05}
        self.all_alphas = {}
        keys = sorted(alphas.keys())
        for i in range(len(alphas)-1):
            curr_key = keys[i]
            next_key = keys[i+1]
            fit = np.polyfit([curr_key, next_key], [alphas[curr_key], alphas[next_key]],1)
            for i in range(curr_key, next_key):
                self.all_alphas[i] = np.polyval(fit, i)
        self.alpha = self.all_alphas[self.num_reps]
            
    def _correct_plot_params(self):
        """
        Will hopefully correct typos in the plotting parameters.
        """
        # Will handle the 'all' keyword in plot params
#        all_poss_plot_params = ['norm', '|c|^2', '|u|^2', 'adiab_state', 
#                                'qm_r', 'norm_traj', 'site_ener', 'fl_fk',
#                                'fl_fk_CC']
 
        self.plot_params = list(set([i.strip().lower() for i in self.plot_params]))
        
        for short_name in ['adiab_state', 'coup']:
            self.plot_params = [short_name if short_name in i else i for i in self.plot_params]
        
        self.plot_paramsC = [i for i in self.plot_params if any(i == j for j in ('|u|^2', '|c|^2'))]
        
    def _set_coeff_params(self):
        self.check_control_coeff = {}
        self.all_coeff_lines = {}
        self.avg_coeff_lines = {}
        self.coeff_widg_axes = {}
        self.coeff_plot_axes = {}

class LoadData(object):
    """
    Will load all data in the folder. Takes 3 inputs:
        * folder      = folder in which to look for data
        * reps        = which replicas to plot (can be 'all' or list/range of integers)
        * plot_params = which parameters to plot (can be a list of strings)
                        valid parameters:
                            - 'norm'
                            - '|C|^2'
                            - '|u|^2'
                            - 'adiab_state'
                            - 'qm'
                            - 'site_ener'
                            - 'fk_fl' or 'fl_fk'
    
    NOTE: Should be in plot_utils
    """
    def __init__(self, folder, reps, plot_params='all', avg_on=True):
        self.folder = folder
        self.reps = reps
        self.plot_params = plot_params
        self.avg_on = avg_on
        
        self.load_timings = OrderedDict()
        
        self.load_all_ham_data()
        self.load_all_di_coeffs()
        self.load_all_ad_coeffs()
        self.load_ad_ener()
        self.load_qm()
        self.load_hist_f()
        self.load_tot_ener()
        self.load_force()
        
        self._average_data()
        self._get_coupling()
        
        self.print_timing_info(self.load_timings, "Timing Data for Reading Data")
        
    def load_all_ham_data(self):
        """
        Will load all the hamiltonian data that can be foun_set_Qlk_controld in 
        the folder specified (dependent on which reps are requested)
        
        N.B. This is always called as it is normally outputted and it is used 
        to find metadata such as how many steps and reps have been loaded. This
        can be much better though. I will improve it when I have some time!
        """
        self.load_timings['H'] = time.time()
        print_step = self.nested_inp_params['FORCE_EVAL']['MIXED']['ADIABATIC']['PRINT']['HAMILTONIAN']['EACH']['MD'][0]
        if type(self.max_step) == str: max_step = self.max_step
        else: max_step = int(self.max_step/print_step)
        self.all_ham_data = load_ham.load_all_ham_in_folder(self.folder, 
                                                            reps=self.reps, 
                                                            max_step=max_step, 
                                                            min_step=self.min_time, 
                                                            stride=self.quick_stride)
        self.avg_ham_data = plot_utils.avg_H_data_dict(self.all_ham_data)
        self.avg_site_ener, self.avg_couplings, self.avg_avg_couplings, self.Stimesteps = plot_utils.get_coup_data(self.avg_ham_data, 'avg_ham')
        self.all_site_ener = [plot_utils.get_coup_data(self.all_ham_data, ham_key) for ham_key in self.all_ham_data]
        self.all_site_ener = [[i[0], i[3]] for i in self.all_site_ener]
        self.load_timings['H'] = time.time() - self.load_timings['H']
        
        self.num_steps = len(self.avg_ham_data['avg_ham'][0])
        self.num_states = len(self.all_ham_data[list(self.all_ham_data.keys())[0]][0][0])
        self.num_reps = len(self.all_ham_data)
        self._get_alpha()
    
    def load_all_di_coeffs(self):
        """ 
        Loads all the diabatic coefficients, no input. Saves diabatic coeffs as self.all_Dcoeff_data
        """
        if any(j in i for j in ('norm','|u|^2') for i in self.plot_params):
            self.load_timings['di coeff'] = time.time()
            print_step = self.nested_inp_params['FORCE_EVAL']['MIXED']['ADIABATIC']['PRINT']['COEFFICIENTS']['EACH']['MD'][0]
            if type(self.max_step) == str: max_step = self.max_step
            else: max_step = int(self.max_step/print_step)
            self.all_Dcoeff_data = load_coeff.load_all_coeff_in_folder(self.folder, 
                                                                       filename_must_contain=['xyz','coeff'], 
                                                                       filename_must_not_contain=['ad'], 
                                                                       reps=self.reps,
                                                                       max_step=max_step, 
                                                                       min_step=self.min_time, 
                                                                       stride=self.quick_stride)
            if not self.all_Dcoeff_data: raise SystemExit("Sorry I can't find any coeff data in folder:\n\n\t$s"%self.folder)
            self.load_timings['di coeff'] = time.time() - self.load_timings['di coeff']
        
    def load_all_ad_coeffs(self):     
        """ 
        Loads all the adiabatic coefficients, no input. Saves adiabatic coeffs as self.all_Acoeff_data.
        """
        if '|c|^2' in self.plot_params or \
           (any(['fk' in j for j in self.plot_params]) and any(["fl" in j for j in self.plot_params])):
            try:
                print_step = self.nested_inp_params['MOTION']['CTMQC']['PRINT']['AD_COEFF']['EACH']['MD'][0]
            except KeyError:
                print_step = self.nested_inp_params['FORCE_EVAL']['MIXED']['ADIABATIC']['PRINT']['COEFFICIENTS']['EACH']['MD'][0]
            if type(self.max_step) == str: max_step = self.max_step
            else: max_step = int(self.max_step/print_step)
            self.load_timings['ad coeff'] = time.time()
            self.all_Acoeff_data = plot_utils.load_Acoeff_data(self.folder, 
                                                               self.reps, 
                                                               self.all_ham_data,
                                                               max_step=max_step, 
                                                               min_step=self.min_time, 
                                                               stride=self.quick_stride)
            self.load_timings['ad coeff'] = time.time() - self.load_timings['ad coeff']
    
    def load_ad_ener(self):
        """
        Loads the adiabatic energy
        """
        if 'adiab_state' in self.plot_params:
            self.load_timings['adiab ener'] = time.time()
            if type(self.max_step) == str: max_step = self.max_step
            else: max_step = self.max_step*self.dt

            self.all_ad_ener_data = load_ener.load_all_ener_ad(self.folder, 
                                                               reps=self.reps, 
                                                               max_time=max_step,
                                                               min_time=self.min_time * self.dt)
            if not self.all_ad_ener_data:
                raise IOError("Can't find any data, please check folder.")
            self.load_timings['adiab ener'] = time.time() - self.load_timings['adiab ener']
 
    def load_tot_ener(self):
        """
        Loads the total energy files. These files contain the kinetic energy, 
        temperature, potential energy, total energy and CPU time taken.
        """
        if all(j in i for i in self.plot_params for j in ['ener', 'cons']):
            self.load_timings['adiab ener'] = time.time()

            self.all_tot_ener = load_ener.load_all_ener_dat(folder   = self.folder, 
                                                            reps     = self.reps, 
                                                            max_time = self.max_time)

            self.load_timings['adiab ener'] = time.time() - self.load_timings['adiab ener']

    def load_qm(self):
        """
        Will load the quantum momentum file into the format in load_QM.
        """
        if any(['qm' in j  for j in self.plot_params]):
            self.load_timings['QM'] = time.time()
            print_step = self.nested_inp_params['MOTION']['CTMQC']['PRINT']['QM']['EACH']['MD'][0]
            if type(self.max_step) == str: max_step = self.max_step
            else: max_step = int(self.max_step/print_step)
            self.all_Qlk_data  = load_QM.load_all_Qlk_in_folder(self.folder, 
                                                                reps=self.reps,
                                                                max_step=max_step, 
                                                                min_step=self.min_time, 
                                                                stride=self.slow_stride)
            if 'qm_r' in self.plot_params:
                print_step = self.nested_inp_params['MOTION']['PRINT']['TRAJECTORY']['EACH']['MD'][0]
                if type(self.max_step) == str: max_step = self.max_step
                else: max_step = int(self.max_step/print_step)
                self.all_pos_data = load_pos.load_all_pos_in_folder(self.folder, 
                                                                    reps=self.reps,
                                                                    max_step=max_step, 
                                                                    min_step=self.min_time, 
                                                                    stride=self.slow_stride)
            self.load_timings['QM'] = time.time() - self.load_timings['QM']
    
    def load_hist_f(self):
        """
        Will load all the time-integrated history forces in a folder.
        """
        if any(['fk' in j for j in self.plot_params]) and any(["fl" in j for j in self.plot_params]):
            self.load_timings['history forces'] = time.time()
            print_step = self.nested_inp_params['MOTION']['CTMQC']['PRINT']['T_INT_FORCE']['EACH']['MD'][0]
            if type(self.max_step) == str: max_step = self.max_step
            else: max_step = int(self.max_step/print_step)
            self.all_tintf_data = load_tintf.load_all_tintf_in_folder(self.folder, 
                                                                      reps=self.reps,
                                                                      max_step=max_step,
                                                                      min_step=self.min_time,
                                                                      stride  = self.slow_stride)
            self.load_timings['history forces'] = time.time() - self.load_timings['history forces']

    def load_force(self):
        """
        Will load all the nuclear forces.
        """
        if any('force' in j for j in self.plot_params):
            self.load_timings['forces'] = time.time()
            print_step = self.nested_inp_params['MOTION']['PRINT']['FORCES']['EACH']['MD'][0]
            if type(self.max_step) == str: max_step = self.max_step
            else: max_step = int(self.max_step/print_step)
            self.all_frc_data = load_frc.load_all_frc_in_folder(self.folder, 
                                                                  reps=self.reps,
                                                                  max_step=max_step,
                                                                  min_step=self.min_time,
                                                                  stride  = self.slow_stride)
            self.load_timings['forces'] = time.time() - self.load_timings['forces']

    def _get_coupling(self):
        """
        Will get the average coupling and print it to the console
        """
        self.coupling = np.mean(np.abs(self.avg_ham_data['avg_ham'][0][:,0,1]))
        self.coupling *= 27e3 #convert to meV
    
    # Will average the coefficient data (ham is averaged by default)
    def _average_data(self):
        """
        Will average the coefficient data and save as new arrays
        """
        self.load_timings['Averaging: '] = OrderedDict()
        if any(j in i for j in ('norm','|u|^2') for i in self.plot_params) and list(self.all_Dcoeff_data.keys())[0] != 0:
            self.load_timings['Averaging: ']['di coeff'] = time.time()
            self.all_Dcoeff_data_avg = plot_utils.avg_coeff_data(self.all_Dcoeff_data)
            self.load_timings['Averaging: ']['di coeff'] = time.time() - self.load_timings['Averaging: ']['di coeff']
        
        if "|c|^2" in self.plot_params and list(self.all_Acoeff_data.keys())[0] != 0 or \
          (any(['fk' in j for j in self.plot_params]) and any(["fl" in j for j in self.plot_params])):
            self.load_timings['Averaging: ']['ad coeff'] = time.time()
            self.all_Acoeff_data_avg = plot_utils.avg_coeff_data(self.all_Acoeff_data)
            self.load_timings['Averaging: ']['ad coeff'] = time.time() - self.load_timings['Averaging: ']['ad coeff']
        
        if 'adiab_state' in self.plot_params:
            self.load_timings['Averaging: ']['adiab_ener'] = time.time()
            self.all_ad_ener_data_avg = plot_utils.avg_E_data_dict(self.all_ad_ener_data)
            self.load_timings['Averaging: ']['adiab_ener'] = time.time() - self.load_timings['Averaging: ']['adiab_ener']
        
        if any(['qm' in j  for j in self.plot_params]):
            self.load_timings['Averaging: ']['qm'] = time.time()
            if 'qm_r' in self.plot_params:
                self.avg_pos_data = plot_utils.avg_pos_data(self.all_pos_data)
            self.avg_Qlk_data = plot_utils.avg_Qlk_data(self.all_Qlk_data)
            self.load_timings['Averaging: ']['qm'] = time.time() - self.load_timings['Averaging: ']['qm']
        
        if any(['fk' in j for j in self.plot_params]) and any(["fl" in j for j in self.plot_params]):
            self.load_timings['Averaging: ']['history forces'] = time.time()
            self.sum_tintf_data = plot_utils.sum_hist_f_data(self.all_tintf_data)#, self.all_Acoeff_data)
            self.sum_tintf_CC_data = plot_utils.sum_hist_f_CC_data(self.all_tintf_data, self.all_Acoeff_data)
            self.load_timings['Averaging: ']['history forces'] = time.time() - self.load_timings['Averaging: ']['history forces']
        
        if any('force' in j for j in self.plot_params):
            self.load_timings['Averaging: ']['force'] = time.time()
            self.avg_frc_data = plot_utils.avg_pos_data(self.all_frc_data)
            # Rename the average key (using average_pos function)
            self.avg_frc_data['avg_frc'] = self.avg_frc_data['avg_pos']
            del self.avg_frc_data['avg_pos']
            self.load_timings['Averaging: ']['force'] = time.time() - self.load_timings['Averaging: ']['force']

        if all(j in i for i in self.plot_params for j in ['ener', 'cons']):
            self.load_timings['Averaging: ']['Tot. Energy'] = time.time()
            self.tot_ener_mean = pd.concat(self.all_tot_ener.values()).groupby('Step').mean()
            self.load_timings['Averaging: ']['Tot. Energy'] = time.time() - self.load_timings['Averaging: ']['Tot. Energy']

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
        print("|"+" "*int(max_len - num_spaces) + title + " "*int(max_len - num_spaces) + " |")
        print("|"+" "*max_len+"|")
        plot_utils.print_timings(timing_dict, max_len=max_len)
        print(" "+"-"*max_len)


class Plot(LoadData, Params, plot_norm.Plot_Norm, plot_coeff.Plot_Coeff, 
           plot_ener.Adiab_States, plot_ham.Coupling, plot_QM.QM_R, 
           plot_QM.QM_t, plot_ham.Site_Ener, plot_tintf.fl_fk, 
           plot_tintf.fl_fk_CC, plot_ener.Energy_Cons, plot_frc.Plot_Frc):
    """
    Will handle plotting of (hopefully) any parameters. Pass a list of string 
    with the parameters that are to be plotted. E.g. Plot(['|u|^2', '|C|^2']) adiab_state
    and this class should plot them
    
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
    """
    
    def __init__(self, plot_params, folder, reps, plot):
        print("Starting on folder ", folder)
        self.plot_params = plot_params
        self._correct_plot_params()
        Params.__init__(self, folder, reps, self.plot_params)
        LoadData.__init__(self, folder, reps, self.plot_params)
        self.reps = reps
        self.xlabel = "Time (fs)"
        self.plot_info = {}
        self.plot = plot
        
        self._create_ax_fig_layout()
#        
        ### REMOVE WHEN ALL CLASS PLOTS ARE CREATED ###
        self.plot_all_reps = True
        self.avg_on = True
        self.fill_between = True
        ###############################################
        
        self.plot_blank()
        if 'norm' in self.plot_params:
            plot_norm.Plot_Norm.__init__(self, self.axes['norm'])
        if '|u|^2' in self.plot_params:
            plot_coeff.Plot_Coeff.__init__(self, self.axes['|u|^2'])
        if '|c|^2' in self.plot_params:
            plot_coeff.Plot_Coeff.__init__(self, self.axes['|c|^2'])            
        if 'adiab_state' in self.plot_params:
            plot_ener.Adiab_States.__init__(self, self.axes['adiab_state'])
        if 'coup' in self.plot_params:
            plot_ham.Coupling.__init__(self, self.axes['coup'])
        if 'qm_r' in self.plot_params:
            plot_QM.QM_R.__init__(self, self.axes['qm_r'])
        if 'qm_t' in self.plot_params:
            plot_QM.QM_t.__init__(self, self.axes['qm_t'])
        if 'site_ener' in self.plot_params:
            plot_ham.Site_Ener.__init__(self, self.axes['site_ener'])
        if 'fl_fk_cc' in self.plot_params:
            plot_tintf.fl_fk_CC.__init__(self, self.axes['fl_fk_cc'])
        if 'fl_fk' in self.plot_params:
            plot_tintf.fl_fk.__init__(self, self.axes['fl_fk'])
        if all(j in i for i in self.plot_params for j in ['ener', 'cons']):
            plot_ener.Energy_Cons.__init__(self, self.axes['energy_cons'])
        if 'force' in self.plot_params:
            plot_frc.Plot_Frc.__init__(self, self.axes['force'])
        
        if self.plot:
            self.__finalise()
        self.print_final_info()
  
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
        for side in ['top','bottom','left','right']:
            axis.spines[side].set_visible(False)
        
        return axis
    
    def _QM_r_axis_special_case(self):
        """
        A special case layout for the qm_r axis
        """
        if 'qm_r' in self.plot_params:
            self.Qlk_widg_f, ax = plt.subplots(2)
            self.axes['qm_r'] = [0,0] 
            self.axes['qm_r'][0] = ax #[plt.subplot2grid( (len(self.plot_params)*2,7),
                                   #               (self.plot_params.index('qm_r')*2,0), 
                                   #               colspan=1),
                                   #plt.subplot2grid( (len(self.plot_params)*2,7),
                                   #               (self.plot_params.index('qm_r')*2+1,0), 
                                   #               colspan=1)]
            plt.figure(self.f.number)
            self.axes['qm_r'][1] =  plt.subplot2grid( (len(self.plot_params),7),
                                                  (self.plot_params.index('qm_r'),0), 
                                                  colspan=7)
            
            self.axes['qm_r'][0][0] = self._clean_widget_axes(self.axes['qm_r'][0][0])
            self.axes['qm_r'][0][1] = self._clean_widget_axes(self.axes['qm_r'][0][1])
    
    def _N_pane_special_case(self, num_panes, ax_name):
        """
        Will create a special case control panel for axes with 2 panes. The 
        slightly odd ax_name input is to enable multiple substring matches
        in a string. It is a list of substrings that are in the axis name. 
        
        E.g. if ax_name = ['fl', 'fk'] all axes with both fl and fk in their 
             name will be found. If ax_name = ['qm_t'] then all axes with qm_t
             in their name will be found etc...        
        """
        if all([any(name in j for j in self.plot_params) for name in ax_name]):
            ax_name = '_'.join(ax_name).strip('_')
            if ax_name not in self.plot_params: return
            widg_ax = []
            for i in range(num_panes):
                widg_ax.append(plt.subplot2grid( (len(self.plot_params)*num_panes,7),
                                          (self.plot_params.index(ax_name)*num_panes+i,0), 
                                          colspan=1))
            self.axes[ax_name] = [widg_ax,
                           plt.subplot2grid( (len(self.plot_params),7),
                                          (self.plot_params.index(ax_name),1), 
                                          colspan=6)]
            for i in range(num_panes):
                self.axes[ax_name][0][i] = self._clean_widget_axes(self.axes[ax_name][0][i])
    
    #Decides what arrangement of axes to use
    def _create_ax_fig_layout(self):
        """
        Will create the layout for the plots and assign an axis to each plotting
        parameter. This will create the self.f and self.axes variables.
        
        The self.axes variable is a dictionary with the plot parameter as a key
        and the axis that has been assigned to it as the value.
        """
        if self.plot:
            self.f = plt.figure()
            self.axes = OrderedDict()
            if len(self.plot_params) <= 4:   
                for i,param in enumerate(self.plot_params):
                    a = []
                    if self._use_control:
                        a.append( plt.subplot2grid((len(self.plot_params),7),(i,0), colspan=1) )
                        a.append( plt.subplot2grid((len(self.plot_params),7),(i,1), colspan=6) )
                        a[0] = self._clean_widget_axes(a[0])
                    else:
                        a.append('')
                        a.append(plt.subplot2grid( (len(self.plot_params),1), (i,0)) )
                    self.axes[param] = a
            else:
                plt.close()
                raise SystemExit("Sorry I don't have any way to handle more than 3 plots at the same time yet!")
    #        self._QM_r_axis_special_case()
            self._N_pane_special_case(3, ['fl','fk', 'cc'])
            self._N_pane_special_case(3, ['fl','fk'])
            self._N_pane_special_case(3, ['qm_t'])
        else:
            self.axes = {i: '' for i in self.plot_params}
#    
    def plot_blank(self):
        """
        A function to quickly plot on a blank axis. This can be editted by 
        anyone wanting to quickly plot something without making it 
        something that is officialy plotted. Things should be plotted on the 
        'blank' axis.
        """
        
        if 'blank' in self.plot_params:
            widgs, ax = self.axes['blank']
            hdata, _ , timesteps = p.avg_ham_data['avg_ham']
            all_Es = np.array([plot_utils.calc_U_matrix(H)[0] for H in hdata])
            ax.plot(timesteps, all_Es[:,0])
            ax.plot(timesteps, all_Es[:,1])
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
            str_sections = {'Date/Time':[], 'Ehrenfest':[], 'CTMQC':[]}
            for i in self.plot_info:
                str_sections[i] = self.plot_info[i]
            
            # Section 1
            strs = str_sections['Date/Time']
            folder = self.folder.strip('/')[self.folder.strip('/').rfind('/'):].strip('/')
            strs.append("Folder = "+ folder)
            strs.append( datetime.datetime.strftime(datetime.datetime.now(), "Time of plot = %H:%M"))
            strs.append( datetime.datetime.strftime(datetime.datetime.now(), "Date = %d/%m/%Y"))
            strs.append( "dt  =  %.3f"%(self.run_inp_params['NUCLEAR_TIMESTEP']))
            strs.append( "Scaling Factor = %.3g"%(self.run_inp_params['SCALING_FACTOR']))
            strs.append( "Coupling  =  %.3g meV"%self.coupling)
            
            # Section 2
            if not self.run_inp_params['USE_QM']:
                strs = str_sections['Ehrenfest']
                if not self.run_inp_params.get('FAST_EHRENFEST'): strs.append("Commutator was ON")
                else: strs.append("Commutator was OFF")
                strs.append("Num Replicas = %i"%self.num_reps)
            
            else:
                strs = str_sections['CTMQC']
                if not self.run_inp_params.get('FAST_EHRENFEST'): strs.append("Commutator was ON")
                else: strs.append("Commutator was OFF")
                strs.append("Num Replicas = %i"%self.num_reps)
                strs.append("Initial Width = %.3g"%self.run_inp_params['INITIAL_SIGMA'])
                strs.append("Tanh Width = %.2g"%(self.run_inp_params['TANH_WIDTH']))
                
            # Find max len of string and add a hash to the end of each line
            max_len_str = np.max([len(i) for j in str_sections for i in str_sections[j]]) + 4
            str_sections = {sect: [i+ ' '*(max_len_str-len(i))+'#' for i in str_sections[sect]] 
                                   for sect in str_sections}
            
            # Do the printing
            tabN = 6
            tab = " "*tabN
            for i in str_sections:
                if not str_sections[i]: continue
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
            #Rep num
            self.title = self.title.replace("**irep**", str(self.num_reps))
            if self.num_reps > 1: self.title = self.title.replace("replica", "replicas")
            
            #Using Comm
            if self.run_inp_params.get('FAST_EHRENFEST'):
                self.title = self.title.replace("**comm**", "off")
            else: self.title = self.title.replace("**comm**", "on")
            
            #Ehrenfest or CTMQC
            if self.run_inp_params['USE_QM']: 
                self.title = self.title.replace("**CT/Eh**", "CTMQC")
                self.title = self.title.replace("**tanh_width**", "tanh width = "+str(self.run_inp_params['TANH_WIDTH']) )
            else: 
                self.title = self.title.replace("**CT/Eh**", "Ehrenfest")
        else:
            self.title = "%s propagation CP2K"%(self.run_inp_params['METHOD_PROPAGATION'].title().replace("_"," "))
        
        self.title = self.title.replace("**coup**", "%.2g"%self.coupling)
        self.title = re.sub("\*\*.*\*\*", "", self.title)
        self.title = re.sub("\( *\)", "", self.title)

    #Will finish off the plots
    def __finalise(self):
        """
        Will finish off the plots by adding necessary (communal) labels etc..
        e.g. will put the time (fs) label on the lowest x axis.
        """
        self._fill_in_the_title()
        if any([j in self.plot_params for j in ('|u|^2', '|c|^2')]):
            # Set legend
            if '|u|^2' in self.plot_params: num_states = len(self.all_Dcoeff_data_avg[3][0])
            elif '|c|^2' in self.plot_params: num_states = len(self.all_Acoeff_data_avg[3][0])
            elif 'adiab_state' in self.plot_params: num_states = len(self.state_cols_AS)
            elif "site_ener" in self.plot_params: num_states = len(self.avg_ham_data['avg_ham'][0][0])
            labels = ["State %i"%(i+1) for i in range(num_states)]
            patches = [mpatches.Patch(color=self.colors[i], label=lab) for i, lab in enumerate(labels)]
            if num_states < 10:
                fit = np.polyfit([2,10], [18, 10], 1)
                self.f.legend(handles=patches, fontsize=np.polyval(fit, num_states), labels=labels, loc='upper left')
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
                AX.set_xlim([self.min_time*0.95*self.dt,self.max_time*1.05])
#            AX.set_ylabel(AX.get_ylabel, fontsize=27)
        # For last axis
        try:
            self.axes[self.non_qlk_params[-1]][1].set_xlabel("Time (fs)", fontsize=27)
        except IndexError:
            pass
        self.axes[self.plot_params[-1]][1].spines['bottom'].set_visible(True)
                        
        if 'qm_r' not in self.plot_params and len(self.non_qlk_params) > 1:
            self.multi = MultiCursor(self.f.canvas, [i[1] for i in self.axes.values()], color='r', lw=1)
        
#        if self._use_control:
#            plt.subplots_adjust(top=0.945,
#                                bottom=0.075,
#                                left=0.14,
#                                right=0.99,
#                                hspace=0.084,
#                                wspace=1.00)
        self.f.tight_layout()
        plt.show()
        



folders = [fold.make_fold_abs(i) for i in folders if os.path.isdir(i)]
folders = [i for i in folders if os.path.isfile(i+'run.inp')]

if len(folders) > 7 and plot:
    raise SystemExit("It seems you are trying to plot %i images!\n\n\nThat is a lot of figures to load, these would probably be better being saved as png files."%len(folders))

def do_1_folder(folder):
    p = Plot(plot_params=plotting_parameters, folder=folder, reps=replicas, plot=plot)
    return p

if plot:
    all_p = []
    for i, f in enumerate(folders):
        all_p.append(do_1_folder(f))
else:
    num_proc = len(folders)
    if num_proc > 21: num_proc = 21
    print("Num Proc = ", num_proc)
    if num_proc > 1:
        p = Pool(num_proc)
        if __name__ == '__main__':
            all_p = p.map(do_1_folder, folders)
    else:
        all_p = [do_1_folder(folders[0])]

#save_place = '/homes/mellis/Documents/Graphs/CTMQC/Testing_tanh'
#for i, p in enumerate(all_p):
#    filename = saveplaces[i] #save_place+'/Scaling_Factor='+str(p.run_inp_params['SCALING_FACTOR'])+'/All/CTMQC_TW='+str(p.run_inp_params['TANH_WIDTH'])    
#    p.f.savefig(filename+'.png')
    
#SAVE = all_p[:] #just in case I overwrite all_p
#p_by_scal = {}
#for p in all_p:
#    scal = p.run_inp_params['SCALING_FACTOR']
#    if scal not in p_by_scal: p_by_scal[scal] = []
#    p_by_scal[scal].append(p)
#tanh_widths = {scal:[p.run_inp_params['TANH_WIDTH'] for p in p_by_scal[scal]] for scal in p_by_scal}
#norm_drifts = {scal:[p.norm_drift for p in p_by_scal[scal]] for scal in p_by_scal}
