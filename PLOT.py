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

from Plot import plot_utils
from Plot import plot_coeff
from Plot import plot_norm
from Plot import plot_QM
from Plot import plot_ham
from Plot import plot_ener

from IO import Folders as fold

# External Modules
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import datetime
from matplotlib.widgets import MultiCursor
from collections import OrderedDict
import time



###############
#CTMQC_low_coup_2mol
folder              = '/scratch/mellis/flavoured-cptk/CTMQC_low_coup_2mol'  
plotting_parameters = ["|C|^2", 'adiab_states']
replicas            = range(47,49)

###############







class Params(object):
    """
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
#        self._get_alpha()
        self.run_inp_params = load_inp.get_all_run_inp_variables(self.folder+'run.inp')

        self._set_title()
#        self.title = r"Adiab coeffs evolution under **CT/Eh** -**irep** reps (with renormalisation)"

        self.colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
                       '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
                       'r','g','b',]
        self.colors = [i for j in range(50) for i in self.colors]
        
        self._use_control = True
#        if self.num_reps == 1: self._use_control = False
        
        self.max_time = 100
    
    def _set_title(self):
        """
        Will set the title of the plot according to the parameters given.
        """
        params_convert = {'|u|^2':'Diabatic Coefficients', 
                          '|c|^2':'Adiabatic Coefficients',
                          'site_ener':'Site Energies', 
                          'qm_t':"Quantum Momentum", 
                          "adiab_states":"Adiabatic States", 
                          "norm":"Diabatic Norm"}
        if len(self.plot_params) == 1:
            params_joined = self.plot_params[0]
        else:
            params_joined = ', '.join([params_convert[i] for i in self.plot_params[:-1]]) + " and " \
                                + params_convert[self.plot_params[-1]]
            
            
        if self.plot_params[0] != "qm_r":
            self.title = "Evolution of %s under **CT/Eh** with **irep** replicas"%(params_joined)
        else:
            self.title = "Spatial distribution of QM averaged over **irep** replicas"
    
    # Will get the number of replicas and the transparency of the lines.
    def _get_alpha(self):
        """
        Will change the alpha value depending on how many reps were used in the 
        simulation.
        """
        alphas = {1: 1, 2: 0.8, 3: 0.6, 10: 0.4, 50: 0.15, 100: 0.125, 200: 0.1, 500: 0.01, 10000:0}
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
        if self.plot_params == 'all':
            self.plot_params = ['norm', '|c|^2', '|u|^2', 'adiab_states', 'qm_r', 'norm_traj']
        else:    
            self.plot_params = [i.strip().lower() for i in self.plot_params]
        
        self.plot_params = ['coup' if 'coup' in i else i for i in self.plot_params]
        
        self.plot_paramsC = [i for i in self.plot_params if any(i == j for j in ('|u|^2', '|c|^2'))]
        
        # Want to add some typo checking using difflib.SequenceMatcher later

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
                            - 'adiab_states'
                            - 'qm'
                            - 'site_ener'
    
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
        
        self._average_data()
        
        self.print_timing_info(self.load_timings, "Timing Data for Reading Data")
        
    def load_all_ham_data(self):
        """
        Will load all the hamiltonian data that can be foun_set_Qlk_controld in the folder 
        specified (dependent on which reps are requested)
        """
        self.load_timings['H'] = time.time()
        self.all_ham_data = load_ham.load_all_ham_in_folder(self.folder, reps=self.reps)
        self.avg_ham_data = plot_utils.avg_H_data_dict(self.all_ham_data)
        self.avg_site_ener, self.avg_couplings, self.avg_avg_couplings, self.Stimesteps = plot_utils.get_coup_data(self.avg_ham_data, 'avg_ham')
        self.all_site_ener = [plot_utils.get_coup_data(self.all_ham_data, ham_key) for ham_key in self.all_ham_data]
        self.all_site_ener = [[i[0], i[3]] for i in self.all_site_ener]
        self.load_timings['H'] = time.time() - self.load_timings['H']
        

        self.num_reps = len(self.all_ham_data)
        self._get_alpha()
    
    def load_all_di_coeffs(self):
        """ 
        Loads all the diabatic coefficients, no input. Saves diabatic coeffs as self.all_Dcoeff_data
        """
        if any(j in i for j in ('norm','|u|^2') for i in self.plot_params):
            self.load_timings['di coeff'] = time.time()
            self.all_Dcoeff_data = load_coeff.load_all_coeff_in_folder(self.folder, filename_must_contain=['xyz','coeff'], filename_must_not_contain=['ad'], reps=self.reps)
            self.load_timings['di coeff'] = time.time() - self.load_timings['di coeff']
        
    def load_all_ad_coeffs(self):     
        """ 
        Loads all the adiabatic coefficients, no input. Saves adiabatic coeffs as self.all_Acoeff_data.
        """
        if '|c|^2' in self.plot_params:
            self.load_timings['ad coeff'] = time.time()
            self.all_Acoeff_data = plot_utils.load_Acoeff_data(self.folder, self.reps, self.all_ham_data)
            self.load_timings['ad coeff'] = time.time() - self.load_timings['ad coeff']
    
    def load_ad_ener(self):
        """
        Loads the adiabatic energy
        """
        if 'adiab_states' in self.plot_params:
            self.load_timings['adiab ener'] = time.time()
            self.all_ad_ener_data = load_ener.load_all_ener_ad(folder, reps=self.reps)
            if not self.all_ad_ener_data:
                raise IOError("Can't find any data, please check folder.")
            self.load_timings['adiab ener'] = time.time() - self.load_timings['adiab ener']
    
    def load_qm(self):
        """
        Will load the quantum momentum file into the format in load_QM.
        """
        if any('qm' in j  for j in self.plot_params):
            self.load_timings['QM'] = time.time()
            self.all_Qlk_data  = load_QM.load_all_Qlk_in_folder(folder, reps=self.reps)
            self.all_pos_data = load_pos.load_all_pos_in_folder(folder, reps=self.reps)
            self.load_timings['QM'] = time.time() - self.load_timings['QM']
        
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
        
        if "|c|^2" in self.plot_params and list(self.all_Acoeff_data.keys())[0] != 0:
            self.load_timings['Averaging: ']['ad coeff'] = time.time()
            self.all_Acoeff_data_avg = plot_utils.avg_coeff_data(self.all_Acoeff_data)
            self.load_timings['Averaging: ']['ad coeff'] = time.time() - self.load_timings['Averaging: ']['ad coeff']
        
        if 'adiab_states' in self.plot_params:
            self.load_timings['Averaging: ']['adiab_ener'] = time.time()
            self.all_ad_ener_data_avg = plot_utils.avg_E_data_dict(self.all_ad_ener_data)
            self.load_timings['Averaging: ']['adiab_ener'] = time.time() - self.load_timings['Averaging: ']['adiab_ener']
        
        if any('qm' in j for j in self.plot_params):
            self.load_timings['Averaging: ']['qm'] = time.time()
            self.avg_pos_data = plot_utils.avg_pos_data(self.all_pos_data)
            self.avg_Qlk_data = plot_utils.avg_Qlk_data(self.all_Qlk_data)
            self.load_timings['Averaging: ']['qm'] = time.time() - self.load_timings['Averaging: ']['qm']

    def print_timing_info(self, timing_dict, title=""):
        """
        Will print any timings info.
        
        Inputs:
            * timing_dict  =>  the dictionary containing timing data
            * title        =>  the title of the timing data
        """
        max_len = 50
        print(" "+"-"*max_len)
        print("|"+" "*((max_len -len(title)-1)/2) +title+" "*((max_len -len(title))/2) + " |")
        print("|"+" "*max_len+"|")
        plot_utils.print_timings(timing_dict, max_len=max_len)
        print(" "+"-"*max_len)


class Plot(LoadData, Params, plot_norm.Plot_Norm, plot_coeff.Plot_Coeff, 
           plot_ener.Adiab_States, plot_ham.Coupling, plot_QM.QM_R, 
           plot_QM.QM_t):
    """
    Will handle plotting of (hopefully) any parameters. Pass a list of string 
    with the parameters that are to be plotted. E.g. Plot(['|u|^2', '|C|^2']) adiab_states
    and this class should plot them
    
    Inputs:
        plot_params    =>  A list containing the parameters needing plotting. 
                           Possible parameters are:
                               * |u|^2        = Diabatic populations
                               * |C|^2        = Adiabatic populations
                               * norm         = The norm of the diabatic coeffs
                               * qm           = The Quantum Momentum
                               * adiab_states = The adiabatic energy levels
                               * site_ener    = The site energies vs time
        folder         =>  The folder containing the data
        reps           =>  Which replica numbers to plot (can be 'all')
    """
    
    def __init__(self, plot_params, folder, reps):
        self.plot_params = plot_params
        self._correct_plot_params()
        Params.__init__(self, folder, reps, self.plot_params)
        LoadData.__init__(self, folder, reps, self.plot_params)
        self.reps = reps
        self.xlabel = "Time (fs)"
        
        self._create_ax_fig_layout()
        
        ### REMOVE WHEN ALL CLASS PLOTS ARE CREATED ###
        self.plot_all_reps = True
        self.avg_on = True
        self.fill_between = True
        ###############################################
        
        if 'norm' in self.plot_params:
            plot_norm.Plot_Norm.__init__(self, self.axes['norm'])
        if '|u|^2' in self.plot_params:
            plot_coeff.Plot_Coeff.__init__(self, self.axes['|u|^2'])
        if '|c|^2' in self.plot_params:
            plot_coeff.Plot_Coeff.__init__(self, self.axes['|c|^2'])            
        if 'adiab_states' in self.plot_params:
            plot_ener.Adiab_States.__init__(self, self.axes['adiab_states'])
        if 'coup' in self.plot_params:
            plot_ham.Coupling.__init__(self, self.axes['coup'])
        if 'qm_r' in self.plot_params:
            plot_QM.QM_R.__init__(self, self.axes['qm_r'])
        if 'qm_t' in self.plot_params:
            plot_QM.QM_t.__init__(self, self.axes['qm_t'])        
        #TODO: These need moving into their own filess
        self._plot_site_ener()
        
        self.__finalise()
  
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
    
    def _Qlk_axis_special_case(self):
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
        
    #Decides what arrangement of axes to use
    def _create_ax_fig_layout(self):
        """
        Will create the layout for the plots and assign an axis to each plotting
        parameter. This will create the self.f and self.axes variables.
        
        The self.axes variable is a dictionary with the plot parameter as a key
        and the axis that has been assigned to it as the value.
        """
        self.f = plt.figure()
        self.axes = {}
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
        self._Qlk_axis_special_case()
    
    #TODO: Move this into it's own file.
    #Will plot site energy difference.
    def _plot_site_ener(self):
        """
        Will plot the site energy difference on the site_ener axis
        """
        if 'site_ener' in self.plot_params:
            ax = self.axes['site_ener'][1]
            #Plot 1: Site Ener
            ax.set_ylabel(r"$\Delta$E (Ha)")
            if self.avg_on:
                ax.plot(self.Stimesteps, self.avg_site_ener)
        
    def print_final_info(self):
        """
        Will print some information about the data just plotted. This will 
        include stuff from the input file such as whether the commutator was 
        used, whether the run was Ehrenfest or not, the data/time of the plot
        , how many replicas were used and how many molecules.
        """
        
        
        
        # Build the strs list
        str_sections = {'Date/Time':[], 'Ehrenfest':[], 'CTMQC':[]}
        
        strs = str_sections['Date/Time']
        strs.append( datetime.datetime.strftime(datetime.datetime.now(), "Time of plot = %H:%M"))
        strs.append( datetime.datetime.strftime(datetime.datetime.now(), "Date = %d/%m/%Y"))
        
        if not self.run_inp_params['USE_QM']:
            strs = str_sections['Ehrenfest']
            if not self.run_inp_params['FAST_EHRENFEST']: strs.append("Commutator was ON")
            else: strs.append("Commutator was OFF")
            strs.append("Num Replicas = %i"%self.num_reps)
        
        else:
            strs = str_sections['CTMQC']
            if not self.run_inp_params['FAST_EHRENFEST']: strs.append("Commutator was ON")
            else: strs.append("Commutator was OFF")
            strs.append("Num Replicas = %i"%self.num_reps)
            strs.append("Initial Width = %.3g"%self.run_inp_params['INITIAL_SIGMA'])
            
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
    'all'
    def _fill_in_the_title(self):
        """
        Will replace certain words in the title with parameters used in the
        graph etc...
        """
        #Rep num
        self.title = self.title.replace("**irep**", str(self.num_reps))
        
        #Using Comm
        if self.run_inp_params['FAST_EHRENFEST']:
            self.title = self.title.replace("**comm**", "off")
        else: self.title = self.title.replace("**comm**", "on")
        
        #Ehrenfest or CTMQC
        if self.run_inp_params['USE_QM']: 
            self.title = self.title.replace("**CT/Eh**", "CTMQC")
        else: self.title = self.title.replace("**CT/Eh**", "Ehrenfest")
    
    #Will finish off the plots
    def __finalise(self):
        """
        Will finish off the plots by adding necessary (communal) labels etc..
        e.g. will put the time (fs) label on the lowest x axis.
        """
        self._fill_in_the_title()
        if any(j in self.plot_params for j in ('|u|^2', '|c|^2')):
            # Set legend
            if '|u|^2' in self.plot_params: num_states = len(self.all_Dcoeff_data_avg[3][0])
            elif '|c|^2' in self.plot_params: num_states = len(self.all_Acoeff_data_avg[3][0])
            elif 'adiab_states' in self.plot_params: num_states = len(self.state_cols_AS)
            labels = ["State %i"%(i+1) for i in range(num_states)]
            patches = [mpatches.Patch(color=self.colors[i], label=lab) for i, lab in enumerate(labels)]
            self.f.legend(handles=patches, fontsize=20, labels=labels)
        self.f.suptitle(self.title, fontsize=20)
        
        # For all axes
        for ax in self.axes:
            AX = self.axes[ax][1]
            AX.spines['top'].set_visible(False)
            AX.spines['right'].set_visible(False)
            AX.spines['bottom'].set_visible(False)
            AX.spines['left'].set_visible(True)
            AX.grid('on', alpha=0.5)
            
        # For last axis
        self.axes[self.non_qlk_params[-1]][1].set_xlabel("Time (fs)")
        self.axes[self.plot_params[-1]][1].spines['bottom'].set_visible(True)
                        
        self.f.tight_layout()
        self.print_final_info()
        if 'qm_r' not in self.plot_params and len(self.non_qlk_params) > 1:
            self.multi = MultiCursor(self.f.canvas, [i[1] for i in self.axes.values()], color='r', lw=1)
#        plt.close()
        plt.show()

#/scratch/mellis/flavoured-cptk/200Rep_3mol
folder = fold.make_fold_abs(folder)


t1 = time.time()
p = Plot(plot_params=plotting_parameters, folder=folder, reps=replicas)
t2 = time.time()

print("Total time taken = %.0e"%(t2-t1))

## Will plot many variations of replica number
#for i in range(2,100,1):
#    replicas = range(1,i)    
#    filename = "/homes/mellis/Documents/Graphs/Testing_Ehrenfest/New/Pop_convergence/Dimer/Diff_reps_to_stitch_w_comm/"
#    if p.run_inp_params['FAST_EHRENFEST']:
#        filename += "%i.png"%(p.num_reps)
#    else:
#        filename += "%i_comm.png"%(p.num_reps)
#    p.f.savefig(filename)

