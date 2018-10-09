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




###############

folder              = '/scratch/mellis/flavoured-cptk/200Rep_3mol'  #'/scratch/mellis/surface_hop/scripts-templates-for-aom-fssh/GENERATOR_FSSH_OS/run-ctmqc-1'
plotting_parameters = ['qm']
replicas            = range(2)

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
        
        self._set_coeff_params()
        
        self._correct_plot_params()
#        self._get_alpha()
        self.run_inp_params = load_inp.get_all_run_inp_variables(self.folder+'run.inp')

        self.title = r"Dimer adiab coeff convergence Ehrenfest -100 reps (without commutator)"
        self.colors = ['r','g','b','#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
        self.colors = [i for j in range(50) for i in self.colors]
        
        self._use_control = True
#        if self.num_reps == 1: self._use_control = False
        
        self.max_time = 100
        
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
            self.plot_params = ['norm', '|c|^2', '|u|^2', 'adiab_states', 'qm', 'norm_traj']
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
                            - 'QM'
                            - 'site_ener'
    
    NOTE: Should be in plot_utils
    """
    def __init__(self, folder, reps, plot_params='all', avg_on=True):
        self.folder = folder
        self.reps = reps
        self.plot_params = plot_params
        self.avg_on = avg_on
        
        self.load_all_ham_data()
        self.load_all_di_coeffs()
        self.load_all_ad_coeffs()
        self.load_ad_ener()
        self.load_qm()
        
        self._average_data()
        
    def load_all_ham_data(self):
        """
        Will load all the hamiltonian data that can be found in the folder 
        specified (dependent on which reps are requested)
        """
        self.all_ham_data = load_ham.load_all_ham_in_folder(self.folder, reps=self.reps)
        self.avg_ham_data = plot_utils.avg_H_data_dict(self.all_ham_data)
        self.avg_site_ener, self.avg_couplings, self.avg_avg_couplings, self.Stimesteps = plot_utils.get_coup_data(self.avg_ham_data, 'avg_ham')
        self.all_site_ener = [plot_utils.get_coup_data(self.all_ham_data, ham_key) for ham_key in self.all_ham_data]
        self.all_site_ener = [[i[0], i[3]] for i in self.all_site_ener]

        self.num_reps = len(self.all_ham_data)
        self._get_alpha()
    
    def load_all_di_coeffs(self):
        """ 
        Loads all the diabatic coefficients, no input. Saves diabatic coeffs as self.all_Dcoeff_data
        """
        if any(j in i for j in ('norm','|u|^2') for i in self.plot_params):
            self.all_Dcoeff_data = load_coeff.load_all_coeff_in_folder(self.folder, filename_must_contain=['xyz','coeff'], filename_must_not_contain=['ad'], reps=self.reps)
        
    def load_all_ad_coeffs(self):     
        """ 
        Loads all the adiabatic coefficients, no input. Saves adiabatic coeffs as self.all_Acoeff_data.
        """
        if '|c|^2' in self.plot_params:
            self.all_Acoeff_data = plot_utils.load_Acoeff_data(self.folder, self.reps, self.all_ham_data)
    
    def load_ad_ener(self):
        """
        Loads the adiabatic energy
        """
        if 'adiab_states' in self.plot_params:
            self.all_ad_ener_data = load_ener.load_all_ener_ad(folder, reps=self.reps)
            if not self.all_ad_ener_data:
                raise IOError("Can't find any data, please check folder.")
    
    def load_qm(self):
        """
        Will load the quantum momentum file into the format in load_QM.
        """
        if 'qm' in self.plot_params:
            self.all_QM_data  = load_QM.load_all_Qlk_in_folder(folder, reps=self.reps)
            self.all_pos_data = load_pos.load_all_pos_in_folder(folder, reps=self.reps)
        
    # Will average the coefficient data (ham is averaged by default)
    def _average_data(self):
        """
        Will average the coefficient data and save as new arrays
        """
        if any(j in i for j in ('norm','|u|^2') for i in self.plot_params) and list(self.all_Dcoeff_data.keys())[0] != 0:
            self.all_Dcoeff_data_avg = plot_utils.avg_coeff_data(self.all_Dcoeff_data)
        if "|c|^2" in self.plot_params and list(self.all_Acoeff_data.keys())[0] != 0:
            self.all_Acoeff_data_avg = plot_utils.avg_coeff_data(self.all_Acoeff_data)
        if 'adiab_states' in self.plot_params:
            self.all_ad_ener_data_avg = plot_utils.avg_E_data_dict(self.all_ad_ener_data)


class Plot(LoadData, Params, plot_norm.Plot_Norm, plot_coeff.Plot_Coeff, 
           plot_ener.Adiab_States, plot_ham.Coupling, plot_QM.QM):
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
        if 'qm' in self.plot_params:
            plot_QM.QM.__init__(self, self.axes['qm'])
        
        #TODO: These need moving into their own filess
        self._plot_site_ener()
#        self._plot_QM()
        
        self.__finalise()
  
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
                
                    # Design the widget axis
                    a[0].set_xticks([])                
                    a[0].set_yticks([])                
                    for side in ['top','bottom','left','right']:
                        a[0].spines[side].set_visible(False)
                else:
                    a.append('')
                    a.append(plt.subplot2grid( (len(self.plot_params),1), (i,0)) )
                self.axes[param] = a
        else:
            plt.close()
            raise SystemExit("Sorry I don't have any way to handle more than 3 plots at the same time yet!")

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
        
        
    #TODO: Move this into it's own file and fix it up so it isn't such a mess!    
    #Will plot the Quantum Momentum term
    def _plot_QM(self):
        """
        Will plot the quantum momentum term
        """
        self.sum_all_atoms_qm = True
        if 'qm' in self.plot_params:
            ax = self.axes['qm'][1]
            if self.plot_all_reps:
#                if self.sum_all_atoms_qm:
#                    for Qlk_filename in self.all_QM_data:
#                        Qlk_data = self.all_QM_data[Qlk_filename]
#                        Qlk_timesteps = Qlk_data[1]
#                        Qlk_data = Qlk_data[0]
#                        
#                        num_atoms = np.shape(Qlk_data[0])[1]/3
#                        QM_mag = load_QM.find_in_Qlk(Qlk_data, params={'at_num':1, 'lk':(1,2), 'cart_dim':1})
#                        QM_mag = np.zeros(np.shape(QM_mag))
#                        
#                        for iatom in range(1,num_atoms+1):
#                            QMX = load_QM.find_in_Qlk(Qlk_data, params={'at_num':iatom, 'lk':(1,2), 'cart_dim':1})
#                            QMY = load_QM.find_in_Qlk(Qlk_data, params={'at_num':iatom, 'lk':(1,2), 'cart_dim':2})
#                            QMZ = load_QM.find_in_Qlk(Qlk_data, params={'at_num':iatom, 'lk':(1,2), 'cart_dim':3})
#                            
#                            QM_mag += QMX**2 + QMY**2 + QMZ**2
#                            
#                        ax.plot(Qlk_timesteps, QM_mag/num_atoms, lw=0.5, alpha=self.alpha)
#                        ylabel = r"$\frac{1}{N_{n}} \sum_{v}$ |Q$_{lk, v}$|$^2$"
##                else:
                    for Qlk_filename in self.all_QM_data:
    #                    Qlk_filename = 'run-QM-1.xyz'
                        Qlk_data = self.all_QM_data[Qlk_filename]
                        Qlk_timesteps = Qlk_data[1]
                        Qlk_data = Qlk_data[0]
                        num_atoms = np.shape(Qlk_data[0])[1]/3
                        for iatom in range(1,num_atoms+1):
    #                        if any(iatom == j for j in (2,4,11,8,1)): continue
                            QMX = load_QM.find_in_Qlk(Qlk_data, params={'at_num':iatom, 'lk':(1,2), 'cart_dim':1})
                            QMY = load_QM.find_in_Qlk(Qlk_data, params={'at_num':iatom, 'lk':(1,2), 'cart_dim':2})
                            QMZ = load_QM.find_in_Qlk(Qlk_data, params={'at_num':iatom, 'lk':(1,2), 'cart_dim':3})
                            
                            QM_mag = QMX**2 + QMY**2 + QMZ**2
                            if len(QM_mag):
                                ax.plot(Qlk_timesteps, QM_mag, label="atom %i"%iatom, lw=0.5, alpha=0.5, color=self.colors[iatom])
                        ylabel = r"|Q$_{lk, v}$|$^2$"
#            if self.avg_on:
#                Qlk_data = self.all_QM_data[self.all_QM_data.keys()[0]][0]
#                QM_mag = load_QM.find_in_Qlk(Qlk_data, params={'at_num':1, 'lk':(1,2), 'cart_dim':1})
#                QM_mag = np.zeros(np.shape(QM_mag))
#                if self.sum_all_atoms_qm:
#                    num_atoms = np.shape(Qlk_data[0])[1]/3
#                    for Qlk_filename in self.all_QM_data:
#                        Qlk_data = self.all_QM_data[Qlk_filename]
#                        Qlk_timesteps = Qlk_data[1]
#                        Qlk_data = Qlk_data[0]
#                        
#                        for iatom in range(1,num_atoms+1):
#                            QMX = load_QM.find_in_Qlk(Qlk_data, params={'at_num':iatom, 'lk':(1,2), 'cart_dim':1})
#                            QMY = load_QM.find_in_Qlk(Qlk_data, params={'at_num':iatom, 'lk':(1,2), 'cart_dim':2})
#                            QMZ = load_QM.find_in_Qlk(Qlk_data, params={'at_num':iatom, 'lk':(1,2), 'cart_dim':3})
#                            
#                            QM_mag += QMX**2 + QMY**2 + QMZ**2
#                            
##                    QM_mag = QM_mag/(self.num_reps*num_atoms)
#                    ax.plot(Qlk_timesteps, QM_mag, lw=1)
#                    ylabel = r"$\frac{1}{N_{n}} \sum_{v}$ |Q$_{lk, v}$|$^2$"
#                else:
#                    pass
            ax.set_ylabel(ylabel)
#            self.title = "Qlk has been normalised for each atom"
#            ax.legend()
    
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
    
    #Will finish off the plots
    def __finalise(self):
        """
        Will finish off the plots by adding necessary (communal) labels etc..
        e.g. will put the time (fs) label on the lowest x axis.
        """
        
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
        self.axes[self.plot_params[-1]][1].set_xlabel("Time (fs)")
        self.axes[self.plot_params[-1]][1].spines['bottom'].set_visible(True)
                        
        self.f.tight_layout()
        self.print_final_info()
#        plt.close()
        plt.show()

#/scratch/mellis/flavoured-cptk/200Rep_3mol
folder = fold.make_fold_abs(folder)


p = Plot(plot_params=plotting_parameters, folder=folder, reps=replicas)


## Will plot many variations of replica number
#for i in range(2,100,1):
#    replicas = range(1,i)    
#    filename = "/homes/mellis/Documents/Graphs/Testing_Ehrenfest/New/Pop_convergence/Dimer/Diff_reps_to_stitch_w_comm/"
#    if p.run_inp_params['FAST_EHRENFEST']:
#        filename += "%i.png"%(p.num_reps)
#    else:
#        filename += "%i_comm.png"%(p.num_reps)
#    p.f.savefig(filename)

