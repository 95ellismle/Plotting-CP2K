#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 16:47:13 2018

@author: mellis
"""
from load import load_coeff
from load import load_ener
from load import load_ham
from load import load_QM

from Plot import plot_utils

from IO import Folders as fold

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.widgets import CheckButtons
import matplotlib.cm as cm
import numpy as np
import os


class Params(object):
    """
    Will store all the parameters needed to plot the graphs. Will also calculate
    parameters such as alpha etc...
    """
    
    def __init__(self, folder, reps, plot_params):
        self.folder = folder
        self.reps = reps
        self.plot_params = plot_params
        
        self._correct_plot_params()
        self._get_alpha()

        self.title = ""
        self.colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
        

#        if self.plot_all_reps: self.fill_between = False
        
        
    # Will get the number of replicas and the transparency of the lines.
    def _get_alpha(self):
        """
        Will get the number of replicas in the folder. Will just count how many 
        position files there are. Will also set the transparency of the replica
        lines.
        """
        if self.folder:
            if self.reps == 'all':    self.num_reps = len([i for i in os.listdir(self.folder) if 'n-pos-' in i])
            else:                     self.num_reps = len(self.reps)
            if self.num_reps == 0: raise SystemExit("Sorry I can't seem to find any replicas! (self.num_reps = %i)"%self.num_reps)
            
            alphas = {1: 1, 2: 0.5, 3: 0.4, 10: 0.2, 50: 0.15, 100: 0.125, 200: 0.1, 500: 0.01, 10000:0}
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
        # Want to add some typo checking using difflib.SequenceMatcher later

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
            self.all_QM_data = load_QM.load_all_Qlk_in_folder(folder, reps=self.reps)
        
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


class Plot_Norm(object):
    """
    Will plot the normalisation graph.
    """
    def __init__(self, axes):
        self.widget_ax = axes[0]
        self.plot_ax = axes[1]
        
        self.all_reps_norm = False
        self.avg_reps_norm = True
        
        self.all_rep_lines_norm = []
        self._plot_all_rep_norm()
        self._plot_norm_graph()

        self._set_norm_control()
        
        self.plot_ax.set_ylabel(r"$\sum_k |u_k^{I}|^2$")
    
    def _check_settings(self, label):
        if label == 'all replicas': #Pressed the all replicas button
            for line in self.all_rep_lines_norm:
                line.set_visible(self.all_reps_norm)
            self.all_reps_norm = not self.all_reps_norm
                
        elif label == 'average':
            self.avg_reps_norm = not self.avg_reps_norm
            self.avg_line_norm.set_visible(self.avg_reps_norm)
            
    #Will set the control panel for norm graph
    def _set_norm_control(self):
        self.check = CheckButtons(self.widget_ax, ('all replicas', 'average'), (self.all_reps_norm, self.avg_reps_norm))
        self.check.on_clicked(self._check_settings)
        
    #Will plot the all replica norms
    def _plot_all_rep_norm(self):
        self.all_reps_norm = True
        self.all_rep_lines_norm = []
        for Dfilename in self.all_Dcoeff_data:
            coeffs, cols, timesteps, pops = self.all_Dcoeff_data[Dfilename]
            norms = np.sum(pops, axis=1)
            self.all_rep_lines_norm.append(self.plot_ax.plot(timesteps, norms, alpha=self.alpha, color='r', lw=1)[0])
        
        #Initialise the replica lines
        for line in self.all_rep_lines_norm:
            line.set_visible(self.all_reps_norm)

    #Will plot normal of diabatic coeffs
    def _plot_norm_graph(self):
        """
        Will plot the normal of the diabatic coefficients. Will sum the 
        populations and plot on the norm axis.
        """
        self.avg_reps_norm = True
        coeffs, cols, timesteps, pops = self.all_Dcoeff_data_avg
        norms = np.sum(pops, axis=1)
        self.avg_line_norm, = self.plot_ax.plot(timesteps, norms, lw=2, color='g')
        self.avg_line_norm.set_visible(self.avg_reps_norm)


class Plot(LoadData, Params, Plot_Norm):
    """
    Will handle plotting of (hopefully) any parameters. Pass a list of string 
    with the parameters that are to be plotted. E.g. Plot(['|u|^2', '|C|^2']) adiab_states
    and this class should plot them
    
    Inputs:
        plot_params    =>  a list containing the parameters needing plotting.    
    """
    
    def __init__(self, plot_params, folder, reps):
        self.plot_params = plot_params
        self._correct_plot_params()
        Params.__init__(self, folder, reps, self.plot_params)
        LoadData.__init__(self, folder, reps, self.plot_params)
        self.reps = reps
        
        self._create_ax_fig_layout()
        
        if 'norm' in self.plot_params:
            Plot_Norm.__init__(self, self.axes['norm'])
        self._plot_site_ener()
        self._plot_diab_pops()
        self._plot_adiab_pops()
        self._plot_adiab_states()
        self._plot_QM()
        
        self.__finalise()
  
    #Decides what arrangement of axes to use
    def _create_ax_fig_layout(self):
        """
        Will create the layout for the plots and assign an axis to each plotting
        parameter. This will create the self.f and self.axes variables.
        
        The self.axes variable is a dictionary with the plot parameter as a key
        and the axis that has been assigned to it as the value.
        """
        self.f = plt.figure(figsize=(50,50))
        self.axes = {}
        if len(self.plot_params) <= 3:   
            for i,param in enumerate(self.plot_params):
                a = []
                a.append( plt.subplot2grid((len(self.plot_params),9),(i,0), colspan=1) )
                a.append( plt.subplot2grid((len(self.plot_params),9),(i,1), colspan=9) )
                a[0].set_xticks([])                
                a[0].set_yticks([])                
                self.axes[param] = a
        else:
            plt.close()
            raise SystemExit("Sorry I don't have any way to handle more than 3 plots at the same time yet!")
                
#            ax1 = plt.subplot2grid((4, 2), (2, 2))
#            ax2 = plt.subplot2grid((4, 2), (2, 1), colspan=1)
#            ax3 = plt.subplot2grid((3, 3), (1, 0), colspan=2, rowspan=2)
#            ax4 = plt.subplot2grid((3, 3), (1, 2), rowspan=2)

#        elif len(self.plot_params) <= 3:        self.f,a = plt.subplots(len(self.plot_params), figsize=figsize)
#        elif len(self.plot_params) == 4:           
#            self.f,a = plt.subplots(2,2, figsize=figsize)
#            a = a.flatten()
#        else: raise SystemExit("Sorry no rule has been set for a arrangement of more than 4 axes!")
#        self.axes = {param:a[i] for i,param in enumerate(self.plot_params)}
    


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
            
    #Will plot diab pops
    def _plot_diab_pops(self):
        """
        Will plot the diabatic populations on the |u|^2 axis.
        """
        if '|u|^2' in self.plot_params:
            ax = self.axes['|u|^2'][1]
            num_states = len(self.all_Dcoeff_data_avg[3][0])
            
            #Plot averages
#            if self.avg_on:
            coeffs, cols, timesteps, pops = self.all_Dcoeff_data_avg
            for i in range(num_states):                    
                ax.plot(timesteps, pops[:,i], color=self.colors[i])
                    
            #Plot faded reps
#            if self.plot_all_reps:
            for Dfilename in self.all_Dcoeff_data:
                coeffs, cols, timesteps, pops = self.all_Dcoeff_data[Dfilename]
                for i in range(num_states):
                    ax.plot(timesteps, pops[:,i], color=self.colors[i], alpha=self.alpha, lw=0.7)
            
            
            ax.set_ylabel(r"$|u_l|^2$")
        
    #Will plot adiab pops
    def _plot_adiab_pops(self):
        """
        Will plot the adiabatic populations on the |C|^2 axis.
        """
        if '|c|^2' in self.plot_params:
            ax = self.axes['|c|^2'][1]
#            if self.avg_on:
            coeffs, cols, timesteps, pops = self.all_Acoeff_data_avg
            for i in range(len(pops[0])):
                ax.plot(timesteps, pops[:,i], color=self.colors[i])
#            if self.plot_all_reps:
            for Afilename in self.all_Acoeff_data:
                coeffs, cols, timesteps, pops = self.all_Acoeff_data[Afilename]
                for i in range(len(pops[0])):
                    ax.plot(timesteps, pops[:,i], color=self.colors[i], alpha=self.alpha, lw=1)
            ax.set_ylabel(r"$|C_l|^2$")
    
    #Will plot the adiabatic states
    def _plot_adiab_states(self):
        """
        Will plot the adiabatic energy levels.
        """
        if 'adiab_states' in self.plot_params:
            ax = self.axes['adiab_states'][1]
            state_cols = [i for i in self.all_ad_ener_data_avg.columns if 'State' in i]
            
            if self.avg_on:
                for i, col in enumerate(state_cols):
                    ax.plot(self.all_ad_ener_data_avg['Time'], self.all_ad_ener_data_avg[col], color=self.colors[i])
            
            if self.plot_all_reps:
                for Efilename in self.all_ad_ener_data:
                    ad_ener_data = self.all_ad_ener_data[Efilename]
                    for i, col in enumerate(state_cols):
                        ax.plot(ad_ener_data['Time'], ad_ener_data[col], alpha=self.alpha, lw=0.7, color=self.colors[i])
            
            if self.fill_between:
                for i, state in enumerate(state_cols[:-1]):
                    y1_y2 = self.all_ad_ener_data_avg[state] - self.all_ad_ener_data_avg[state_cols[i+1]]
                    max_diff, min_diff = np.max(y1_y2), np.min(y1_y2)
                    normed_diffs = (y1_y2 - min_diff)/(max_diff - min_diff)
                    for dt, timestep in enumerate(self.all_ad_ener_data_avg['Time'][:-1]):
                        T = [timestep, self.all_ad_ener_data_avg['Time'][dt+1]]
                        y1 = self.all_ad_ener_data_avg[state].iloc[[dt, dt+1]]
                        y2 = self.all_ad_ener_data_avg[state_cols[i+1]].iloc[[dt, dt+1]]
                        diffy = 1-normed_diffs[dt]
                        color = cm.hot(diffy*2)
                        ax.fill_between(T, y1, y2, alpha=0.1, color=color, lw=0)

            ax.set_ylabel(r"$E^{ad}_{l}$")
    
    #Will plot the Quantum Momentum term
    def _plot_QM(self):
        """
        Will plot the quantum momentum term
        """
        if 'qm' in self.plot_params:
            ax = self.axes['qm'][1]
            if self.plot_all_reps:
                for Qlk_filename in self.all_QM_data:
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
                        QM_mag /= np.max(QM_mag)
                        ax.plot(Qlk_timesteps, QM_mag, label="atom %i"%iatom)
            ax.set_ylabel(r"|Q$_{lk}$|$^2$")
#            ax.legend()
            
    #Will finish off the plots
    def __finalise(self):
        """
        Will finish off the plots by adding necessary (communal) labels etc..
        e.g. will put the time (fs) label on the lowest x axis.
        """
        
        if any(j in self.plot_params for j in ('|u|^2', '|c|^2')):
            # Set legend
            if '|u|^2' in self.plot_params: num_states = len(self.all_Dcoeff_data_avg[3][0])
            if '|c|^2' in self.plot_params: num_states = len(self.all_Acoeff_data_avg[3][0])
            leg_dict = {'State %i'%i:self.colors[i] for i in range(num_states)}
            patches = [mpatches.Patch(color=leg_dict[i], label=i) for i in leg_dict]
            self.f.legend(handles=patches, fontsize=20, labels=leg_dict.keys())
        self.f.suptitle(self.title, fontsize=20)
        
        # For all axes
        for ax in self.axes:
            self.axes[ax][1].spines['top'].set_visible(False)
            self.axes[ax][1].spines['right'].set_visible(False)
            self.axes[ax][1].spines['bottom'].set_visible(False)
            self.axes[ax][1].spines['left'].set_visible(True)
            self.axes[ax][1].grid('on', alpha=0.5)

        # For last axis
        self.axes[self.plot_params[-1]][1].set_xlabel("Time (fs)")
        self.axes[self.plot_params[-1]][1].spines['bottom'].set_visible(True)
        
        self.f.tight_layout()
        
folder = fold.make_fold_abs('../Data/200Rep_3mol') #The folder to look in for the data
p = Plot(['|u|^2','norm', '|c|^2'], folder, 'all')
#plt.subplots_adjust(hspace=0.1)
