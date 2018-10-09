#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  3 14:35:47 2018

@author: mellis
"""
from matplotlib.widgets import CheckButtons
import matplotlib.pyplot as plt


class Plot_Coeff(object):
    """
    Will plot the coefficient data. 
        
    Inputs:
        axes  => a list of the axes to plot on. The first item should be the 
                 widget axis the second will be the axis to plot the data.  
    
    NOTE: This really is a bit of a hack, I should really go back and have a 
          think of a more elegant to code this later.
    """
    def __init__(self, axes):
        self._set_di_or_ad()

        self.coeff_widg_axes[self.di_or_ad] = axes[0]        
        self.coeff_plot_axes[self.di_or_ad] = axes[1]        
        
        self.all_reps_coeff = False
        self.avg_reps_coeff = True
        
        if self._use_control:
            self._set_coeff_control()
        self._set_coeff_data()
        
        self._plot_all_reps_coeff()
        self._plot_avg_reps_coeff()
        
        self.coeff_plot_axes[self.di_or_ad].set_ylabel(self.ylabel)
        
    def _set_coeff_data(self):
        """
        Will set which data is being used (adiabatic or diabatic)
        """
        if self.di_or_ad == "|u|^2":
            self.coeff_data = self.all_Dcoeff_data
            self.all_coeff_data_avg = self.all_Dcoeff_data_avg
            self.num_states = len(self.all_Dcoeff_data_avg[0][0])
            self.ylabel     = r"$|u_l|^2$"
            
        elif self.di_or_ad == "|c|^2":
            self.coeff_data = self.all_Acoeff_data
            self.all_coeff_data_avg = self.all_Acoeff_data_avg
            self.num_states = len(self.all_Acoeff_data_avg[0][0])
            self.ylabel     = r"$|C_l|^2$"
    
    def _set_di_or_ad(self):
        """
        Will determine whether we are working with adiabatic data or diabatic data
        """
        if '|u|^2' in self.plot_paramsC: 
            self.di_or_ad = '|u|^2'
            self.plot_paramsC.remove("|u|^2")
        elif '|c|^2' in self.plot_paramsC:
            self.di_or_ad = '|c|^2'
        
        
    def _plot_all_reps_coeff(self):
        """
        Will plot the graph showing all the replicas
        """
        self.all_coeff_lines[self.di_or_ad] = []
        for filename in self.coeff_data:
            coeffs, cols, timesteps, pops = self.coeff_data[filename]
            for i in range(self.num_states):
                self.all_coeff_lines[self.di_or_ad].append(self.coeff_plot_axes[self.di_or_ad].plot(timesteps, pops[:,i], color=self.colors[i], alpha=self.alpha, lw=0.7)[0])
        
        # Set initial visibility
        for line in self.all_coeff_lines[self.di_or_ad]:
            line.set_visible(self.all_reps_coeff)

    def _plot_avg_reps_coeff(self):
        """
        Will plot the average coefficient lines for each state.
        """
        self.avg_coeff_lines[self.di_or_ad] = []
        coeffs, cols, timesteps, pops = self.all_coeff_data_avg
        for i in range(self.num_states):                    
            self.avg_coeff_lines[self.di_or_ad].append(self.coeff_plot_axes[self.di_or_ad].plot(timesteps, pops[:,i], '.', color=self.colors[i])[0])
        
        # Set initial visibility
        for line in self.avg_coeff_lines[self.di_or_ad]:
            line.set_visible(self.avg_reps_coeff)
                
    def _set_coeff_control(self):
        """
        Will set the control panel for the coeff graphs
        """
        if self.di_or_ad == '|u|^2':
            self.check_control_coeff[self.di_or_ad] = CheckButtons(self.coeff_widg_axes[self.di_or_ad], ('all Dreplicas', 'Daverage'), (self.all_reps_coeff, self.avg_reps_coeff))
            self.check_control_coeff[self.di_or_ad].on_clicked(self._check_settings_coeff)
        
        elif self.di_or_ad == '|c|^2':
            self.check_control_coeff[self.di_or_ad] = CheckButtons(self.coeff_widg_axes[self.di_or_ad], ('all replicas', 'average'), (self.all_reps_coeff, self.avg_reps_coeff))
            self.check_control_coeff[self.di_or_ad].on_clicked(self._check_settings_coeff)

    def _check_settings_coeff(self, label):
        """
        Will handle the switching on and off for the coeff graphs.
        """
        if label == 'all replicas': #Pressed the all replicas button
            for line in self.all_coeff_lines['|c|^2']:
                line.set_visible(not line.get_visible())
                
        elif label == 'average':
            for line in self.avg_coeff_lines['|c|^2']:
                line.set_visible(not line.get_visible())
        
        elif label == 'all Dreplicas': #Pressed the all replicas button
            for line in self.all_coeff_lines['|u|^2']:
                line.set_visible(not line.get_visible())
                
        elif label == 'Daverage':
            for line in self.avg_coeff_lines['|u|^2']:
                line.set_visible(not line.get_visible())
        plt.draw()     