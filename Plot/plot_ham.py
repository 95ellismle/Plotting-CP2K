#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 12:23:18 2018

@author: mellis
"""

import matplotlib.pyplot as plt
from matplotlib.widgets import CheckButtons

import numpy as np


class Coupling(object):
    """
    Will plot the couplings in the hamiltonian
    
    Inputs:
        * axes => axes to plot upon. First in list is the control panel axis.
                  Second is the plotting axis [list of plt.axes]
    """
    def __init__(self, axes):
        self.coup_widg_ax, self.coup_plot_ax = axes
        
        #Setting initial default values
        self.avg_reps_coup = True
        self.all_reps_coup = False
        
        #Plotting
        self.all_rep_lines_coup = []
#        self._plot_all_rep_coup()
        self._plot_avg_coup()
        
        #Connect checkboxes to plot control
        if self._use_control:    self._set_coup_control()
        
        self.coup_plot_ax.set_ylabel(r"$H_{ab}$")
    
    def _check_settings_coup(self, label):
        if label == 'all replicas': #Pressed the all replicas button
            for line in self.all_rep_lines_coup:
                line.set_visible(not line.get_visible())
                
        elif label == 'average':
            self.avg_line_coup.set_visible(not self.avg_line_coup.get_visible())
        plt.draw()
            
    #Will set the control panel for coup graph
    def _set_coup_control(self):
        self.check_coup = CheckButtons(self.coup_widg_ax, ('all replicas', 'average'), (self.all_reps_coup, self.avg_reps_coup))
        self.check_coup.on_clicked(self._check_settings_coup)
        
#    #Will plot the all replica couplings
#    def _plot_all_rep_coup(self):
#        self.all_rep_lines_coup = []
#        for Dfilename in self.all_Dcoeff_data:
#            coeffs, cols, timesteps, pops = self.all_Dcoeff_data[Dfilename]
#            coups = np.sum(pops, axis=1)
#            self.all_rep_lines_coup.append(self.plot_ax.plot(timesteps, coups, alpha=self.alpha, color='r', lw=1)[0])
#        
#        #Initialise the replica lines
#        for line in self.all_rep_lines_coup:
#            line.set_visible(self.all_reps_coup)


    #Will plot coup of diabatic coeffs
    def _plot_avg_coup(self):
        """
        Will plot the off-diagonal hamiltonian elements.
        """
        avg_Hs,timesteps = self.avg_ham_data['avg_ham']
#        avg_coups = np.tril(avg_Hs[0])
        avg_coups = avg_Hs[:,0,1]
        self.coup_plot_ax.plot(timesteps, avg_coups)        