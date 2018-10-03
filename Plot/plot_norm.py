#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  3 15:53:35 2018

@author: mellis
"""
from matplotlib.widgets import CheckButtons
import matplotlib.pyplot as plt
import numpy as np

class Plot_Norm(object):
    """
    Will plot the normalisation graph.
    
    Inputs:
        axes  => a list of the axes to plot on. The first item should be the 
                 widget axis the second will be the axis to plot the data.
    """
    def __init__(self, axes):
        self.widget_ax = axes[0]
        self.plot_ax = axes[1]
        
        #Setting initial default values
        self.avg_reps_norm = True
        self.all_reps_norm = False
        
        #Plotting
        self.all_rep_lines_norm = []
        self._plot_all_rep_norm()
        self._plot_norm_graph()
        
        #Connect checkboxes to plot control
        if self._use_control:    self._set_norm_control()
        
        self.plot_ax.set_ylabel(r"$\sum_k |u_k^{I}|^2$")
    
    def _check_settings_norm(self, label):
        if label == 'all replicas': #Pressed the all replicas button
            for line in self.all_rep_lines_norm:
                line.set_visible(not line.get_visible())
                
        elif label == 'average':
            self.avg_line_norm.set_visible(not self.avg_line_norm.get_visible())
        plt.draw()
            
    #Will set the control panel for norm graph
    def _set_norm_control(self):
        self.check_norm = CheckButtons(self.widget_ax, ('all replicas', 'average'), (self.all_reps_norm, self.avg_reps_norm))
        self.check_norm.on_clicked(self._check_settings_norm)
        
    #Will plot the all replica norms
    def _plot_all_rep_norm(self):
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
        coeffs, cols, timesteps, pops = self.all_Dcoeff_data_avg
        norms = np.sum(pops, axis=1)
        self.avg_line_norm, = self.plot_ax.plot(timesteps, norms, lw=2, color='g')
        
        self.avg_line_norm.set_visible(self.avg_reps_norm)