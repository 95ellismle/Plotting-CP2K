#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  3 15:53:35 2018

@author: mellis
"""
from matplotlib.widgets import CheckButtons
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

def linear_fit(x, m ,c):
    return m*x + c

class Plot_Norm(object):
    """
    Will plot the normalisation graph.
    
    Inputs:
        axes  => a list of the axes to plot on. The first item should be the 
                 widget axis the second will be the axis to plot the data.
    """
    def __init__(self, axes):
        self.plot_info['Norm'] = []
        self.widget_ax = axes[0]
        self.plot_ax = axes[1]
        
        #Setting initial default values
        self.avg_reps_norm = True
        self.all_reps_norm = False
        
        #Plotting
        self.all_rep_lines_norm = []
        self._plot_all_rep_norm()
        self._plot_norm_graph()
        Plot_Norm.put_drift_annotation(self)
        
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
        """
        Will plot all the replicas as thin red lines.
        """
        self.all_rep_lines_norm = []
        for Dfilename in self.all_Dcoeff_data:
            coeffs, cols, timesteps, pops = self.all_Dcoeff_data[Dfilename]
            norms = np.sum(pops, axis=1)
            self.all_rep_lines_norm.append(self.plot_ax.plot(timesteps, norms, alpha=self.alpha, color='r', lw=1)[0])
        
        #Initialise the replica lines
        for line in self.all_rep_lines_norm:
            line.set_visible(self.all_reps_norm)
    
    @staticmethod
    def put_drift_annotation(self):
        """
        Will put an annotation of the drift value on the norm graph
        """
        coeffs, cols, timesteps, pops = self.all_Dcoeff_data_avg
        norms = np.sum(pops, axis=1)
        fit = np.polyfit(timesteps, norms, 1)
        errs = [1e-10]+[1e-6]*(len(timesteps)-1)
        fit2, pcov = curve_fit(linear_fit, timesteps, norms, p0=[fit[0], 1], sigma=errs)
        step = np.mean(np.diff(timesteps))
        text = r"Avg drift per rep = %.2g ps$^{-1}$"%(fit2[0]*(1000/step))
        
        y1 = np.polyval(fit2, timesteps)
        self.plot_ax.plot(timesteps, y1, 'k--', lw=0.5)
        
        all_norms = [np.sum(self.all_Dcoeff_data[i][3], axis=1) for i in self.all_Dcoeff_data]
        min_x, min_y = min(timesteps), np.min(all_norms)
        range_x, range_y = max(timesteps) - min(timesteps), np.max(all_norms) - np.min(all_norms)
        loc = (min_x+0.05*range_x, min_y + 0.85*range_y)
        
        self.norm_drift = fit[0]
        self.plot_ax.annotate(text, loc, fontsize=18)
        self.plot_info['Norm'].append(text[:text.find('ps')])

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