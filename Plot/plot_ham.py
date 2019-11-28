#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 12:23:18 2018

@author: mellis
"""

import matplotlib.pyplot as plt
from matplotlib.widgets import CheckButtons

#import numpy as np
import itertools as IT

class Site_Ener(object):
    """
    Will plot the Qlk against time graph.
    
    Inputs:
        axes  => a list of the axes to plot on. The first item_plot_site_ener should be the 
                 widget axis the second will be the axis to plot the data.
    """
    def __init__(self, axes):
        if self.plot:
            self.widget_ax = axes[0]
            self.plot_ax = axes[1]
            
            #Setting initial default values
            Site_Ener.avg_reps = True
            Site_Ener.all_reps = True
            
            #Plotting
            self._plot_all_rep_site_ener()
            self._plot_avg_site_ener()
            
            #Connect checkboxes to plot control
            if self._use_control:    self._set_site_ener_control()
            
            self.plot_ax.set_ylabel(r"$\Delta E_{i,j}$")
    
    @staticmethod
    def _check_settings_site_ener(label):
        if label == 'all replicas': #Pressed the all replicas button
            for line in Site_Ener.all_rep_lines:
                line.set_visible(not line.get_visible())
                
        elif label == 'average':
            for line in Site_Ener.avg_rep_lines:
                line.set_visible(not line.get_visible())
        plt.draw()
            
    #Will set the control panel for site_ener graph
    def _set_site_ener_control(self):
        self.check_site_ener = CheckButtons(self.widget_ax, ('all replicas', 'average'), (Site_Ener.all_reps, Site_Ener.avg_reps))
        self.check_site_ener.on_clicked(Site_Ener._check_settings_site_ener)
  
    def _plot_all_rep_site_ener(self):
       """
       Will plot all replicas site energies
       """
       Site_Ener.all_rep_lines = []
       for Hrep in self.all_ham_data:
           Hs, cols, timesteps = self.all_ham_data[Hrep]
           all_combs = IT.combinations(range(len(Hs[0])),2)
           for coli, (i,j) in enumerate(all_combs):
               site_ener = Hs[:,i,i] - Hs[:,j,j]
               ln, = self.plot_ax.plot(timesteps, site_ener, color=self.colors[coli], lw=0.7)
               Site_Ener.all_rep_lines.append(ln)
       
       for line in Site_Ener.all_rep_lines:
           line.set_visible(Site_Ener.all_reps)
           line.set_visible(Site_Ener.all_reps)
           
    def _plot_avg_site_ener(self):
        """
        Will plot the average site energy for all replicas.
        """
        Site_Ener.avg_rep_lines = []
        for Hrep in self.avg_ham_data:
           Hs, cols, timesteps = self.avg_ham_data[Hrep]
           all_combs = IT.combinations(range(len(Hs[0])),2)
           for coli, (i,j) in enumerate(all_combs):
               site_ener = Hs[:,i,i] - Hs[:,j,j]
               ln, = self.plot_ax.plot(timesteps, site_ener, color=self.colors[coli], label="%i-%i"%(i,j), lw=1.5)
               Site_Ener.avg_rep_lines.append(ln)
       
        for line in Site_Ener.all_rep_lines:
           line.set_visible(Site_Ener.avg_reps)
        self.plot_ax.legend()


class Coupling(object):
    """
    Will plot the couplings in the hamiltonian

    Inputs:
        * axes => axes to plot upon. First in list is the control panel axis.
                  Second is the plotting axis [list of plt.axes]
    """
    def __init__(self, axes):
        if self.plot:
            self.coup_widg_ax, self.coup_plot_ax = axes

            # Setting initial default values
            self.avg_reps_coup = True
            self.all_reps_coup = False

            # Plotting
            self.all_rep_lines_coup = []
            self._plot_all_rep_coup()
            self._plot_avg_coup()

            # Connect checkboxes to plot control
            if self._use_control:
                self._set_coup_control()

            self.coup_plot_ax.set_ylabel(r"$H_{ab}$ [meV]")

    def _check_settings_coup(self, label):
        if label == 'all replicas':  # Pressed the all replicas button
            for line in self.all_rep_lines_coup:
                line.set_visible(not line.get_visible())

        elif label == 'average':
            self.avg_line_coup.set_visible(
                                           not self.avg_line_coup.get_visible()
                                          )
        plt.draw()

    # Will set the control panel for coup graph
    def _set_coup_control(self):
        self.check_coup = CheckButtons(self.coup_widg_ax,
                                       ('all replicas', 'average'),
                                       (self.all_reps_coup,
                                        self.avg_reps_coup))

        self.check_coup.on_clicked(self._check_settings_coup)

    def _plot_all_rep_coup(self):
        """
        Will plot the coupling for each replica and save the lines in
        all_rep_lines_coup.
        """
        self.all_rep_lines_coup = []
        for Hkey in self.all_meta_ham:
            ln, = self.coup_plot_ax.plot(self.all_meta_ham[Hkey]['tsteps'],
                                         self.all_meta_ham[Hkey]['coup'],
                                         'k-', alpha=self.alpha)
            self.all_rep_lines_coup.append(ln)

        # Initialise the replica lines
        for line in self.all_rep_lines_coup:
            line.set_visible(self.all_reps_coup)

    def _plot_avg_coup(self):
        """
        Will plot the off-diagonal hamiltonian elements average over all reps.
        """
        avg_Hs, _, timesteps = self.avg_ham_data['avg_ham']
        avg_coups = avg_Hs[:, 0, 1]
        ln, = self.coup_plot_ax.plot(timesteps, avg_coups)
        self.avg_line_coup = ln

class H(object):
    """
    Will plot the couplings in the hamiltonian

    Inputs:
        * axes => axes to plot upon. First in list is the control panel axis.
                  Second is the plotting axis [list of plt.axes]
    """
    def __init__(self, axes):
        if self.plot:
            H.widg_ax, H.plot_ax = axes

            # Setting initial default values
            H.avg_reps = True
            H.all_reps = False

            # Plotting
            H.all_rep_lines = []
            H._plot_all_rep(self)
            #self._plot_avg_coup()

            ## Connect checkboxes to plot control
            #if self._use_control:
            #    self._set_coup_control()

            H.plot_ax.set_ylabel(r"$H$ [meV]")

    def _check_settings(self, label):
        if label == 'all replicas':  # Pressed the all replicas button
            for line in self.all_rep_lines_H:
                line.set_visible(not line.get_visible())

        elif label == 'average':
            self.avg_line_H.set_visible(
                                           not self.avg_line_H.get_visible()
                                          )
        plt.draw()

    # Will set the control panel for H graph
    def _set_control(self):
        self.check_H = CheckButtons(self.H_widg_ax,
                                       ('all replicas', 'average'),
                                       (self.all_reps_H,
                                        self.avg_reps_H))

        self.check_H.on_clicked(self._check_settings_H)

    @staticmethod
    def _plot_all_rep(self):
        """
        Will plot the Hling for each replica and save the lines in
        all_rep_lines_H.
        """
        for f in self.all_ham_data:
            data, _, time = self.all_ham_data[f]
            count = 0
            for i in range(data.shape[1]):
               for j in range(data.shape[1]):
                  H.plot_ax.plot(time, data[:, i, j], color=self.colors[count],
                                 alpha=self.alpha, lw=1)
                  count += 1


    def _plot_avg(self):
        """
        Will plot the off-diagonal hamiltonian elements average over all reps.
        """
        avg_Hs, _, timesteps = self.avg_ham_data['avg_ham']
        avg_Hs = avg_Hs[:, 0, 1]
        ln, = self.H_plot_ax.plot(timesteps, avg_Hs)
        self.avg_line_H = ln
