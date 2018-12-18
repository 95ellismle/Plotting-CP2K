#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 09:30:31 2018

@author: Sangeya
"""
import numpy as np
import matplotlib.cm as cm
from matplotlib.widgets import CheckButtons
import matplotlib.pyplot as plt




class Adiab_States(object):
    """
    Will plot the adiabatic states and handle the checkbuttons turning on and 
    off the 'fill between', 'all reps' and 'avg reps' options
    
    Inputs:
        * axes => The axes on which to plot (first is the checkbox axis, 
                                             second is the plotting axis)
    """
    
    def __init__(self, axes):
        if self.plot:
            self.AS_widg_ax, self.AS_plot_ax = axes
            self.state_cols_AS = [i for i in self.all_ad_ener_data_avg.columns if 'State' in i]
    
            
            # Set defaults
            self.plot_all_AS = False
            self.plot_avg_AS = True
            self.fill_between_AS = False
            
            # Plot everything
            self._plot_all_reps_AS()
            self._plot_avg_rep_AS()
            self._fill_between_AS()
            
            # Set the control
            if self._use_control: self._set_control_AS()
            
            # Finish up (make things look pretty)
            self.AS_plot_ax.set_ylabel(r"$E^{ad}_{l}$ [Ha]")
    
    
    def _plot_all_reps_AS(self):
        """
        Will plot all replicas as light thin lines on the graph
        """
        #Plot graphs
        self.AS_all_lines = []
        for Efilename in self.all_ad_ener_data:
            ad_ener_data = self.all_ad_ener_data[Efilename]
            for i, col in enumerate(self.state_cols_AS):
                ln, = self.AS_plot_ax.plot(ad_ener_data['Time'], 
                                           ad_ener_data[col], 
                                           alpha=self.alpha, 
                                           lw=0.7, 
                                           color=self.colors[i])
                self.AS_all_lines.append(ln)
        
        #Set initial visibility
        for line in self.AS_all_lines:
            line.set_visible(self.plot_all_AS)

    def _plot_avg_rep_AS(self):
        """
        Will plot the average replica as a thick line on the graph
        """
        #Plot the lines
        self.AS_avg_lines = []
        for i, col in enumerate(self.state_cols_AS):
            ln, = self.AS_plot_ax.plot(self.all_ad_ener_data_avg['Time'], 
                                       self.all_ad_ener_data_avg[col], 
                                       color=self.colors[i])
            self.AS_avg_lines.append(ln)
        
        #Set initial visibility
        for line in self.AS_avg_lines:
            line.set_visible(self.plot_avg_AS)
    
    def _fill_between_AS(self):
        """
        Will fill the adiabatic states with a color dependant on how close 
        they are to each other.
        """
        scaler = 2.5
        self.all_fill_bars = []
        #Fill colors originally
        for i, state in enumerate(self.state_cols_AS[:-1]):
            y1_y2 = self.all_ad_ener_data_avg[state] - self.all_ad_ener_data_avg[self.state_cols_AS[i+1]]
            max_diff, min_diff = np.max(y1_y2), np.min(y1_y2)
            normed_diffs = (y1_y2 - min_diff)/(max_diff - min_diff)
            for dt, timestep in enumerate(self.all_ad_ener_data_avg['Time'][:-1]):
                T = [timestep, self.all_ad_ener_data_avg['Time'][dt+1]]
                y1 = self.all_ad_ener_data_avg[state].iloc[[dt, dt+1]]
                y2 = self.all_ad_ener_data_avg[self.state_cols_AS[i+1]].iloc[[dt, dt+1]]
                diffy = scaler*(1-normed_diffs[dt])
                if diffy < 1:
                    color = cm.hot(diffy)
                    
                    ln = self.AS_plot_ax.fill_between(T, y1, y2, alpha=0.4, color=color, lw=0)
                    self.all_fill_bars.append(ln)
        
        # Set initial visibility
        for line in self.all_fill_bars:
            line.set_visible(self.fill_between_AS)
        
    def _set_control_AS(self):
        """
        Will connect the checkbutton control panel with the plotting axis.
        """
        self.check_AS = CheckButtons(self.AS_widg_ax, 
                                     ('all rep', 'avg', 'fill'), 
                                     (self.plot_all_AS, self.plot_avg_AS, self.fill_between_AS))
        self.check_AS.on_clicked(self._on_click_AS)
    
    def _on_click_AS(self, label):
        """
        Will handle button presses from the check buttons
        """
        if label == 'all rep':
            for line in self.AS_all_lines:
                line.set_visible(not line.get_visible())
        elif label == 'avg':
            for line in self.AS_avg_lines:
                line.set_visible(not line.get_visible())
        elif label == 'fill':
            for line in self.all_fill_bars:
                line.set_visible(not line.get_visible())
        plt.draw()
        
        


class Energy_Cons(object):
    """
    Will plot the kinetic, potential and total energy of the system.
    
    Inputs:
        * axes => The axes on which to plot (first is the checkbox axis, 
                                             second is the plotting axis)
    """
    
    def __init__(self, axes):
        if self.plot:
            Energy_Cons.widg_ax, Energy_Cons.plot_ax = axes
            
            Energy_Cons.plot_all_energy_drifts(self)
            Energy_Cons.plot_avg_energy_drift(self)
            
            Energy_Cons.plot_ax.set_ylabel("Energy Drift [Ha]")
    
    @staticmethod
    def plot_all_energy_drifts(self):
        """
        Will plot the all energy drifts on the Energy_Cons.plot_ax 
        """
        # Total energy
        for irep in self.all_tot_ener:
            data = self.all_tot_ener[irep]
            Energy_Cons.plot_ax.plot(data['Time'], 
                                     data['E_cons'], 
                                     'k-', 
                                     lw=0.5, 
                                     alpha=self.alpha)
    
    @staticmethod
    def plot_avg_energy_drift(self):
        """
        Will plot the average energy drift on the Energy_Cons.plot_ax 
        """
        Energy_Cons.plot_ax.plot(self.tot_ener_mean['Time'], 
                                 self.tot_ener_mean['E_cons'],
                                 'k-',
                                 lw=1.5)
    
