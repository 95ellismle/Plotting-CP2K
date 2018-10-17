#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 16 16:42:54 2018

@author: mellis
"""
import matplotlib.pyplot as plt
from matplotlib.widgets import CheckButtons, RadioButtons
import itertools as IT
import numpy as np

class fl_fk(object):
    """
    Will plot the values of f_l - f_k for all l,k combinations.
    
    Inputs:
        axes  => a list of the axes to plot on. The first item_plot_site_ener should be the 
                 widget axis the second will be the axis to plot the data.
    """
    def __init__(self, axes):
        fl_fk.widget_ax = axes[0]
        fl_fk.plot_ax = axes[1]
        
        #Setting initial default values
        fl_fk.xyz = [True, False, False]
        fl_fk.ats_to_plot = [True]
        fl_fk.all_reps = False
        
        #Plotting
#        fl_fk._plot_all_rep(self)
        fl_fk._plot_sum(self)
        
        #Connect checkboxes to plot control
        fl_fk.colors = self.colors
        if self._use_control:    fl_fk._set_control()
        
        fl_fk.plot_ax.set_ylabel(r"$|C_{k}|^2 |C_{l}|^2 ( f_{k} - f_{l} )$")
    
    @staticmethod
    def _plot_all_rep(self):
        """
        Will plot all replica f_l - f_k data.
        """
        fl_fk.all_rep_lines = []
        for Trep in self.all_tintf_data:
            (data, cols), timesteps = self.all_tintf_data[Trep]
            nstates = max(cols[0,:,1].astype(int))
            fl_fk.natom   = int(len(data[0])/nstates)
            all_combs = IT.combinations(range(nstates), 2)
            for l, k in all_combs:
                for iat in range(fl_fk.natom):
                    fk = data[cols[:,:,1] == str(k+1)]
                    fl = data[cols[:,:,1] == str(l+1)]
                    
                    fk = np.array([fk[12*dt:12*(dt+1)] for dt in range(len(data))])
                    fl = np.array([fl[12*dt:12*(dt+1)] for dt in range(len(data))])
                    
                    Fl_Fk = fl-fk
                    ln, = fl_fk.plot_ax.plot(timesteps, Fl_Fk[:,iat,0], color=self.colors[iat])
                    fl_fk.all_rep_lines.append(ln)
    
    @staticmethod
    def _plot_sum(self):
        """
        Will plot the average (over replicas) of the f_l - f_k data.
        """
        fl_fk.sum_rep_lines = [[],[],[]]
        for Trep in self.avg_tintf_data:
            (data, cols), timesteps = self.avg_tintf_data[Trep]
            nstates = max(cols[0,:,1].astype(int))
            fl_fk.natom   = int(len(data[0])/nstates)
            all_combs = IT.combinations(range(nstates), 2)
            for l, k in all_combs:
                for iat in range(fl_fk.natom):
                    fk = data[cols[:,:,1] == str(k+1)]
                    fl = data[cols[:,:,1] == str(l+1)]
                    
                    fk = np.array([fk[12*dt:12*(dt+1)] for dt in range(len(data))])
                    fl = np.array([fl[12*dt:12*(dt+1)] for dt in range(len(data))])
                    
                    Fl_Fk = fl-fk
                    
                    cl_ck = self.all_Acoeff_data_avg[-1][:,l] \
                          * self.all_Acoeff_data_avg[-1][:,k] \
                          * len(self.all_tintf_data)
                    
                    for idim in range(3):
                        ln, = fl_fk.plot_ax.plot(timesteps, 
                                             Fl_Fk[:,iat,idim]*cl_ck, 
                                             color=self.colors[idim])
                        fl_fk.sum_rep_lines[idim].append(ln)
        
        for idim, lines in enumerate(fl_fk.sum_rep_lines):
            for line in lines:
                line.set_visible(fl_fk.xyz[idim])
        
    @staticmethod
    def _set_control():
        """
        Will connect the controls to the plotting
        
        """
        fl_fk.ats_to_plot = fl_fk.ats_to_plot*fl_fk.natom
        
        fl_fk.xyz_control = CheckButtons(fl_fk.widget_ax[0], ['X','Y','Z'], fl_fk.xyz)
        fl_fk.xyz_control.on_clicked(fl_fk._cart_control)
        
#        fl_fk.widget_ax[1].set_title("Color by:", fontsize=15)
        fl_fk.color_control = RadioButtons(fl_fk.widget_ax[1], ['XYZ', 'atom'], [True, False])
        fl_fk.color_control.on_clicked(fl_fk._color_control)
        
        fl_fk.atom_control = CheckButtons(fl_fk.widget_ax[2], 
                                          ['atom %i'%(i+1) for i in range(fl_fk.natom)], 
                                          fl_fk.ats_to_plot)
        fl_fk.atom_control.on_clicked(fl_fk._atom_control)
        
    
    @staticmethod
    def _atom_control(label):
        atom_num = int(label.replace("atom", ""))-1
        fl_fk.ats_to_plot[atom_num] = not fl_fk.ats_to_plot[atom_num]
        for idim, lines in enumerate(fl_fk.sum_rep_lines):
            lines[atom_num].set_visible(fl_fk.ats_to_plot[atom_num] and fl_fk.xyz[idim])            
        plt.draw()
        
    @staticmethod
    def _color_control(label):
        """
        Will control what coloring to use in the plots
        """
        if label == "XYZ":
            for idim, lines in enumerate(fl_fk.sum_rep_lines):
                for line in lines:
                    line.set_color(fl_fk.colors[idim])
        if label == "atom":
            for lines in fl_fk.sum_rep_lines:
                for iat, line in enumerate(lines):
                    line.set_color(fl_fk.colors[iat])
        
        plt.draw()
    
    @staticmethod
    def _cart_control(label):
        """
        Will determine what coloring scheme to use in the plots
        """
        if label == 'X':
            fl_fk.xyz[0] = not fl_fk.xyz[0]
            fl_fk.xyz[0] = not fl_fk.xyz[0]
            for iat, line in enumerate(fl_fk.sum_rep_lines[0]):
                line.set_visible(fl_fk.ats_to_plot[iat] and fl_fk.xyz[0])
        elif label == 'Y':
            fl_fk.xyz[1] = not fl_fk.xyz[1]
            for iat, line in enumerate(fl_fk.sum_rep_lines[1]):
                line.set_visible(fl_fk.ats_to_plot[iat] and fl_fk.xyz[1])
        elif label == 'Z':
            fl_fk.xyz[2] = not fl_fk.xyz[2]
            for iat, line in enumerate(fl_fk.sum_rep_lines[2]):
                line.set_visible(fl_fk.ats_to_plot[iat] and fl_fk.xyz[2])
        plt.draw()
        