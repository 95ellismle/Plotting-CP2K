#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 16 16:42:54 2018

@author: mellis
"""
import matplotlib.pyplot as plt
from matplotlib.widgets import CheckButtons, RadioButtons, TextBox
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
        fl_fk.get_metadata(self)
        
        #Plotting
        fl_fk._plot_sum_CC(self)
        
        #Connect checkboxes to plot control
        fl_fk.colors = self.colors
        if self._use_control:    fl_fk._set_control()
        
        fl_fk.plot_ax.set_ylabel(r"$|C_{l}|^2 |C_{k}|^2 ( f_{k} - f_{l} )$", fontsize=28)
    
    @staticmethod
    def _plot_sum(self):
        """
        Will plot the average (over replicas) of the f_l - f_k data.
        """
        fl_fk.sum_rep_lines = [[],[],[]]
        for Trep in self.sum_tintf_data:
            (data, cols), timesteps = self.sum_tintf_data[Trep]
            all_combs = IT.combinations(range(fl_fk.nstates), 2)
            for l, k in all_combs:
                fk = np.array([i[cols[0,:,1] == str(k+1)] for i in data])
                fl = np.array([i[cols[0,:,1] == str(l+1)] for i in data])
                
                self.Fl_Fk = (fl-fk)
                for iat in range(fl_fk.natom):
                    for idim in range(3):
                        ln, = fl_fk.plot_ax.plot(timesteps, 
                                             self.Fl_Fk[:,iat,idim], 
                                             color=self.colors[idim], 
                                             lw=1.2)
                        fl_fk.sum_rep_lines[idim].append(ln)
        
        for idim, lines in enumerate(fl_fk.sum_rep_lines):
            for line in lines:
                line.set_visible(fl_fk.xyz[idim])
                
    @staticmethod
    def get_metadata(self):
        """
        Will get num atoms, num states etc...
        """
        for Trep in self.sum_tintf_data:
            (data, cols), timesteps = self.sum_tintf_data[Trep]
            fl_fk.nstates = max(cols[0,:,1].astype(int))
            fl_fk.natom   = int(len(data[0])/fl_fk.nstates)
            return 

    @staticmethod
    def _plot_sum_CC(self):
        """
        Will plot the \sum_{J}^{Nrep} C_k C_l *(f_l - f_k) term
        """
        data_dict, timesteps = self.sum_tintf_CC_data
        fl_fk.sum_rep_lines = [[],[],[]]
        for Tkey in data_dict:
            data = data_dict[Tkey]
            for iat in range(fl_fk.natom):
                for idim in range(3):
                    ln, = fl_fk.plot_ax.plot(timesteps, 
                                             data[:,iat,idim], 
                                             color=self.colors[iat], 
                                             lw=1.2)
                    fl_fk.sum_rep_lines[idim].append(ln)
        fl_fk.Xlines, fl_fk.Ylines, fl_fk.Zlines = fl_fk.sum_rep_lines
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
        fl_fk.color_control = RadioButtons(fl_fk.widget_ax[1], ['atom', 'xyz'])
        fl_fk.color_control.on_clicked(fl_fk._color_control)
        
        fl_fk.textbox = TextBox(fl_fk.widget_ax[2], 'Atoms:', initial="all")
        fl_fk.textbox.on_submit(fl_fk._submit_text)
        
    
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
        if label == "xyz":
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
    
    @staticmethod
    def _turn_on_atoms(Min, Max, type_switch='opposite'):
        """
        Will turn on atoms in a given range dependent on the current state of 
        the X, Y and Z checkboxes
        """
        poss_lines   = [fl_fk.Xlines, fl_fk.Ylines, fl_fk.Zlines]
        all_switches = fl_fk.xyz
        for i, test in enumerate(all_switches):
            if test:
                for ln in poss_lines[i][ Min:Max ]:
                    if type_switch == 'opposite':
                        ln.set_visible(not ln.set_visible)
                    else:
                        ln.set_visible(type_switch)
        plt.draw()
        
    @staticmethod
    def _submit_text(text):
        """
        Will decide which atoms to plot from the submitted text in the text box
        """
        if text == 'all':
            fl_fk._turn_on_atoms(0,len(fl_fk.Xlines), True)
        else:
            text = text.split(',')
            for item in text:
                minmax = item.split('-')
                if len(minmax) == 2:
                    Min, Max = minmax
                    if Max == ':':
                        Max = len(fl_fk.Xlines)
                    try:
                        Min = int(Min)
                        Max = int(Max)
                    except:
                        print("Tried to convert '%s' to ints but couldn't"%item)
                elif len(minmax) == 1:
                    try:
                        Min = int(minmax[0])
                        Max = int(minmax[0])+1
                    except:
                        print("Tried to convert '%s' to ints but couldn't"%item)
                if Max > len(fl_fk.Xlines):
                    Max = len(fl_fk.Xlines)
                fl_fk._turn_on_atoms(0,len(fl_fk.Xlines), False)
                fl_fk._turn_on_atoms(Min, Max, True)
        
