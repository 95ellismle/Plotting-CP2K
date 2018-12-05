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
        if self.plot:
            fl_fk.widget_ax = axes[0]
            fl_fk.plot_ax = axes[1]
            fl_fk.plot_ax.set_autoscale_on(True)
            
            #Setting initial default values
            fl_fk.xyz = [True, False, False]
            fl_fk.ats_to_plot = [True]
            fl_fk.all_reps = False
            get_metadata(self, fl_fk)
            
            
            #Plotting
            fl_fk._set_DATA(self)
            plot_sum(self, fl_fk)
            
            #Connect checkboxes to plot control
            fl_fk.colors = self.colors
            if self._use_control:    set_control(fl_fk)
            
            fl_fk.plot_ax.axhline(0, color='k', lw=0.5)
            fl_fk.plot_ax.set_ylabel(r"$\sum_{J} ( f_{k}^{J} - f_{l}^{J} )$ [$\frac{Ha}{bohr}$]", fontsize=28)

    @staticmethod
    def _set_DATA(self):
        """
        Will set the fl_fk.DATA variable before plotting. This should be changed
        in child classes to the relevant variable.
        """
        fl_fk.DATA = self.sum_tintf_data  
            
    
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
    def scaley_ax(ax):
        """
        Will get the ydata from an axis and scale the lims to fit the data
        """
        print(ax.get_ydata())

        
    @staticmethod
    def _submit_text(text):
        """
        Will decide which atoms to plot from the submitted text in the text box
        """
        if text == 'all':
            turn_on_atoms(0,len(fl_fk.Xlines), fl_fk, True)
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
                        Max = 1e10
                        Min = 0
                if Max > len(fl_fk.Xlines):
                    Max = len(fl_fk.Xlines)
                turn_on_atoms(0,len(fl_fk_CC.Xlines), fl_fk, False)
                turn_on_atoms(Min, Max, fl_fk, True)
        

class fl_fk_CC(object):
    """
    Will plot the values of |C_l|^2 |C_k|^2 (f_l - f_k) for all l,k combinations.
    
    Inputs:
        axes  => a list of the axes to plot on. The first item_plot_site_ener should be the 
                 widget axis the second will be the axis to plot the data.
    """
    
    def __init__(self, axes):
        if self.plot:
            fl_fk_CC.widget_ax = axes[0]
            fl_fk_CC.plot_ax = axes[1]
            fl_fk_CC.plot_ax.set_autoscale_on(True)
            
            #Setting initial default values
            fl_fk_CC.xyz = [True, False, False]
            fl_fk_CC.ats_to_plot = [True]
            fl_fk_CC.all_reps = False
            get_metadata(self, fl_fk_CC)
            
            
            #Plotting
            fl_fk_CC._set_DATA(self)
            plot_sum(self, fl_fk_CC)
            
            #Connect checkboxes to plot control
            fl_fk_CC.colors = self.colors
            if self._use_control:    set_control(fl_fk_CC)
            fl_fk_CC.plot_ax.axhline(0, color='k', lw=0.5)
            fl_fk_CC.plot_ax.set_ylabel(r"$\sum_{J} |C_{l}^{J}|^2 |C_{k}^{J}|^2 ( f_{k}^{J} - f_{l}^{J} )$ [$\frac{Ha}{bohr}$]", fontsize=28)

    @staticmethod
    def _set_DATA(self):
        """
        Will set the fl_fk_CC.DATA variable before plotting. This should be changed
        in child classes to the relevant variable.
        """
        fl_fk_CC.DATA = self.sum_tintf_CC_data  

    @staticmethod
    def _plot_sum(self, DATA):
        """
        Will plot the \sum_{J}^{Nrep} C_k C_l *(f_l - f_k) term
        """
        data_dict, timesteps = DATA
        fl_fk_CC.sum_rep_lines = [[],[],[]]
        for Tkey in data_dict:
            data = data_dict[Tkey]
            for iat in range(fl_fk_CC.natom):
                for idim in range(3):
                    ln, = fl_fk_CC.plot_ax.plot(timesteps, 
                                             data[:,iat,idim], 
                                             color=self.colors[iat], 
                                             lw=1.2)
                    fl_fk_CC.sum_rep_lines[idim].append(ln)
        fl_fk_CC.Xlines, fl_fk_CC.Ylines, fl_fk_CC.Zlines = fl_fk_CC.sum_rep_lines
        for idim, lines in enumerate(fl_fk_CC.sum_rep_lines):
            for line in lines:
                line.set_visible(fl_fk_CC.xyz[idim])

    @staticmethod
    def _set_control():
        """
        Will connect the controls to the plotting
        
        """
        fl_fk_CC.ats_to_plot = fl_fk_CC.ats_to_plot*fl_fk_CC.natom
        
        fl_fk_CC.xyz_control = CheckButtons(fl_fk_CC.widget_ax[0], ['X','Y','Z'], fl_fk_CC.xyz)
        fl_fk_CC.xyz_control.on_clicked(fl_fk_CC._cart_control)
        
#        fl_fk_CC.widget_ax[1].set_title("Color by:", fontsize=15)
        fl_fk_CC.color_control = RadioButtons(fl_fk_CC.widget_ax[1], ['atom', 'xyz'])
        fl_fk_CC.color_control.on_clicked(fl_fk_CC._color_control)
        
        fl_fk_CC.textbox = TextBox(fl_fk_CC.widget_ax[2], 'Atoms:', initial="all")
        fl_fk_CC.textbox.on_submit(fl_fk_CC._submit_text)
        
    
    @staticmethod
    def _atom_control(label):
        atom_num = int(label.replace("atom", ""))-1
        fl_fk_CC.ats_to_plot[atom_num] = not fl_fk_CC.ats_to_plot[atom_num]
        for idim, lines in enumerate(fl_fk_CC.sum_rep_lines):
            lines[atom_num].set_visible(fl_fk_CC.ats_to_plot[atom_num] and fl_fk_CC.xyz[idim])            
        plt.draw()
        
    @staticmethod
    def _color_control(label):
        """
        Will control what coloring to use in the plots
        """
        if label == "xyz":
            for idim, lines in enumerate(fl_fk_CC.sum_rep_lines):
                for line in lines:
                    line.set_color(fl_fk_CC.colors[idim])
        if label == "atom":
            for lines in fl_fk_CC.sum_rep_lines:
                for iat, line in enumerate(lines):
                    line.set_color(fl_fk_CC.colors[iat])
        plt.draw()
    
    @staticmethod
    def _cart_control(label):
        """
        Will determine what coloring scheme to use in the plots
        """
        if label == 'X':
            fl_fk_CC.xyz[0] = not fl_fk_CC.xyz[0]
            for iat, line in enumerate(fl_fk_CC.sum_rep_lines[0]):
                line.set_visible(fl_fk_CC.ats_to_plot[iat] and fl_fk_CC.xyz[0])
        elif label == 'Y':
            fl_fk_CC.xyz[1] = not fl_fk_CC.xyz[1]
            for iat, line in enumerate(fl_fk_CC.sum_rep_lines[1]):
                line.set_visible(fl_fk_CC.ats_to_plot[iat] and fl_fk_CC.xyz[1])
        elif label == 'Z':
            fl_fk_CC.xyz[2] = not fl_fk_CC.xyz[2]
            for iat, line in enumerate(fl_fk_CC.sum_rep_lines[2]):
                line.set_visible(fl_fk_CC.ats_to_plot[iat] and fl_fk_CC.xyz[2])
        plt.draw()
    
    @staticmethod
    def scaley_ax(ax):
        """
        Will get the ydata from an axis and scale the lims to fit the data
        """
        print(ax.get_ydata())
            
    @staticmethod
    def _submit_text(text):
        """
        Will decide which atoms to plot from the submitted text in the text box
        """
        if text == 'all':
            turn_on_atoms(0,len(fl_fk_CC.Xlines), fl_fk_CC, True)
        else:
            text = text.split(',')
            for item in text:
                minmax = item.split('-')
                if len(minmax) == 2:
                    Min, Max = minmax
                    if Max == ':':
                        Max = len(fl_fk_CC.Xlines)
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
                        Max = 1e10
                        
                if Max > len(fl_fk_CC.Xlines):
                    Max = len(fl_fk_CC.Xlines)
                turn_on_atoms(0,len(fl_fk_CC.Xlines), fl_fk_CC, False)
                turn_on_atoms(Min, Max, fl_fk_CC, True)


"""
Defining some shared methods between the 2 fl_fk classes. This is a quick hack,
they should really have a shared parent... I might get round to making one soon
"""

def plot_sum(self, class_obj):
    """
    Will plot the \sum_{J}^{Nrep} C_k C_l *(f_l - f_k) term
    """
    data_dict, timesteps = class_obj.DATA
    class_obj.sum_rep_lines = [[],[],[]]
    for Tkey in data_dict:
        data = data_dict[Tkey]
        for iat in range(class_obj.natom):
            for idim in range(3):
                ln, = class_obj.plot_ax.plot(timesteps, 
                                         data[:,iat,idim], 
                                         color=self.colors[iat], 
                                         lw=1.2)
                class_obj.sum_rep_lines[idim].append(ln)
    class_obj.Xlines, class_obj.Ylines, class_obj.Zlines = class_obj.sum_rep_lines
    for idim, lines in enumerate(class_obj.sum_rep_lines):
        for line in lines:
            line.set_visible(class_obj.xyz[idim])

def turn_on_atoms(Min, Max, class_obj, type_switch='opposite'):
    """
    Will turn on atoms in a given range dependent on the current state of 
    the X, Y and Z checkboxes
    """
    poss_lines   = [class_obj.Xlines, class_obj.Ylines, class_obj.Zlines]
    all_switches = class_obj.xyz
    for i, test in enumerate(all_switches):
        if test:
            for ln in poss_lines[i][ Min:Max ]:
                if type_switch == 'opposite':
                    ln.set_visible(not ln.set_visible)
                else:
                    ln.set_visible(type_switch)
#        class_obj.plot_ax.relim()
#        class_obj.plot_ax.autoscale_view(True,True,True)
    plt.draw()

def get_metadata(self, class_obj):
    """
    Will get num atoms, num states etc...
    """
    class_obj.nsteps = len(self.sum_tintf_data[1])
    class_obj.nstates = len(list(set([int(j) for i in self.sum_tintf_data[0].keys() for j in i])))
    class_obj.natom   = np.shape(self.sum_tintf_data[0]['01'])[1]
    return 


def set_control(class_obj):
    """
    Will connect the controls to the plotting
    
    """
    class_obj.ats_to_plot = class_obj.ats_to_plot*class_obj.natom
    
    class_obj.xyz_control = CheckButtons(class_obj.widget_ax[0], ['X','Y','Z'], class_obj.xyz)
    class_obj.xyz_control.on_clicked(class_obj._cart_control)
    
#        class_obj.widget_ax[1].set_title("Color by:", fontsize=15)
    class_obj.color_control = RadioButtons(class_obj.widget_ax[1], ['atom', 'xyz'])
    class_obj.color_control.on_clicked(class_obj._color_control)
    
    class_obj.textbox = TextBox(class_obj.widget_ax[2], 'Atoms:', initial="all")
    class_obj.textbox.on_submit(class_obj._submit_text)