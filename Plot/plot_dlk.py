#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 16:11:09 2018

@author: mellis
"""

from matplotlib.widgets import Slider, CheckButtons
import matplotlib.pyplot as plt
import numpy as np
from collections import OrderedDict
import time

from load import load_dlk


class dlk(object):
    """
    Will plot the dlk against time graph.

    Inputs:
        axes  => a list of the axes to plot on. The first item should be the
                 widget axis the second will be the axis to plot the data.
    """
    def __init__(self, axes):
        if self.plot:
            dlk.widget_ax = axes[0]
            dlk.plot_ax = axes[1]

            # Setting initial default values
            dlk.X, dlk.Y, dlk.Z, dlk.Mag = [False, False, False, True]
            dlk._show_all, dlk._show_avg = True, False

            # Plotting
            dlk._plot_avg_dlk_t(self)
            dlk._plot_all(self)

#            # Connect checkboxes to plot control
#            if self._use_control:
#                dlk._set_dlk_t_control()

            dlk.plot_ax.set_ylabel(r"$\mathbf{d}_{lk, nu}^{(I)}$")

    @staticmethod
    def _check_settings_dlk_t(label):
        if label == 'X':
            dlk.X = not dlk.X
            for ln in dlk.Xlines:
                ln.set_visible(dlk.X)
        if label == 'Y':
            dlk.Y = not dlk.Y
            for ln in dlk.Ylines:
                ln.set_visible(dlk.Y)
        if label == 'Z':
            dlk.Z = not dlk.Z
            for ln in dlk.Zlines:
                ln.set_visible(dlk.Z)
        if label == 'Magnitude':
            dlk.Mag = not dlk.Mag
            for ln in dlk.Maglines:
                ln.set_visible(dlk.Mag)
        plt.draw()

    @staticmethod
    def _turn_on_atoms(Min, Max, type_switch='opposite'):
        """
        Will turn on atoms in a given range dependent on the current state of
        the X, Y, Z and Mag checkboxes

        Inputs:
            * Min => min atom to plot in range
            * Max => max atom to plot in range
            * type_switch => the visibility of the line (opposite means
                             opposite of the current state).
        """
        # First handle the X, Y, Z, Mag lines (averages)
        poss_lines = [dlk.Xlines, dlk.Ylines, dlk.Zlines, dlk.Maglines]
        all_switches = [dlk.X, dlk.Y, dlk.Z, dlk.Mag]
        for i, test in enumerate(all_switches):
            if test:
                for ln in poss_lines[i][Min: Max]:
                    if type_switch == 'opposite':
                        ln.set_visible(not ln.set_visible)
                    else:
                        ln.set_visible(type_switch)

        # Then handle the all replica lines
        for iat, linelist in enumerate(dlk.all_rep_lines[Min: Max]):
            for ln in linelist:
                if type_switch == 'opposite':
                    ln.set_visible(not ln.get_visible())
                else:
                    ln.set_visible(type_switch)
        plt.draw()

    @staticmethod
    def _submit_text(text):
        """
        Will decide which atoms to plot from the submitted text in the text box
        """
        if text == 'all':
            dlk._turn_on_atoms(0, len(dlk.Xlines), True)
        else:
            text = text.split(',')
            for item in text:
                minmax = item.split('-')
                if len(minmax) == 2:
                    Min, Max = minmax
                    if Max == ':':
                        Max = len(dlk.Xlines)
                    try:
                        Min = int(Min)
                        Max = int(Max)
                    except ValueError:
                        print("Tried to convert '%s' to ints" +
                              " but couldn't" % item)
                elif len(minmax) == 1:
                    try:
                        Min = int(minmax[0])
                        Max = int(minmax[0])+1
                    except ValueError:
                        print("Tried to convert '%s' to ints but" +
                              " couldn't" % item)
                if Max > len(dlk.Xlines):
                    Max = len(dlk.Xlines)
                dlk._turn_on_atoms(0, len(dlk.Xlines), False)  # Turn all off
                dlk._turn_on_atoms(Min, Max, True)  # Turn the correct ones on

    # Will set the control panel for dlk graph
    @staticmethod
    def _set_dlk_t_control():
        dlk.check_dlk_t = CheckButtons(dlk.widget_ax[0],
                                        ('X', 'Y', 'Z', 'Magnitude'),
                                        [dlk.X, dlk.Y, dlk.Z, dlk.Mag])
        dlk.check_dlk_t.on_clicked(dlk._check_settings_dlk_t)

#        dlk.textbox = TextBox(dlk.widget_ax[1], 'Atoms:', initial="all")
#        dlk.textbox.on_submit(dlk._submit_text)

    # Will plot the Quantum Momentum term
    @staticmethod
    def _plot_avg_dlk_t(self):
        """
        Will plot the dlk vs time for the average replica
        """
        ax = dlk.plot_ax

        dlk.Xlines = []
        dlk.Ylines = []
        dlk.Zlines = []
        dlk.Maglines = []

        dlk_timesteps = self.avg_dlk_data['time']  # grab timesteps from dlk_data list

        for iatom in self.atoms_to_plot:
            mask = self.avg_dlk_data['v'] == iatom
            mask = mask & (self.avg_dlk_data['l'] == 1)
            mask = mask & (self.avg_dlk_data['k'] == 2)
            
            dlkX = self.avg_dlk_data[mask]['dlk(x)']
            dlk_mag = dlkX**2

            if 'dlk(y)' in self.avg_dlk_data.columns:
               dlkY = self.avg_dlk_data[mask]['dlk(y)']
               dlk_mag += dlkY**2
            if 'dlk(z)' in self.avg_dlk_data.columns:
               dlkZ = self.avg_dlk_data[mask]['dlk(z)']
               dlk_mag += dlkZ**2

            dlk_mag = np.sqrt(dlk_mag)

            # Plot Mag
            ln, = ax.plot(dlk_timesteps,
                          dlk_mag,  # /np.max(dlk_mag[:, 0]),
                          '--',
                          color=self.colors[iatom-1])
            ln.set_visible(dlk.Mag and dlk._show_avg)
            dlk.Maglines.append(ln)

    # Will plot the all replica dlk_ts
    @staticmethod
    def _plot_all(self):
        """
        Will plot all the dlk vs time lines for each replica and save all the
        line in a 2D list called
        """
        dlk.all_rep_lines = [[] for i in range(self.num_active_atoms)]

        ax = dlk.plot_ax
        for dlk_filename in self.all_dlk_data:
            dlk_data = self.all_dlk_data[dlk_filename]  # list of all dlk_data

            for iatom in self.atoms_to_plot:
               mask = dlk_data['v'] == iatom
               mask = mask & (dlk_data['l'] == 1)
               mask = mask & (dlk_data['k'] == 2)
               
               dlkX = dlk_data[mask]['dlk(x)']
               dlk_timesteps = dlk_data[mask]['time']  # grab timesteps from dlk_data list

               sqrtOn = False
               cols = dlk_data.columns 
               if "dlk(y)" in cols or 'dlk(z)' in cols:
                  dlk_mag = dlkX**2
                  sqrtOn = True

               else: dlk_mag = dlkX

               if 'dlk(y)' in dlk_data.columns:
                  dlkY = dlk_data[mask]['dlk(y)']
                  dlk_mag += dlkY**2
               if 'dlk(z)' in dlk_data.columns:
                  dlkZ = dlk_data[mask]['dlk(z)']
                  dlk_mag += dlkZ**2

               if sqrtOn: dlk_mag = np.sqrt(dlk_mag)
               
               ln, = ax.plot(dlk_timesteps,
                             dlk_mag,  # /np.max(dlk_mag[:, 0]),
                             '-',
                             color=self.colors[iatom-1],
                             alpha=self.alpha,
                             lw=0.7)
               dlk.all_rep_lines[iatom-1].append(ln)
            ax.set_ylabel(r"|Q$_{lk}^{(I)}$|$^2$")

        # Initialise the replica lines
        for atlist in dlk.all_rep_lines:
            for line in atlist:
                line.set_visible(dlk._show_all)

