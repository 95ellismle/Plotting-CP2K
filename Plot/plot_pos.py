#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 30 13:15:54 2019

@author: mellis
"""

from matplotlib.widgets import CheckButtons
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


def linear_fit(x, m, c):
    return m * x + c


class PlotPos(object):
    """
    Will plot the normalisation graph.

    Inputs:
        axes  => a list of the axes to plot on. The first item should be the
                 widget axis the second will be the axis to plot the data.
    """
    def __init__(self, axes):
        if self.plot:
            PlotPos.widget_ax = axes[0]
            PlotPos.plot_ax = axes[1]

            # Setting initial default values
            PlotPos.__show_avg_reps = True
            PlotPos.__show_all_reps = True

            # Will do the plotting
            PlotPos.__plot_all(self)
            # Connect checkboxes to plot control
            if self._use_control:
                PlotPos.__set_control(self)

            PlotPos.plot_ax.set_ylabel(r"Pos")

    @staticmethod
    def __button_control(self, label):
        """
        Will handle the button events
        """

    # Will set the control panel for norm graph
    @staticmethod
    def __set_control(self):
        """
        Will put the control panel on the axis and connect the off and on
        functions to it.
        """

    # Will plot normal of diabatic coeffs
    @staticmethod
    def __plot_all(self):
        """
        Will plot the positions of atoms for all trajectories. Will sum the
        populations and plot on the norm axis.
        """
        for fileName in self.all_pos_data:
            (data, cols), timesteps = self.all_pos_data[fileName]

            # This is a horrible hack should be improved!
            #  * Make more robust way to get active atoms
            #  * Use numpy fancy indexing instead of list comprehension
            mask = cols != 'Ne'
            activeAtoms = [[pos for cond, pos in zip(mask[step], data[step])
                            if cond]
                           for step in range(len(data))]
            activeAtoms = np.array(activeAtoms)

            # Plot atoms
            for iat in self.atoms_to_plot:
                x = activeAtoms[:, iat-1, 0]
                y = activeAtoms[:, iat-1, 1]
                z = activeAtoms[:, iat-1, 2]
                mag = np.sqrt(x**2 + y**2 + z**2)

                PlotPos.plot_ax.plot(timesteps,
                                     mag,
                                     color=self.colors[iat-1],
                                     alpha=self.alpha,
                                     lw=0.5)

    # Will plot normal of diabatic coeffs
    @staticmethod
    def __plot_avg(self):
        """
        Will plot the positions of atoms for all trajectories. Will sum the
        populations and plot on the norm axis.
        """


class PlotPosSig(object):
    """
    Will plot the normalisation graph.

    Inputs:
        axes  => a list of the axes to plot on. The first item should be the
                 widget axis the second will be the axis to plot the data.
    """
    def __init__(self, axes):
        if self.plot:
            PlotPosSig.widget_ax = axes[0]
            PlotPosSig.plot_ax = axes[1]

            # Setting initial default values
            PlotPosSig.__show_avg_reps = True
            PlotPosSig.__show_all_reps = True

            # Will do the plotting
            PlotPosSig.__plot_all(self)
            # Connect checkboxes to plot control
            if self._use_control:
                PlotPosSig.__set_control(self)

            PlotPosSig.plot_ax.set_ylabel(r"$\alpha_{\nu}^{(I)} R_{\nu}^{(I)}$")
#            PlotPosSig.plot_ax.set_ylim([5, 1050])

    @staticmethod
    def __button_control(self, label):
        """
        Will handle the button events
        """

    # Will set the control panel for norm graph
    @staticmethod
    def __set_control(self):
        """
        Will put the control panel on the axis and connect the off and on
        functions to it.
        """

    # Will pos/2(sigma**2)
    @staticmethod
    def __plot_all(self):
        """
        Will plot the positions of atoms for all trajectories. Will sum the
        populations and plot on the norm axis.
        """
        for pFName, sFName in zip(self.all_pos_data, self.all_sigma):
            (posData, cols), timesteps = self.all_pos_data[pFName]
            sigma, timesteps = self.all_sigma[sFName]

            # This is a horrible hack should be improved!
            #  * Make more robust way to get active atoms
            #  * Use numpy fancy indexing instead of list comprehension
            mask = cols != 'Ne'
            activeAtoms = [[pos for cond, pos in zip(mask[step], posData[step])
                            if cond]
                           for step in range(len(posData))]
            activeAtoms = np.array(activeAtoms)

            # Plot atoms
            for iat in self.atoms_to_plot:
                x = activeAtoms[:, iat-1, 0]
                y = activeAtoms[:, iat-1, 1]
                z = activeAtoms[:, iat-1, 2]
                mag = np.sqrt(x**2 + y**2 + z**2)
                sigma_prefactor = 1 / (2*(sigma[:, iat]**2))

                PlotPosSig.plot_ax.plot(timesteps,
                                        mag * sigma_prefactor,
                                        color=self.colors[iat-1],
                                        alpha=self.alpha,
                                        lw=0.5)

    # Will plot normal of diabatic coeffs
    @staticmethod
    def __plot_avg(self):
        """
        Will plot the positions of atoms for all trajectories. Will sum the
        populations and plot on the norm axis.
        """
