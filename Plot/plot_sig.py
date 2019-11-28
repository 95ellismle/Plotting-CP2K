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


class Sig(object):
    """
    Will plot the Qlk against time graph.

    Inputs:
        axes  => a list of the axes to plot on. The first item should be the
                 widget axis the second will be the axis to plot the data.
    """
    def __init__(self, axes):
        if self.plot:
            Sig.widget_ax = axes[0]
            Sig.plot_ax = axes[1]

            # Setting initial default values
            Sig.X, Sig.Y, Sig.Z, Sig.Mag = [False, False, False, True]
            Sig._show_all, Sig._show_avg = True, False

            # Plotting
            #Sig._plot_avg_Qlk_t(self)
            Sig._plot_all(self)

#            # Connect checkboxes to plot control
#            if self._use_control:
#                Sig._set_Qlk_t_control()

            Sig.plot_ax.set_ylabel(r"$\sigma_{\nu}^{(I)}$ [bohr]")


    # Will plot the all replica sigmas
    @staticmethod
    def _plot_all(self):
        """
        Will plot all the Qlk vs time lines for each replica and save all the
        line in a 2D list called
        """
        Sig.all_rep_lines = [[] for i in range(self.num_active_atoms)]

        ax = Sig.plot_ax
        for f in self.all_sigma:
            data, time = self.all_sigma[f]
            ax.plot(time, data, 'r.')

