#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 12:17:40 2019

@author: oem
"""
from matplotlib.widgets import CheckButtons
import matplotlib.pyplot as plt


class K(object):
    """
    Will plot the K data.

    Inputs:
        axes  => a list of the axes to plot on. The first item should be the
                 widget axis the second will be the axis to plot the data.

    NOTE: This really is a bit of a hack, I should really go back and have a
          think of a more elegant to code this later.
    """
    def __init__(self, axes):
        if self.plot:
            K.widgAx = axes[0]
            K.plotAx = axes[1]

            K._plotAll(self)

            K.plotAx.set_ylabel("$X_{ctmqc, l}^{(I)}$")

    @staticmethod
    def _plotAll(self):
        """
        Will plot K showing all the replicas as thin lines (different colors
        are different states)
        """
        for filename in self.all_K:
            data, timesteps = self.all_K[filename]
            for i in range(self.num_states):
                K.plotAx.plot(timesteps, data[:, i],
                              alpha=0.7, lw=1,
                              color=self.colors[i])
