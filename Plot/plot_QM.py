#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 16:11:09 2018

@author: mellis
"""

from matplotlib.widgets import CheckButtons
import matplotlib.pyplot as plt
import numpy as np


class QM(object):
    """
    Will plot the Quantum Momentum.
            
    Inputs:
        axes  => a list of the axes to plot on. The first item should be the 
                 widget axis the second will be the axis to plot the data.  
    """
    def __init__(self, axes):
        self.qm_widg_ax, self.qm_plot_ax = axes
        
        self._plot_avg_qm_vs_pos()
    
    def _plot_avg_qm_vs_pos(self):
        """
        Will plot the quantum momentum vs position.
        """
        avg_pos_mag = np.zeros([len(self.all_pos_data[self.all_pos_data.keys()[0]]   [0][0])])
        for pos_key in self.all_pos_data:
            self.all_pos_data[pos_key]
#            print(np.linalg.norm(self.all_pos_data[pos_key][0][0], axis=2))
        avg_pos_mag /= len(self.all_pos_data)
#        np.linalg.norm(p.all_pos_data['run-pos-1-1.xyz'][0][0], axis=2)
#        self.qm_plot_ax.plot([1,2,3],[2,3,4])
        plt.close()