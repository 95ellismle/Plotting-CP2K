#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 16 16:42:54 2018

@author: mellis
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.widgets import CheckButtons, RadioButtons  # ,TextBox
import itertools as IT
import numpy as np

from load import load_tintf


class AdMom(object):
    """
    Will plot history dependent force term for each state
    """
    def __init__(self, axes):
        if self.plot:
            AdMom.widget_ax = axes[0]
            AdMom.plot_ax = axes[1]
            AdMom.plot_ax.set_autoscale_on(True)
            
            AdMom._show_all = True
    
            AdMom._plot_all(self)
            
            AdMom.plot_ax.set_ylabel(r"$\mathbf{f}_{l, \nu}^{(I)}$ [bohr$^{-1}$]")

    # Will plot the all replica adiab momenta
    @staticmethod
    def _plot_all(self):
        """
        Will plot all the f vs time lines for each replica and save all the
        line in a 2D list called
        """
        AdMom.all_rep_lines = [[] for i in range(self.num_active_atoms)]
 
        ax = AdMom.plot_ax
        for Fname in self.all_adMom_data:
            data = self.all_adMom_data[Fname]  # list of all data
 
            for iatom in self.atoms_to_plot:
               mask = data['v'] == iatom
               mask = mask & (data['l'] == 1)
               masked_data = data[mask]
               time = masked_data["time"]
 
               adMomX = masked_data['f(x)']
               if all(j not in data.columns for j in ('f(y)','f(z)')):
                  moreThan1Dim = False
                  adMom_mag = adMomX
               else:
                  moreThan1Dim = True
                  adMom_mag = adMomX**2
 
               if 'f(y)' in data.columns:
                  adMomY = masked_data['f(y)']
                  adMom_mag += adMomY**2
               if 'f(z)' in data.columns:
                  adMomZ = masked_data['f(z)']
                  adMom_mag += adMomZ**2
 
               if moreThan1Dim:
                  adMom_mag = np.sqrt(adMom_mag)
 
               ln, = ax.plot(time,
                             adMom_mag,
                             '-',
                             color=self.colors[iatom-1],
                             alpha=self.alpha,
                             lw=0.7)
               AdMom.all_rep_lines[iatom-1].append(ln)
            ax.set_ylabel(r"|$\mathbf{f}_{l, \nu}^{(I)}$|$^2$")
 
        # Initialise the replica lines
        for atlist in AdMom.all_rep_lines:
            for line in atlist:
                line.set_visible(AdMom._show_all)
