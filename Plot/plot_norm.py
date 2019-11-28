#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  3 15:53:35 2018

@author: mellis
"""
from matplotlib.widgets import CheckButtons
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


def linear_fit(x, m, c):
    return m * x + c


class Plot_Norm(object):
    """
    Will plot the normalisation graph.

    Inputs:
        axes  => a list of the axes to plot on. The first item should be the
                 widget axis the second will be the axis to plot the data.
    """
    def __init__(self, axes):
        self.plot_info['Norm'] = []
        if self.plot:
            if '|c|^2' in self.load_params:
               Plot_Norm.coeff_data = self.all_Acoeff_data
               Plot_Norm.avg_coeff_data = self.all_Acoeff_data_avg
               label = "C"
            else:
               Plot_Norm.coeff_data = self.all_Dcoeff_data
               Plot_Norm.avg_coeff_data = self.all_Dcoeff_data_avg
               label = "u"
            self.widget_ax = axes[0]
            self.plot_ax = axes[1]

            # Setting initial default values
            self.avg_reps_norm = True
            self.all_reps_norm = True

            # Plotting
            self.all_rep_lines_norm = []

            self._plot_all_rep_norm()
            self._plot_norm_graph()
            # Connect checkboxes to plot control
            if self._use_control:
                self._set_norm_control()

            self.plot_ax.set_ylabel(r"$\sum_k |%s_k^{I}|^2$" % label)

        Plot_Norm._get_all_rep_norms(self)
        Plot_Norm.put_drift_annotation(self)

    def _check_settings_norm(self, label):
        if label == 'all replicas':  # Pressed the all replicas button
            for line in self.all_rep_lines_norm:
                line.set_visible(not line.get_visible())

        elif label == 'average':
            self.avg_line_norm.set_visible(
                                           not self.avg_line_norm.get_visible()
                                          )
        plt.draw()

    # Will set the control panel for norm graph
    def _set_norm_control(self):
        self.check_norm = CheckButtons(self.widget_ax,
                                       ('all replicas', 'average'),
                                       (self.all_reps_norm, self.avg_reps_norm)
                                       )
        self.check_norm.on_clicked(self._check_settings_norm)

    # Will plot the all replica norms
    def _plot_all_rep_norm(self):
        """
        Will plot all the replicas as thin red lines.
        """

        self.all_rep_lines_norm = []
        for Dfilename in Plot_Norm.coeff_data:
            coeffs, cols, timesteps, pops = Plot_Norm.coeff_data[Dfilename]
            norms = np.sum(pops, axis=1)
            self.all_rep_lines_norm.append(self.plot_ax.plot(timesteps,
                                                             norms,
                                                             alpha=self.alpha,
                                                             color='r',
                                                             lw=1
                                                             )[0])

        # Initialise the replica lines
        for line in self.all_rep_lines_norm:
            line.set_visible(self.all_reps_norm)

    @staticmethod
    def _get_all_rep_norms(self):
        """
        Will get all the norm drifts from each replica.
        """
        unit_conv = 1 
        if self.units == "au":
             unit_conv = 0.02418884

        self.norm_drift_per_rep = []
        for Dfilename in Plot_Norm.coeff_data:
            coeffs, cols, timesteps, pops = Plot_Norm.coeff_data[Dfilename]
            norms = np.sum(pops, axis=1)
            # Convert the AU to fs (if required)
            fit = np.polyfit(timesteps * unit_conv, norms, 1)
            errs = [1e-10]+[1e-6]*(len(timesteps)-1)
            fit2, pcov = curve_fit(linear_fit,
                                   timesteps,
                                   norms,
                                   p0=[fit[0], 1],
                                   sigma=errs)
            # Convert per fs to per ps
            self.norm_drift_per_rep.append(fit2[0]*1000)

        self.norm_drift_per_rep = np.array(self.norm_drift_per_rep)
        self.worst_reps['norm'] = np.argmax(
                                           np.absolute(self.norm_drift_per_rep)
                                           ) + 1
        self.best_reps['norm'] = np.argmin(
                                           np.absolute(self.norm_drift_per_rep)
                                          ) + 1

    @staticmethod
    def put_drift_annotation(self):
        """
        Will put an annotation of the drift value on the norm graph
        """
        unit_conv = 1 
        if self.units == "au":
             unit_conv = 0.02418884

        coeffs, cols, timesteps, pops = Plot_Norm.avg_coeff_data
        norms = np.sum(pops, axis=1)
        fit = np.polyfit(timesteps * unit_conv, norms, 1)
        errs = [1e-10]+[1e-6]*(len(timesteps)-1)
        fit2, pcov = curve_fit(linear_fit,
                               timesteps * unit_conv,
                               norms,
                               p0=[fit[0], 1],
                               sigma=errs)
        text = r"Avg drift per rep = %.2g ps$^{-1}$" % (fit2[0] * 1000)

        if self.plot:
            # y1 = np.polyval(fit2, timesteps)
            # self.plot_ax.plot(timesteps, y1, 'k--', lw=0.5)
            all_norms = [np.sum(Plot_Norm.coeff_data[i][3], axis=1)
                         for i in Plot_Norm.coeff_data]
            min_len = np.min([len(i) for i in all_norms])
            all_norms = [i[:min_len] for i in all_norms]
            timesteps = timesteps[:min_len]

            min_x, min_y = min(timesteps), np.min(all_norms)
            range_x = max(timesteps) - min(timesteps)
            range_y = np.max(all_norms) - np.min(all_norms)
            loc = (min_x + 0.05 * range_x,
                   min_y + 0.85 * range_y)

            self.plot_ax.annotate(text, loc, fontsize=18)

        self.norm_drift = fit2[0]*1000
        self.plot_info['Norm'].append(text[:text.find('ps')])

    # Will plot normal of diabatic coeffs
    def _plot_norm_graph(self):
        """
        Will plot the normal of the diabatic coefficients. Will sum the
        populations and plot on the norm axis.
        """
        coeffs, cols, timesteps, pops = Plot_Norm.avg_coeff_data
        norms = np.sum(pops, axis=1)
        self.avg_line_norm, = self.plot_ax.plot(timesteps,
                                                norms,
                                                lw=2,
                                                color='g')

        self.avg_line_norm.set_visible(self.avg_reps_norm)
