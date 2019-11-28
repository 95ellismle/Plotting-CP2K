#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 09:30:31 2018

@author: Sangeya
"""
import numpy as np
import matplotlib.cm as cm
from matplotlib.widgets import CheckButtons
import matplotlib.pyplot as plt


def linear_fit(x, m, c):
    """
    Will return a value on the line given by y = mx + c
    """
    return m * x + c


class Ad_State(object):
    """
    Will plot the adiabatic states and handle the checkbuttons turning on and
    off the 'fill between', 'all reps' and 'avg reps' options

    Inputs:
        * axes => The axes on which to plot (first is the checkbox axis,
                                             second is the plotting axis)
    """
    def __init__(self, axes):
        if self.plot:
            Ad_State.widg_ax, Ad_State.plot_ax = axes
            self.state_cols_AS = [i for i in self.all_ad_ener_data_avg.columns
                                  if 'State' in i]

            # Set defaults
            self.plot_all_AS = False
            self.plot_avg_AS = True
            self.fill_between_AS = False

            # Plot everything
            self._plot_all_reps_AS()
            self._plot_avg_rep_AS()
            self._fill_between_AS()

            # Set the control
            if self._use_control:
                self._set_control_AS()

            # Finish up (make things look pretty)
            Ad_State.plot_ax.set_ylabel(r"$E^{ad}_{l}$ [Ha]")

    def _plot_all_reps_AS(self):
        """
        Will plot all replicas as light thin lines on the graph
        """
        # Plot graphs
        self.AS_all_lines = []
        for Efilename in self.all_ad_ener_data:
            ad_ener_data = self.all_ad_ener_data[Efilename]
            for i, col in enumerate(self.state_cols_AS):
                ln, = Ad_State.plot_ax.plot(ad_ener_data['Time'],
                                           ad_ener_data[col],
                                           alpha=self.alpha,
                                           lw=0.7,
                                           color=self.colors[i])
                self.AS_all_lines.append(ln)
            Ad_State.plot_ax.plot(ad_ener_data['Time'],
                                  ad_ener_data['Pot'], 'k--',
                                  lw = 0.5, alpha=self.alpha)

        # Set initial visibility
        for line in self.AS_all_lines:
            line.set_visible(self.plot_all_AS)

    def _plot_avg_rep_AS(self):
        """
        Will plot the average replica as a thick line on the graph
        """
        # Plot the lines
        self.AS_avg_lines = []
        for i, col in enumerate(self.state_cols_AS):
            ln, = Ad_State.plot_ax.plot(self.all_ad_ener_data_avg['Time'],
                                       self.all_ad_ener_data_avg[col],
                                       color=self.colors[i])
            self.AS_avg_lines.append(ln)

        # Set initial visibility
        for line in self.AS_avg_lines:
            line.set_visible(self.plot_avg_AS)

    def _fill_between_AS(self):
        """
        Will fill the adiabatic states with a color dependant on how close
        they are to each other.
        """
        scaler = 2.5
        self.all_fill_bars = []
        # Fill colors originally
        for i, state in enumerate(self.state_cols_AS[:-1]):
            y1_y2 = self.all_ad_ener_data_avg[state] - \
                    self.all_ad_ener_data_avg[self.state_cols_AS[i+1]]
            max_diff, min_diff = np.max(y1_y2), np.min(y1_y2)
            normed_diffs = (y1_y2 - min_diff)/(max_diff - min_diff)

            for dt, timestep in enumerate(
                                         self.all_ad_ener_data_avg['Time'][:-1]
                                          ):
                T = [timestep, self.all_ad_ener_data_avg['Time'][dt+1]]
                y1 = self.all_ad_ener_data_avg[state].iloc[[dt, dt+1]]
                y2 = self.all_ad_ener_data_avg[
                                               self.state_cols_AS[i+1]
                                              ].iloc[[dt, dt+1]]
                diffy = scaler*(1-normed_diffs[dt])
                if diffy < 1:
                    color = cm.hot(diffy)

                    ln = Ad_State.plot_ax.fill_between(T,
                                                      y1,
                                                      y2,
                                                      alpha=0.4,
                                                      color=color,
                                                      lw=0)
                    self.all_fill_bars.append(ln)

        # Set initial visibility
        for line in self.all_fill_bars:
            line.set_visible(self.fill_between_AS)

    def _set_control_AS(self):
        """
        Will connect the checkbutton control panel with the plotting axis.
        """
        self.check_AS = CheckButtons(Ad_State.widg_ax,
                                     ('all rep', 'avg', 'fill'),
                                     (self.plot_all_AS,
                                      self.plot_avg_AS,
                                      self.fill_between_AS))
        self.check_AS.on_clicked(self._on_click_AS)

    def _on_click_AS(self, label):
        """
        Will handle button presses from the check buttons
        """
        if label == 'all rep':
            for line in self.AS_all_lines:
                line.set_visible(not line.get_visible())
        elif label == 'avg':
            for line in self.AS_avg_lines:
                line.set_visible(not line.get_visible())
        elif label == 'fill':
            for line in self.all_fill_bars:
                line.set_visible(not line.get_visible())
        plt.draw()


class Energy_Cons(object):
    """
    Will plot the kinetic, potential and total energy of the system.

    Inputs:
        * axes => The axes on which to plot (first is the checkbox axis,
                                             second is the plotting axis)
    """

    def __init__(self, axes):
        Energy_Cons.__calc_avg_ener_drift(self)
        if self.plot:
            Energy_Cons.widg_ax, Energy_Cons.plot_ax = axes

            Energy_Cons.plot_all = True
            Energy_Cons.plot_avg = True
            Energy_Cons.kinetic = True
            Energy_Cons.potential = True
            Energy_Cons.total = True

            Energy_Cons.plot_all_energy_drifts(self)
            Energy_Cons.plot_avg_energy_drift(self)

            # Set the control panelfrom scipy.optimize import curve_fit
            if self._use_control:
                Energy_Cons.__set_control()

            Energy_Cons.plot_ax.set_ylabel("E [Ha]")
            Energy_Cons.__put_avg_ener_drift(self)

    @staticmethod
    def plot_all_energy_drifts(self):
        """
        Will plot the all energy drifts on the Energy_Cons.plot_ax
        """
        Energy_Cons.kin_lines = []
        Energy_Cons.pot_lines = []
        Energy_Cons.tot_lines = []
        # Total energy
        for irep in self.all_tot_ener:
            data = self.all_tot_ener[irep]
            alpha = self.alpha + 0.1
            if alpha > 1:
                alpha = 1
            ln, = Energy_Cons.plot_ax.plot(data['Time'],
                                           data['E_cons'],
                                           'k-',
                                           lw=0.9,
                                           alpha=alpha)
            Energy_Cons.tot_lines.append(ln)
            ln.set_visible(Energy_Cons.plot_all and Energy_Cons.total)

            ln, = Energy_Cons.plot_ax.plot(data['Time'],
                                           data['Kin'],
                                           'r-',
                                           lw=0.9,
                                           alpha=alpha)
            Energy_Cons.kin_lines.append(ln)
            ln.set_visible(Energy_Cons.plot_all and Energy_Cons.kinetic)

            ln, = Energy_Cons.plot_ax.plot(data['Time'],
                                           data['Pot'],
                                           'g-',
                                           lw=0.9,
                                           alpha=alpha)
            Energy_Cons.pot_lines.append(ln)
            ln.set_visible(Energy_Cons.plot_all and Energy_Cons.potential)

    @staticmethod
    def plot_avg_energy_drift(self):
        """
        Will plot the average energy drift on the Energy_Cons.plot_ax
        """
        Energy_Cons.avg_lines = {}
        ln, = Energy_Cons.plot_ax.plot(self.tot_ener_mean['Time'],
                                       self.tot_ener_mean['E_cons'],
                                       'k-',
                                       lw=1.3)
        Energy_Cons.avg_lines['Total'] = ln
        ln.set_visible(Energy_Cons.plot_all and Energy_Cons.total)

        ln, = Energy_Cons.plot_ax.plot(self.tot_ener_mean['Time'],
                                       self.tot_ener_mean['Kin'],
                                       'r-',
                                       lw=1.3)
        Energy_Cons.avg_lines['Kinetic'] = ln
        ln.set_visible(Energy_Cons.plot_all and Energy_Cons.kinetic)

        ln, = Energy_Cons.plot_ax.plot(self.tot_ener_mean['Time'],
                                       self.tot_ener_mean['Pot'],
                                       'g-',
                                       lw=1.3)
        Energy_Cons.avg_lines['Potential'] = ln
        ln.set_visible(Energy_Cons.plot_all and Energy_Cons.potential)

        Energy_Cons.__calc_avg_ener_drift(self)

        plt.draw()

    @staticmethod
    def __set_control():
        """
        Will connect the checkbutton control panel with the plotting axis.
        """
        Energy_Cons.checkButs = CheckButtons(Energy_Cons.widg_ax,
                                             ('all reps',
                                              'avg',
                                              'Kinetic',
                                              'Potential',
                                              'Total Energy'),
                                             (Energy_Cons.plot_all,
                                              Energy_Cons.plot_avg,
                                              Energy_Cons.kinetic,
                                              Energy_Cons.potential,
                                              Energy_Cons.total))
        Energy_Cons.checkButs.on_clicked(Energy_Cons.__on_click)

    @staticmethod
    def __on_click(label):
        if label == 'all reps':
            Energy_Cons.plot_all = not Energy_Cons.plot_all
        elif label == 'avg':
            Energy_Cons.plot_avg = not Energy_Cons.plot_avg
        elif label == 'Kinetic':
            Energy_Cons.kinetic = not Energy_Cons.kinetic
        elif label == 'Potential':
            Energy_Cons.potential = not Energy_Cons.potential
        elif label == 'Total Energy':
            Energy_Cons.total = not Energy_Cons.total

        # Handle the all rep lines
        for kln, tln, pln in zip(Energy_Cons.kin_lines,
                                 Energy_Cons.tot_lines,
                                 Energy_Cons.pot_lines):
            kln.set_visible(Energy_Cons.plot_all and Energy_Cons.kinetic)
            tln.set_visible(Energy_Cons.plot_all and Energy_Cons.total)
            pln.set_visible(Energy_Cons.plot_all and Energy_Cons.potential)

        # Handle the average lines
        Energy_Cons.avg_lines['Total'].set_visible(Energy_Cons.plot_avg and
                                                 Energy_Cons.total)
        Energy_Cons.avg_lines['Kinetic'].set_visible(Energy_Cons.plot_avg and
                                                 Energy_Cons.kinetic)
        Energy_Cons.avg_lines['Potential'].set_visible(Energy_Cons.plot_avg and
                                                 Energy_Cons.potential)

        plt.draw()

    @staticmethod
    def __get_all_rep_drifts(self):
        """
        Will get the energy drift for each replica
        """
        lab_to_name_map = {'Kin': 'Kinetic', 'Pot': 'Potential', 'E_cons': 'Total'}
        
        unit_conv = 1
        if self.units == "au":
           unit_conv = 0.02418884

        # Find drifts
        self.ener_drift_per_rep = {'Kinetic': [], 'Potential': [], 'Total': []}
        for irep in self.all_tot_ener:
            data = self.all_tot_ener[irep]
            for lab in ('E_cons', 'Kin', 'Pot'):
                # Convert the time units from AU -> fs.
                fit = np.polyfit(data['Time']*unit_conv, data[lab], 1) 
                name = lab_to_name_map[lab]
                # Timestep in fs and Energy in Hartree so x1000 to convert to ps
                conv = (self.dt * 1000) / self.num_active_atoms
                self.ener_drift_per_rep[name].append(np.array(fit[0]) * conv)

        # Find largest and smallest drifts rep indices
        for lab in ('E_cons', 'Kin', 'Pot'):
            name = lab_to_name_map[lab]
            self.worst_reps[name] = np.argmax(self.ener_drift_per_rep[name])+1
            self.best_reps[name] = np.argmin(self.ener_drift_per_rep[name])+1

    @staticmethod
    def __calc_avg_ener_drift(self):
        """
        Will fit a linear line of best fit to the average 
        the average drift per ps.
        """
        Energy_Cons.__get_all_rep_drifts(self)

        self.Tot_Avg_Energy_Drift = np.mean(self.ener_drift_per_rep['Total'])

    @staticmethod
    def __put_avg_ener_drift(self):
        """
        Will put the average total energy drift on the plot as an annotation
        """
        minX = np.min(self.tot_ener_mean['Time'])
        minY = np.min(self.tot_ener_mean['E_cons'])
        allTotEner = [self.all_tot_ener[irep]['E_cons']
                      for irep in self.all_tot_ener]
        rangeX = np.max(self.tot_ener_mean['Time']) - minX
        rangeY = np.max([np.max(i) for i in allTotEner]) - minY

        x, y = minX + (rangeX * 0.1), minY + (rangeY * 0.9)
        e = self.Tot_Avg_Energy_Drift
        ann_txt = "Average Total Energy Drift = " + \
            r"%.2g $\frac{Ha}{ps \ atom}$" % e
        Energy_Cons.plot_ax.annotate(ann_txt,
                                     (x, y), fontsize=18)
