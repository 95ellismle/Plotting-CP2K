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

class Pos3D(object):
    """
    Will plot the positions in a plane with trajectories shown as trails.
    
    Inputs:
        axes => the 3D axes to plot on
    """
    def __init__(self, axis):
        if self.plot:
            Pos3D.ax = axis[0]
            
            Pos3D._plotAll(self)
            
            # Get rid of the ugly panes
            Pos3D.ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
            Pos3D.ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
            Pos3D.ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
            # Get rid of the ugly spines
            Pos3D.ax.w_xaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
            Pos3D.ax.w_yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
            Pos3D.ax.w_zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
            # Get rid of the ticks
            xticks = Pos3D.ax.get_xticks()
            yticks = Pos3D.ax.get_yticks()
            zticks = Pos3D.ax.get_zticks()

            xticks = np.linspace(xticks[1], xticks[-2], 2)
            yticks = np.linspace(yticks[1], yticks[-2], 2)
            zticks = np.linspace(zticks[1], zticks[-2], 2)

            Pos3D.ax.set_xticks(xticks)
            Pos3D.ax.set_yticks(yticks)
            Pos3D.ax.set_zticks(zticks)

            Pos3D.ax.set_xlabel(r"$Pos_{x}$ [bohr]")
            Pos3D.ax.set_ylabel(r"$Pos_{y}$ [bohr]")
            Pos3D.ax.set_zlabel(r"$Pos_{z}$ [bohr]")

    @staticmethod
    def _plotAll(self):
        """
        Will plot all positions the from each replica on a 3D axis. The traj is
        shown by a thin line.
        """
        for pKey in self.all_pos_data:
            pos, cols = self.all_pos_data[pKey][0]

            molAts = [np.arange(mol*6, (mol+1)*6)
                      for mol in range(self.num_states)]
            molCCols, molHCols = Pos3D.__get_mol_col(self)
            CAtoms, HAtoms = [], []
            for mol in molAts:
                CAtoms.append(mol[cols[0, mol] == 'C'])
                HAtoms.append(mol[cols[0, mol] == 'H'])

            for mol in range(self.num_states):
                # Plot last timestep
                for v in CAtoms[mol]:
                    Pos3D.ax.plot([pos[-1, v, 0]],
                                  [pos[-1, v, 1]],
                                  [pos[-1, v, 2]], 'o', alpha=0.6,
                                  color=molCCols[mol])
                for v in HAtoms[mol]:
                    Pos3D.ax.plot([pos[-1, v, 0]],
                                  [pos[-1, v, 1]],
                                  [pos[-1, v, 2]], 'o', alpha=0.6,
                                  color=molHCols[mol])
                # Plot Carbons
                for v in CAtoms[mol]:
                    Pos3D.ax.plot(pos[:, v, 0],
                                  pos[:, v, 1],
                                  pos[:, v, 2], '-', alpha=0.6,
                                  color=molCCols[mol])
                for v in HAtoms[mol]:
                    Pos3D.ax.plot(pos[:, v, 0],
                                  pos[:, v, 1],
                                  pos[:, v, 2], '-', alpha=0.6,
                                  color=molHCols[mol])

    @staticmethod
    def __get_mol_col(self):
        """
        Will return the colors of the atoms of each mol.
        """
        molCCols = [(i, i, i) for i in np.linspace(0, 0.4, self.num_states)]
        molHCols = [(i, i, 0) for i in np.linspace(1, 0.6, self.num_states)]
        return molCCols, molHCols


class PlotPos(object):
    """
    Will plot the positions vs time graph.

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

            axLab = r"$|\mathbf{R}_{\nu}^{(I)} - \mathbf{R}_{0}|$ [bohr]"
            PlotPos.plot_ax.set_ylabel(axLab)

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
            data, cols, timesteps = self.all_pos_data[fileName]

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


class PosStd(object):
    """
    Will plot the normalisation graph.

    Inputs:
        axes  => a list of the axes to plot on. The first item should be the
                 widget axis the second will be the axis to plot the data.
    """
    def __init__(self, axes):
        if self.plot:
            PosStd.widget_ax = axes[0]
            PosStd.plot_ax = axes[1]

            # Setting initial default values
            PosStd.__show_avg_reps = True
            PosStd.__show_all_reps = True

            # Will do the plotting
            PosStd.__plot_all(self)

            PosStd.plot_ax.set_ylabel(r"$\sigma(\mathbf{R})$")
#            PosStd.plot_ax.set_ylim([5, 1050])

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
        pKeys = list(self.all_pos_data.keys())
        mask = self.all_pos_data[pKeys[0]][1] != 'Ne'
        activeAtoms = [self.all_pos_data[fname][0][mask] for fname in self.all_pos_data]
        activeAtoms = np.reshape(activeAtoms, (self.num_reps,
                                               self.num_pos_steps,
                                               self.num_active_atoms,
                                               3)
                                 )
        
        self.allPos = activeAtoms

        timesteps = self.all_pos_data[pKeys[0]][2]

        distFromOrigin = np.linalg.norm(activeAtoms, axis=3)
        stdDevPerRepPerStep = np.std(distFromOrigin, axis=2)
        avgStdDev = np.mean(stdDevPerRepPerStep, axis=0)
        
        for stddev in stdDevPerRepPerStep:
            PosStd.plot_ax.plot(timesteps, stddev, lw=0.7)

        PosStd.plot_ax.plot(timesteps, avgStdDev, 'k--', lw=2)
            


    # Will plot normal of diabatic coeffs
    @staticmethod
    def __plot_avg(self):
        """
        Will plot the positions of atoms for all trajectories. Will sum the
        populations and plot on the norm axis.
        """
