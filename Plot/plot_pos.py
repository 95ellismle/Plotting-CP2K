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


class COM(object):
    """
    Will plot the center of mass of the positions
     
     Inputs:
         axes => the 3D axes to plot on
    """
    def __init__(self, axis):
        if self.plot:
            COM_widg_ax, COM.plot_ax = axis
            
            COM._plotAll(self)
            COM.plot_ax.set_ylabel(r"|$\mathbf{COM}^{(I)} - \mathcal{O}$| [%s]" % self.unitsLength[self.units])
             
    @staticmethod
    def _plotAll(self):
        """
        Will plot every replicas COM 
        """
        for rep in self.all_COM:
            com, timesteps = self.all_COM[rep]
            comDistFromOrigin = np.linalg.norm(com, axis=1)
            COM.plot_ax.plot(timesteps, comDistFromOrigin,
                             alpha=self.alpha, lw=0.7)
 


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

            Pos3D.ax.set_xlabel(r"$Pos_{x}$ [%s]" % self.unitsLength[self.units])
            Pos3D.ax.set_ylabel(r"$Pos_{y}$ [%s]" % self.unitsLength[self.units])
            Pos3D.ax.set_zlabel(r"$Pos_{z}$ [%s]" % self.unitsLength[self.units])

    @staticmethod
    def _plotAll(self):
        """
        Will plot all positions the from each replica on a 3D axis. The traj is
        shown by a thin line.
        """
        for pKey in self.all_pos_data:
            pos, cols, _ = self.all_pos_data[pKey]

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

            axLab = r"$\mathbf{R}_{\nu}^{(I)}$ [%s]" % self.unitsLength[self.units]
            PlotPos.plot_ax.set_ylabel(axLab)

    @staticmethod
    def __button_control(self, label):
        """
        Will handle the button events
        """

    # Will set the control panel for positions
    @staticmethod
    def __set_control(self):
        """
        Will put the control panel on the axis and connect the off and on
        functions to it.
        """

    # Will plot all positions
    @staticmethod
    def __plot_all(self):
        """
        Will plot the positions of atoms for all trajectories.
        """
        for fileName in self.all_pos_data:
            data, cols, timesteps = self.all_pos_data[fileName]

            # Plot atoms
            for iat in self.atoms_to_plot:
                v = self.active_atoms[iat-1]
                x = data[:, v, 0]
                #y = data[:, v, 1]
                #z = data[:, v, 2]
                #mag = np.sqrt(x**2 + y**2 + z**2)
                #mag = x

                #PlotPos.plot_ax.plot(timesteps,
                #                     mag,
                #                     color=self.colors[iat],
                #                     #color=self.colors[iat],
                #                     alpha=self.alpha,
                #                     lw=0.5)
                PlotPos.plot_ax.plot(timesteps,
                                     x,
                                     color='r',  # self.colors[iat-1],
                                     alpha=self.alpha,
                                     lw=0.5)

                #PlotPos.plot_ax.plot(timesteps,
                #                     y,
                #                     color='g',  # self.colors[iat-1],
                #                     alpha=self.alpha,
                #                     lw=0.5)

                #PlotPos.plot_ax.plot(timesteps,
                #                     z,
                #                     color='b',  # self.colors[iat-1],
                #                     alpha=self.alpha,
                #                     lw=0.5)

    # Will plot normal of diabatic coeffs
    @staticmethod
    def __plot_avg(self):
        """
        Will plot the positions of atoms for all trajectories. Will sum the
        populations and plot on the norm axis.
        """


class PlotVel(object):
    """
    Will plot the velocities vs time graph.

    Inputs:
        axes  => a list of the axes to plot on. The first item should be the
                 widget axis the second will be the axis to plot the data.
    """
    def __init__(self, axes):
        if self.plot:
            PlotVel.widget_ax = axes[0]
            PlotVel.plot_ax = axes[1]

            # Setting initial default values
            PlotVel.__show_avg_reps = True
            PlotVel.__show_all_reps = True

            # Will do the plotting
            PlotVel.__plot_all(self)
            # Connect checkboxes to plot control
            if self._use_control:
                PlotVel.__set_control(self)

            axLab = r"$|\mathbf{R}_{\nu}^{(I)}|$ [%s]" % self.unitsLength[self.units]
            PlotVel.plot_ax.set_ylabel(axLab)

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
        Will plot the velocities of atoms for all trajectories. Will sum the
        populations and plot on the norm axis.
        """
        for fileName in self.all_vel_data:
            data, cols, timesteps = self.all_vel_data[fileName]

            # This is a horrible hack should be improved!
            #  * Make more robust way to get active atoms
            #  * Use numpy fancy indexing instead of list comprehension
            mask = cols != 'Ne'
            activeAtoms = [[vel for cond, vel in zip(mask[step], data[step])
                            if cond]
                           for step in range(len(data))]
            activeAtoms = np.array(activeAtoms)

            # Plot atoms
            for iat in self.atoms_to_plot:
                x = activeAtoms[:, iat-1, 0]
                #y = activeAtoms[:, iat-1, 1]
                #z = activeAtoms[:, iat-1, 2]
                mag = x #np.sqrt(x**2 + y**2 + z**2)

                PlotVel.plot_ax.plot(timesteps,
                                     mag,
                                     color=self.colors[iat-1],
                                     alpha=self.alpha,
                                     lw=0.5)

    # Will plot normal of diabatic coeffs
    @staticmethod
    def __plot_avg(self):
        """
        Will plot the velocities of atoms for all trajectories. Will sum the
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
            PosStd.__plot(self)

            PosStd.plot_ax.set_ylabel(r"$\sigma(|\mathbf{R}_{\nu}^{(I)}|)$ [%s]" % self.unitsLength[self.units])
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
    def __plot(self):
        """
        Will plot the positions of atoms for all trajectories. Will sum the
        populations and plot on the norm axis.
        """
        pKeys = list(self.all_pos_data.keys())
        allPos = [self.all_pos_data[fName][0] for fName in self.all_pos_data]
        allPos = [np.sqrt(p[:, :, 0]**2 + p[:, :, 1]**2 + p[:, :, 2]**2) for p in allPos]
        std = np.std(allPos, axis=0)

        timesteps = self.all_pos_data[pKeys[0]][2]
        for iatom in self.atoms_to_plot:
           PosStd.plot_ax.plot(timesteps, std[:, iatom-1], color=self.colors[iatom], lw=0.7)

            


    # Will plot normal of diabatic coeffs
    @staticmethod
    def __plot_avg(self):
        """
        Will plot the positions of atoms for all trajectories. Will sum the
        populations and plot on the norm axis.
        """
