#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 12:52:33 2018

@author: mellis
"""

import matplotlib.pyplot as plt
from matplotlib.widgets import CheckButtons


class Ad_Frc(object):
   """
   Will plot the adiabatic forces on each state.

    Inputs:
        axes  => a list of the axes to plot on. The first item should be the 
                 widget axis the second will be the axis to plot the data.
   """
   def __init__(self, axes):
      if self.plot:
         Ad_Frc.widg_ax = axes[0]
         Ad_Frc.plot_ax = axes[1]
         
         Ad_Frc.plot_all(self)

         Ad_Frc.plot_ax.set_ylabel(r"$-\nabla_{\nu} E_{l}^{ad, (I)}$ [$\frac{Ha}{bohr}$]")

   @staticmethod
   def plot_all(self):
      """
      Will plot a line for every replica.
      """
      for repKey in self.all_ad_frc_data:
         data, cols, timesteps = self.all_ad_frc_data[repKey]



class QM_Frc(object):
    """
    Will plot the Quantum Momentum forces against time.
    
    Inputs:
        axes  => a list of the axes to plot on. The first item should be the 
                 widget axis the second will be the axis to plot the data.
    """
    def __init__(self, axes):
        if self.plot: 
            QM_Frc.widget_ax = axes[0]
            QM_Frc.plot_ax = axes[1]
        
            #Setting initial default values
            QM_Frc.avg_reps = True
            QM_Frc.all_reps = False
            QM_Frc.atom = 6
        
            QM_Frc.all_rep_lines = []

            QM_Frc.plot_all(self)

            QM_Frc.plot_ax.set_ylabel(r"$|\mathbf{F}_{qm, \nu}^{(I)}|^2$ [au_f]")

    @staticmethod
    def plot_all(self):
      """
      Will plot the quantum momentum forces for each replica.
      """
      for key in self.all_qm_frc_data:
         data, cols, timesteps = self.all_qm_frc_data[key]
         carbonAtoms = [0, 3, 6, 9]
         hydrogenAtoms = [1, 2, 4, 5, 7, 8, 10, 11]
         for v in carbonAtoms:  # self.num_active_atoms:
            X = data[:, v, 0]  # X force
            Y = data[:, v, 1]  # Y force
            Z = data[:, v, 2]  # Z force
            mag = X**2 + Y**2 + Z**2

            QM_Frc.plot_ax.plot(timesteps,
                                mag, color='k',
                                lw=0.7, alpha=self.alpha)
         
         for v in hydrogenAtoms:  # self.num_active_atoms:
            X = data[:, v, 0]  # X force
            Y = data[:, v, 1]  # Y force
            Z = data[:, v, 2]  # Z force
            mag = X**2 + Y**2 + Z**2

            QM_Frc.plot_ax.plot(timesteps,
                                mag, color='y',
                                lw=0.7, alpha=self.alpha)




class Plot_Frc(object):
    """
    Will plot the nuclear force graph.
    
    Inputs:
        axes  => a list of the axes to plot on. The first item should be the 
                 widget axis the second will be the axis to plot the data.
    """
    def __init__(self, axes):
        if self.plot: 
            Plot_Frc.widget_ax = axes[0]
            Plot_Frc.plot_ax = axes[1]
        
            #Setting initial default values
            Plot_Frc.avg_reps = True
            Plot_Frc.all_reps = False
            Plot_Frc.atom = 6
        
            Plot_Frc.all_rep_lines = []
            Plot_Frc.avg_rep_lines = []
        
            #Plot_Frc.plot_avg_reps(self)
            Plot_Frc.plot_all_reps(self)
            
            

            #Connect checkboxes to plot control
            if self._use_control:    Plot_Frc._set_control()
            
            Plot_Frc.plot_ax.set_ylabel(r"|F$_{v}^{I}|^2$")

    #@staticmethod
    #def plot_avg_reps(self):
    #    """
    #    Will plot the average replica force
    #    """
    #    timesteps = self.avg_frc_data['avg_frc'][1]
    #    data = self.avg_frc_data['avg_frc'][0][0] #get all frc data
    #    for idim in range(3):
    #        Fdata = data[:,Plot_Frc.atom, idim]
    #        ln, = Plot_Frc.plot_ax.plot(timesteps, Fdata, lw=2)
    #        ln.set_visible(Plot_Frc.avg_rep_lines) #Decide inital visibility
    #        Plot_Frc.avg_rep_lines.append(ln)
    
    @staticmethod
    def plot_all_reps(self):
        """
        Will plot all the nuclear force data for each replica
        """
        for fname in self.all_frc_data:
            data, _, timesteps = self.all_frc_data[fname]
            carbonAtoms = [0, 3, 6, 9]
            hydrogenAtoms = [1, 2, 4, 5, 7, 8, 10, 11]
            for v in carbonAtoms:
               X = data[:, v, 0]  # X force
               Y = data[:, v, 1]  # Y force
               Z = data[:, v, 2]  # Z force
               mag = X**2 + Y**2 + Z**2

               Plot_Frc.plot_ax.plot(timesteps,
                                    mag, color='k',
                                    lw=0.7, alpha=self.alpha)
            
            for v in hydrogenAtoms:  # self.num_active_atoms:
               X = data[:, v, 0]  # X force
               Y = data[:, v, 1]  # Y force
               Z = data[:, v, 2]  # Z force
               mag = X**2 + Y**2 + Z**2

               Plot_Frc.plot_ax.plot(timesteps,
                                    mag, color='y',
                                    lw=0.7, alpha=self.alpha)


        #for idim in range(3):
        #    Fdata = data[:, Plot_Frc.atom, idim]
        #    ln, = Plot_Frc.plot_ax.plot(timesteps, Fdata, alpha=self.alpha, lw=1, color=self.colors[idim])
        #    ln.set_visible(Plot_Frc.all_rep_lines) #Decide inital visibility
        #    Plot_Frc.all_rep_lines.append(ln)
        
            
    @staticmethod
    def _set_control():
        """
        Will create the control panel for the force graph
        """
        Plot_Frc.check_site_ener = CheckButtons(Plot_Frc.widget_ax, ('all replicas', 'average'), (Plot_Frc.all_reps, Plot_Frc.avg_reps))
        Plot_Frc.check_site_ener.on_clicked(Plot_Frc._checkboxes_clicked)
    
    @staticmethod
    def _checkboxes_clicked(label):
        if 'all replicas' == label:
            for line in Plot_Frc.all_rep_lines:
                line.set_visible(not line.get_visible())
        if 'average' == label:
            for line in Plot_Frc.avg_rep_lines:
                line.set_visible(not line.get_visible())
        plt.draw()
  




# The following class is a template for a plot class, if a new thing needs 
#  plotting repeatedly (can use blank for quick temporary plots) then this 
#  class can be used as a template for the class
#class Plot_Frc(object):
#    """
#    Will plot the ... graph. 
#    
#    Inputs:
#        axes  => a list of the axes to plot on. The first item should be the 
#                 widget axis the second will be the axis to plot the data.
#    """
#    def __init__(self, axes):
#        if self.plot: 
#            Plot_Frc.widget_ax = axes[0]
#            Plot_Frc.plot_ax = axes[1]
#        
#            #Setting initial default values
#            Plot_Frc.avg_reps = True
#            Plot_Frc.all_reps = False
#            Plot_Frc.atom = 6
#        
#            Plot_Frc.all_rep_lines = []
#            Plot_Frc.avg_rep_lines = []
#        
#            Plot_Frc.plot_avg_reps(self)
#            Plot_Frc.plot_all_reps(self)
#            
#            
#
#            #Connect checkboxes to plot control
#            if self._use_control:    Plot_Frc._set_control()
#            
#            Plot_Frc.plot_ax.set_ylabel(r"F$_{%i}^{I}$"%(Plot_Frc.atom+1))
#
#    @staticmethod
#    def plot_avg_reps(self):
#        """
#        Will plot the average replica force
#        """
#        timesteps = self.avg_frc_data['avg_frc'][1]
#        data = self.avg_frc_data['avg_frc'][0][0] #get all frc data
#        for idim in range(3):
#            Fdata = data[:,Plot_Frc.atom, idim]
#            ln, = Plot_Frc.plot_ax.plot(timesteps, Fdata, lw=2)
#            ln.set_visible(Plot_Frc.avg_rep_lines) #Decide inital visibility
#            Plot_Frc.avg_rep_lines.append(ln)
#    
#    @staticmethod
#    def plot_all_reps(self):
#        """
#        Will plot all the nuclear force data for each replica
#        """
#        for fname in self.all_frc_data:
#            (data,_), timesteps = self.all_frc_data[fname]
#            for idim in range(3):
#                Fdata = data[:, Plot_Frc.atom, idim]
#                ln, = Plot_Frc.plot_ax.plot(timesteps, Fdata, alpha=self.alpha, lw=1, color=self.colors[idim])
#                ln.set_visible(Plot_Frc.all_rep_lines) #Decide inital visibility
#                Plot_Frc.all_rep_lines.append(ln)
#        
#            
#    @staticmethod
#    def _set_control():
#        """
#        Will create the control panel for the force graph
#        """
#        Plot_Frc.check_site_ener = CheckButtons(Plot_Frc.widget_ax, ('all replicas', 'average'), (Plot_Frc.all_reps, Plot_Frc.avg_reps))
#        Plot_Frc.check_site_ener.on_clicked(Plot_Frc._checkboxes_clicked)
#    
#    @staticmethod
#    def _checkboxes_clicked(label):
#        if 'all replicas' == label:
#            print("ALL")
#            for line in Plot_Frc.all_rep_lines:
#                line.set_visible(not line.get_visible())
#        if 'average' == label:
#            print("AVG")
#            for line in Plot_Frc.avg_rep_lines:
#                line.set_visible(not line.get_visible())
#        plt.draw()
