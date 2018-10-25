#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 16:11:09 2018

@author: mellis
"""

from matplotlib.widgets import Slider, CheckButtons
import matplotlib.pyplot as plt
import numpy as np

from load import load_QM

class QM_R(object):
    """
    Will plot the Quantum Momentum.
            
    Inputs:
        axes  => a list of the axes to plot on. The first item should be the 
                 widget axis the second will be the axis to plot the data.  
    """
    def __init__(self, axes):
        self.qm_widg_ax, self.qm_plot_ax = axes
        
        self.current_dt = 0
        self.Qlk_cart_dims  = [True, True, True]
        
        self._plot_avg_qm_vs_pos()
        
#        self._set_Qlk_control()
        
        self.qm_plot_ax.set_xlabel("R [au]", fontsize=28)
        self.qm_plot_ax.set_ylabel(r"${Q^{J}_{12, \nu}}$",
                                   fontsize=28)
    
    def _match_timesteps_Qlk_pos(self, all_Qlk_data, all_pos_data):
        """
        Will find the Qlk data and the posistion data that share the same 
        timesteps.
        """
        for Qlk_rep, pos_rep in zip(all_Qlk_data, all_pos_data):
            Qlk_data = all_Qlk_data[Qlk_rep]
            pos_data = all_pos_data[pos_rep]
            
            Qlk_timesteps = Qlk_data[1]
            pos_timesteps = pos_data[1]
            
            # Find pos that are definitely in Qlk
            pos_mask = np.arange(len(pos_timesteps))
            pos_mask = pos_mask[[i in Qlk_timesteps for i in pos_timesteps]]
            pos_data = [  [pos_data[0][0][pos_mask],
                           pos_data[0][1][pos_mask]],
                        pos_timesteps[pos_mask]]
            
            # Find Qlk that are definitely in the pos data (that is definitely in Qlk data)
            Qlk_mask = np.arange(len(Qlk_timesteps))
            Qlk_mask = Qlk_mask[[i in pos_timesteps for i in Qlk_timesteps]]
            Qlk_data = [  [Qlk_data[0][0][Qlk_mask],
                           Qlk_data[0][1][Qlk_mask]],
                        Qlk_timesteps[Qlk_mask]]
    
            # Pass back to dictionary
            all_Qlk_data[Qlk_rep]  = Qlk_data
            all_pos_data[pos_rep] = pos_data
        
        return all_Qlk_data, all_pos_data
    
    def _plot_avg_qm_vs_pos(self):
        """
        Will plot the quantum momentum vs position.

        1) Use masks to select Qlk and pos data that share same timesteps
        2) Find the pos and Qlk data corresponding to atom v
        3) Plot the pos (xyz) data vs Qlk (xyz) for each state. 
        """
        
        # First get shared timesteps between Qlk and pos
        all_QM_data, all_pos = self._match_timesteps_Qlk_pos(self.all_Qlk_data, self.all_pos_data)
        
        Qkeys = list(self.all_Qlk_data.keys())
        cart_dim = 1
        self.Qlk_ntimesteps = len(self.all_Qlk_data[Qkeys[0]][0][0])
        natom = np.max(self.all_Qlk_data[Qkeys[0]][0][1][0,:,0])
        self.cart_lines = [[],[],[]]
        self.at_lines   = [[] for i in range(natom)]
        cart_cols = ['r','g','b']
        
        #Best thing to do here is store everything in a 1D list and create 3 
        # mappings. 1 to map lines to cartesian dimension and 1 to map the line
        # to the atom index and 1 to map the line to the timestep. These can be
        # used to control the plotting.
#        for Qrep, Rrep in zip(self.all_Qlk_data, all_pos):
#            for cart_dim in range(3):
##                for at_num in range(1,natom+1):
#                    at_num = 1
#                    # Find all the QM_data for each atom
#                    QM_data = load_QM.find_in_Qlk(all_QM_data[Qrep][0], 
#                                                  params  =  {'at_num':at_num, 
#                                                              'cart_dim':cart_dim+1,
#                                                              'lk':(1,2)})
#                    line, = self.qm_plot_ax.plot(all_pos[Rrep][0][0][:,at_num, cart_dim], 
#                                                 QM_data[:,0], 'o', color=self.colors[at_num])
#                    
#                    self.cart_lines[cart_dim].append(line)
#                    self.at_lines[at_num-1].append(line)
#            
        for irep, (Qrep, Rrep) in enumerate(zip(self.all_Qlk_data, all_pos)):
            print("Plotted rep %i"%irep)
            for cart_dim in range(3):
                for dt in range(self.Qlk_ntimesteps):
                    # Find all the QM_data for each atom
                    QM_data = [load_QM.find_in_Qlk(all_QM_data[Qrep][0], 
                                                  params  =  {'at_num':at_num,
                                                          'cart_dim':cart_dim+1,
                                                          'lk':(1,2)})[dt] 
                                                    for at_num in range(1,2)]
                    line, = self.qm_plot_ax.plot(all_pos[Rrep][0][0][dt,1, cart_dim], 
                                                 QM_data, 'o', color=self.colors[cart_dim])
                    
#                    self.cart_lines[cart_dim].append(line)
#                    self.at_lines[at_num-1].append(line)       
        
#        self.qm_plot_ax.set_ylim([np.min(all_QM_data), np.max(all_QM_data)])
#        self.qm_plot_ax.autoscale()#set_ylim([np.min(all_QM_data), np.max(all_QM_data)])
        
#        for cart_dim in range(3):
#            for line in self.Qlk_pos_avg_lines[cart_dim]:
#                line.set_visible(False)
#            if self.Qlk_cart_dims[cart_dim]:
#                self.Qlk_pos_avg_lines[cart_dim][self.current_dt].set_visible(True)
        
    
    def _set_Qlk_control(self):
        """
        Will set the control panel for the Qlk graphs
        """
        self.MD_dt = self.run_inp_params['NUCLEAR_TIMESTEP']
        
        self.Qlk_cart_dim_butt  = CheckButtons(self.qm_widg_ax[0], ["X","Y","Z"], self.Qlk_cart_dims)
        self.Qlk_cart_dim_butt.on_clicked(self._Qlk_cart_dim_control)
        
        self.Qlk_dt_slider = Slider(self.qm_widg_ax[1], 'Time step', 0, self.Qlk_ntimesteps-1, valinit=0)#, valstep=1)
        self.Qlk_dt_slider.on_changed(self._set_Qlk_slider_control)
        
        self.Qlk_vlines = [self.axes[param][1].axvline(self.current_dt*self.MD_dt) for param in self.non_qlk_params]
    
    def _Qlk_cart_dim_control(self, label):
        """
        Will control which cartesian dimensions are plotted. 
        (see plt.CheckButton docs)
        
        Inputs:
            *label   =>  The label of the checkbutton pressed
        """
        for i, dim in enumerate(["X","Y","Z"]):
            if label == dim: self.Qlk_cart_dims[i] = not self.Qlk_cart_dims[i]
        self._set_Qlk_slider_control(self.Qlk_dt_slider.val)
    
    def _set_Qlk_slider_control(self, val):
        """
        Will control what happens when the slider is interacted with
        """
        # Need to be an int as it indexes
        val = int(val)
        # Set the required timestep to be visible
        for cart_dim in range(3):
            self.Qlk_pos_avg_lines[cart_dim][self.current_dt].set_visible(False)
        for cart_dim, cart_val in enumerate(self.Qlk_cart_dims):
            self.Qlk_pos_avg_lines[cart_dim][val].set_visible(cart_val)
            self.current_dt = val
        
        # Draw the time lines
        for line in self.Qlk_vlines: line.remove()
        self.Qlk_vlines = [self.axes[param][1].axvline(val*self.MD_dt) for param in self.non_qlk_params]
        
        plt.figure(self.f.number)
        plt.draw()


class QM_t(object):
    """
    Will plot the Qlk against time graph.
    
    Inputs:
        axes  => a list of the axes to plot on. The first item_plot_site_ener should be the 
                 widget axis the second will be the axis to plot the data.
    """
    def __init__(self, axes):
        self.widget_ax = axes[0]
        self.plot_ax = axes[1]
        
        #Setting initial default values
        self.avg_reps_Qlk_t = True
        self.all_reps_Qlk_t = False
        
        #Plotting
        self._plot_all_rep_Qlk_t()
        self._plot_avg_QM_t()
        
        #Connect checkboxes to plot control
        if self._use_control:    self._set_Qlk_t_control()
        
        self.plot_ax.set_ylabel(r"$|Q_{12,\nu}^{(I)}|^{2}$")
    
    def _check_settings_Qlk_t(self, label):
        if label == 'all replicas': #Pressed the all replicas button
            for line in self.all_Qlk_t_lines:
                line.set_visible(not line.get_visible())
                
        elif label == 'average':
            for line in self.avg_qlk_t_lines:
                line.set_visible(not line.get_visible())
        plt.draw()
            
    #Will set the control panel for Qlk_t graph
    def _set_Qlk_t_control(self):
        self.check_Qlk_t = CheckButtons(self.widget_ax, ('all replicas', 'average'), (self.all_reps_Qlk_t, self.avg_reps_Qlk_t))
        self.check_Qlk_t.on_clicked(self._check_settings_Qlk_t)
        
    #Will plot the all replica Qlk_ts
    def _plot_all_rep_Qlk_t(self):
        """
        Will plot all the Qlk vs time lines for each replica
        """
        self.all_Qlk_t_lines = []
        
        ax = self.axes['qm_t'][1]
        for Qlk_filename in self.all_Qlk_data:
#                Qlk_filename = 'run-QM-1.xyz'
            Qlk_data = self.all_Qlk_data[Qlk_filename]
            Qlk_timesteps = Qlk_data[1]
            Qlk_data = Qlk_data[0]
            num_atoms = int(np.shape(Qlk_data[0])[1]/3)
            for iatom in range(1,num_atoms+1):
#                        if any(iatom == j for j in (2,4,11,8,1)): continue
                QMX = load_QM.find_in_Qlk(Qlk_data, params={'at_num':iatom, 'lk':(1,2), 'cart_dim':1})
                QMY = load_QM.find_in_Qlk(Qlk_data, params={'at_num':iatom, 'lk':(1,2), 'cart_dim':2})
                QMZ = load_QM.find_in_Qlk(Qlk_data, params={'at_num':iatom, 'lk':(1,2), 'cart_dim':3})
                
                QM_mag = np.sqrt(QMX**2 + QMY**2 + QMZ**2)
#                        QM_mag /= np.max(QM_mag)
                ln, = ax.plot(Qlk_timesteps, QMX, label="atom %i"%iatom, color=self.colors[iatom])
                self.all_Qlk_t_lines.append(ln)
            ax.set_ylabel(r"|Q$_{lk}^{(I)}$|$^2$")
            
        #Initialise the replica lines
        for line in self.all_Qlk_t_lines:
            line.set_visible(self.all_reps_Qlk_t)

    #Will plot the Quantum Momentum term
    def _plot_avg_QM_t(self):
        """
        Will plot the Qlk vs time for the average replica
        """
        self.avg_qlk_t_lines = []
        ax = self.axes['qm_t'][1]
        Qlk_filename = "avg_Qlk"
        Qlk_data = self.avg_Qlk_data[Qlk_filename]
        Qlk_timesteps = Qlk_data[1]
        Qlk_data = Qlk_data[0]
        num_atoms = int(np.shape(Qlk_data[0])[1]/3)
        for iatom in range(1,num_atoms+1):
            QMX = load_QM.find_in_Qlk(Qlk_data, params={'at_num':iatom, 'lk':(1,2), 'cart_dim':1})
            QMY = load_QM.find_in_Qlk(Qlk_data, params={'at_num':iatom, 'lk':(1,2), 'cart_dim':2})
            QMZ = load_QM.find_in_Qlk(Qlk_data, params={'at_num':iatom, 'lk':(1,2), 'cart_dim':3})
            QM_mag = np.sqrt(QMX**2 + QMY**2 + QMZ**2)
#                        QM_mag /= np.max(QM_mag)
            ln, = ax.plot(Qlk_timesteps, QMX, label="atom %i"%iatom, color=self.colors[iatom])
            self.avg_qlk_t_lines.append(ln)


