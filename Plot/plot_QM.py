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

class QM(object):
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
        
        self._set_Qlk_control()
        
        self.qm_plot_ax.set_xlabel("R [au]")
        self.qm_plot_ax.set_ylabel(r"$\langle$ Q$_{lk, \nu}^{J}$ $\rangle_{J}$")
    
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
        """
        
        """
        1) Use masks to select Qlk and pos data that share same timesteps
        2) Find the pos and Qlk data corresponding to atom v
        3) Plot the pos (xyz) data vs Qlk (xyz) for each state. 
    
        """
        
        # First get shared timesteps between Qlk and pos
        self.avg_QM_data, avg_pos = self._match_timesteps_Qlk_pos(self.avg_Qlk_data, self.avg_pos_data)
        
        cart_dim = 1
        self.Qlk_ntimesteps = len(self.avg_QM_data['avg_Qlk'][0][0])
        print(self.Qlk_ntimesteps)
        natom = np.max(self.avg_Qlk_data['avg_Qlk'][0][1][0,:,0])
        self.Qlk_pos_avg_lines = [[],[],[]]
        
        cart_cols = ['r','g','b']
        
        all_QM_data = [[],[],[]]
        for cart_dim in range(3):
            for timestep in range(self.Qlk_ntimesteps):
                # Find all the QM_data for each atom
                QM_data = [load_QM.find_in_Qlk(self.avg_QM_data['avg_Qlk'][0], params={'at_num':at_num, 
                                                                'cart_dim':cart_dim,
                                                                'lk':(1,2)})[timestep][0]
                                    for at_num in range(natom)]
                line, = self.qm_plot_ax.plot(avg_pos['avg_pos'][0][0][timestep,:natom, cart_dim], 
                                             QM_data, 'o', color=cart_cols[cart_dim])
                self.Qlk_pos_avg_lines[cart_dim].append(line)
                all_QM_data[cart_dim].append(QM_data)
            
#        self.qm_plot_ax.set_ylim([np.min(all_QM_data), np.max(all_QM_data)])
#        self.qm_plot_ax.autoscale()#set_ylim([np.min(all_QM_data), np.max(all_QM_data)])
        
        for cart_dim in range(3):
            for line in self.Qlk_pos_avg_lines[cart_dim]:
                line.set_visible(False)
            if self.Qlk_cart_dims[cart_dim]:
                self.Qlk_pos_avg_lines[cart_dim][self.current_dt].set_visible(True)
        
    
    def _set_Qlk_control(self):
        """
        Will set the control panel for the Qlk graphs
        """
        self.MD_dt = self.run_inp_params['NUCLEAR_TIMESTEP']
        
        self.Qlk_cart_dim_butt  = CheckButtons(self.qm_widg_ax[0], ["X","Y","Z"], self.Qlk_cart_dims)
        self.Qlk_cart_dim_butt.on_clicked(self._Qlk_cart_dim_control)
        
        self.Qlk_dt_slider = Slider(self.qm_widg_ax[1], 'Time step', 0, self.Qlk_ntimesteps-1, valinit=0, valstep=1)
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
