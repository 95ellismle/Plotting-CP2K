#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 16:11:09 2018

@author: mellis
"""

from matplotlib.widgets import Slider, CheckButtons, TextBox
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

        for irep, (Qrep, Rrep) in enumerate(zip(self.all_Qlk_data, all_pos)):
            print("Plotted rep %i"%irep)
            for cart_dim in range(3):
                for dt in range(self.Qlk_ntimesteps):
                    # Find all the QM_data for each atom
                    QM_data = [load_QM.find_in_Qlk(all_QM_data[Qrep][0], 
                                                  params  =  {'at_num':at_num,
                                                          'cart_dim':cart_dim+1,
                                                          'lk':(1,2)})[dt] 
                                                    for at_num in range(natom)]
                    line, = self.qm_plot_ax.plot(all_pos[Rrep][0][0][dt,1, cart_dim], 
                                                 QM_data, 'o', color=self.colors[cart_dim])
                    
#                    self.cart_lines[cart_dim].append(line)
#                    self.at_lines[at_num-1].append(line)       
        
#        all_QM_data = [[],[],[]]
#        for cart_dim in range(3):
#            for timestep in range(self.Qlk_ntimesteps):
#                # Find all the QM_data for each atom
#                QM_data = [load_QM.find_in_Qlk(self.avg_QM_data['avg_Qlk'][0], 
#                                               params   =      {'at_num':at_num, 
#                                                                'cart_dim':cart_dim+1,
#                                                                'lk':(1,2)})[timestep][0]
#                                    for at_num in range(1,natom+1)]
#    
#                line, = self.qm_plot_ax.plot(avg_pos['avg_pos'][0][0][timestep,:natom, cart_dim], 
#                                             QM_data, 'o', color=cart_cols[cart_dim])
#                self.Qlk_pos_avg_lines[cart_dim].append(line)
#                all_QM_data[cart_dim].append(QM_data)
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
        QM_t.widget_ax = axes[0]
        QM_t.plot_ax = axes[1]
        
        #Setting initial default values
        QM_t.X, QM_t.Y, QM_t.Z, QM_t.Mag = [False, False, False, True]
        
        #Plotting
#        self._plot_all_rep_Qlk_t()
        QM_t._plot_avg_QM_t(self)
        
        #Connect checkboxes to plot control
        if self._use_control:    QM_t._set_Qlk_t_control()
        
        QM_t.plot_ax.set_ylabel(r"$Q_{12,\nu}^{(I)}$")
    
    @staticmethod
    def _check_settings_Qlk_t(label):
        if label == 'X':
            QM_t.X = not QM_t.X
            for ln in QM_t.Xlines:
                ln.set_visible(QM_t.X)
        if label == 'Y':
            QM_t.Y = not QM_t.Y
            for ln in QM_t.Ylines:
                ln.set_visible(QM_t.Y)
        if label == 'Z':
            QM_t.Z = not QM_t.Z
            for ln in QM_t.Zlines:
                ln.set_visible(QM_t.Z)
        if label == 'Magnitude':
            QM_t.Mag = not QM_t.Mag
            for ln in QM_t.Maglines:
                ln.set_visible(QM_t.Mag)
        plt.draw()
    
    @staticmethod
    def _turn_on_atoms(Min, Max, type_switch='opposite'):
        """
        Will turn on atoms in a given range dependent on the current state of 
        the X, Y, Z and Mag checkboxes
        """
        poss_lines   = [QM_t.Xlines, QM_t.Ylines, QM_t.Zlines, QM_t.Maglines]
        all_switches = [QM_t.X, QM_t.Y, QM_t.Z, QM_t.Mag]
        for i, test in enumerate(all_switches):
            if test:
                for ln in poss_lines[i][ Min:Max ]:
                    if type_switch == 'opposite':
                        ln.set_visible(not ln.set_visible)
                    else:
                        ln.set_visible(type_switch)
        plt.draw()
        
    @staticmethod
    def _submit_text(text):
        """
        Will decide which atoms to plot from the submitted text in the text box
        """
        if text == 'all':
            QM_t._turn_on_atoms(0,len(QM_t.Xlines), True)
        else:
            text = text.split(',')
            for item in text:
                minmax = item.split('-')
                if len(minmax) == 2:
                    Min, Max = minmax
                    if Max == ':':
                        Max = len(QM_t.Xlines)
                    try:
                        Min = int(Min)
                        Max = int(Max)
                    except:
                        print("Tried to convert '%s' to ints but couldn't"%item)
                elif len(minmax) == 1:
                    try:
                        Min = int(minmax[0])
                        Max = int(minmax[0])+1
                    except:
                        print("Tried to convert '%s' to ints but couldn't"%item)
                if Max > len(QM_t.Xlines):
                    Max = len(QM_t.Xlines)
                QM_t._turn_on_atoms(0,len(QM_t.Xlines), False)
                QM_t._turn_on_atoms(Min, Max, True)
        
    #Will set the control panel for Qlk_t graph
    @staticmethod
    def _set_Qlk_t_control():
        QM_t.check_Qlk_t = CheckButtons(QM_t.widget_ax[0], ('X', 'Y', 'Z', 'Magnitude'), [QM_t.X, QM_t.Y, QM_t.Z, QM_t.Mag])
        QM_t.check_Qlk_t.on_clicked(QM_t._check_settings_Qlk_t)
        
        QM_t.textbox = TextBox(QM_t.widget_ax[1], 'Atoms:', initial="all")
        QM_t.textbox.on_submit(QM_t._submit_text)
        
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
                ln, = ax.plot(Qlk_timesteps, 
                              QM_mag, 
                              label="atom %i"%iatom, 
                              color=self.colors[iatom-1])
                self.all_Qlk_t_lines.append(ln)
            ax.set_ylabel(r"|Q$_{lk}^{(I)}$|$^2$")
            
        #Initialise the replica lines
        for line in self.all_Qlk_t_lines:
            line.set_visible(self.all_reps_Qlk_t)

    #Will plot the Quantum Momentum term
    @staticmethod
    def _plot_avg_QM_t(self):
        """
        Will plot the Qlk vs time for the average replica
        """
        ax = QM_t.plot_ax
        Qlk_filename = "avg_Qlk"
        Qlk_data = self.avg_Qlk_data[Qlk_filename]
        Qlk_timesteps = Qlk_data[1]
        Qlk_data = Qlk_data[0]
        num_atoms = int(np.shape(Qlk_data[0])[1]/3)
        
        QM_t.Xlines = []
        QM_t.Ylines = []
        QM_t.Zlines = []
        QM_t.Maglines = []
        
        for iatom in range(1,num_atoms+1):
            QMX = load_QM.find_in_Qlk(Qlk_data, params={'at_num':iatom, 'lk':(1,2), 'cart_dim':1})
            QMY = load_QM.find_in_Qlk(Qlk_data, params={'at_num':iatom, 'lk':(1,2), 'cart_dim':2})
            QMZ = load_QM.find_in_Qlk(Qlk_data, params={'at_num':iatom, 'lk':(1,2), 'cart_dim':3})
            QM_mag = np.sqrt(QMX**2 + QMY**2 + QMZ**2)
            
            ln, = ax.plot(Qlk_timesteps, 
                          QMX, 
                          color=self.colors[0])
            ln.set_visible(QM_t.X)
            QM_t.Xlines.append(ln)
            
            ln, = ax.plot(Qlk_timesteps, 
                          QMY, 
                          color=self.colors[1])
            ln.set_visible(QM_t.Y)
            QM_t.Ylines.append(ln)
            
            ln, = ax.plot(Qlk_timesteps, 
                          QMZ, 
                          color=self.colors[2])
            ln.set_visible(QM_t.Z)
            QM_t.Zlines.append(ln)
            
            ln, = ax.plot(Qlk_timesteps, 
                          QM_mag, 
                          '--', 
                          color=self.colors[iatom-1])
            ln.set_visible(QM_t.Mag)
            QM_t.Maglines.append(ln)
