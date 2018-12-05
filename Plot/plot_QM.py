#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 16:11:09 2018

@author: mellis
"""

from matplotlib.widgets import Slider, CheckButtons, TextBox
import matplotlib.pyplot as plt
import numpy as np
from collections import OrderedDict
import time

from load import load_QM

class QM_R(object):
    """
    Will plot the Quantum Momentum.
            
    Inputs:
        axes  => a list of the axes to plot on. The first item should be the 
                 widget axis the second will be the axis to plot the data.  
    """
    def __init__(self, axes):
        if self.plot:
            QM_R.widg_ax, QM_R.plot_ax = axes
            QM_R.timing_dict = OrderedDict()
            
            QM_R.lines = []
            QM_R.cart_dims  = [True, True, True]
            
            QM_R.timing_dict['Plot'] = time.time()
            QM_R._plot_qm_vs_pos_all_reps(self, 1000)
            QM_R.timing_dict['Plot'] = time.time() - QM_R.timing_dict['Plot']
            
            QM_R._set_Qlk_control(self)
            
            QM_R.plot_ax.set_xlabel("R [angstrom]", fontsize=28)
            QM_R.plot_ax.set_ylabel(r"${Q^{J}_{12, \nu}}$ [$\frac{Ha \cdot s}{l}$]",
                                       fontsize=28)
        
            self.print_timing_info(QM_R.timing_dict, "Qlk vs Pos timings")
    
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
    
    @staticmethod
    def put_time_annotation(self):
        """
        Will put an annotation of the current time on the graph
        """
        
    @staticmethod
    def _plot_qm_vs_pos_all_reps(self, istep):
        """
        Will plot the quantum momentum vs position.

        1) Use masks to select Qlk and pos data that share same timesteps
        2) Find the pos and Qlk data corresponding to atom v
        3) Plot the pos (xyz) data vs Qlk (xyz) for each state. 
        """
        
        # First get shared timesteps between Qlk and pos
        all_QM_data, all_pos = self._match_timesteps_Qlk_pos(self.all_Qlk_data, self.all_pos_data)
        
        Qkeys, Pkeys = list(all_QM_data.keys()), list(all_pos.keys())
        
        self.Qlk_ntimesteps = len(self.all_Qlk_data[Qkeys[0]][0][0])
        natom = np.max(self.all_Qlk_data[Qkeys[0]][0][1][0,:,0])
        self.cart_lines = [[],[],[]]
        self.at_lines   = [[] for i in range(natom)]

        idim   = 0
        iatom  = 0
        for irep in range(len(all_QM_data)):
            for iatom in range(natom):
                rep_key = Qkeys[irep]
                Qlk_data = load_QM.find_in_Qlk(all_QM_data[rep_key][0], {
                                                                     'at_num':iatom+1, 
                                                                     'cart_dim':idim+1, 
                                                                     'lk':(1,2)      })
            
                if self.run_inp_params['NUMBER_ATOMS_PER_SITE'] > 6:
                    print("""The QM vs pos graph won't work properly for anything but 
                          Ethylene. This is because the code isn't matching the correct
                          active atoms position to the quantum momentum. This can be 
                          done using the AOM_COEFF.include file. It has been 
                          implemented in the movie maker. There will be some useful 
                          code there to help with this.""")
                    raise SystemExit("UNSUPPORTED MOLECULE TYPE")
                else:
                    act_atoms = self.run_inp_params['NUMBER_ATOMS_PER_SITE'] * \
                                self.run_inp_params['NUMBER_DIABATIC_STATES']
                    
                single_rep_pos = all_pos[Pkeys[irep]][0][0]
                single_rep_pos = single_rep_pos[:,:act_atoms]
                
                ln, = QM_R.plot_ax.plot(single_rep_pos[istep,iatom, idim], 
                                  Qlk_data[istep,0], 
                                  '.', color=self.colors[iatom])
                QM_R.lines.append(ln)

        
        QM_R.put_time_annotation(self)
    
    @staticmethod
    def _set_qm_vs_pos_all_reps(self, istep):
        """
        Will reset the line ydata according to the timestep.
        """
        
        # First get shared timesteps between Qlk and pos
        all_QM_data, all_pos = self._match_timesteps_Qlk_pos(self.all_Qlk_data, self.all_pos_data)
        
        Qkeys, Pkeys = list(all_QM_data.keys()), list(all_pos.keys())
        
        self.Qlk_ntimesteps = len(self.all_Qlk_data[Qkeys[0]][0][0])
        natom = np.max(self.all_Qlk_data[Qkeys[0]][0][1][0,:,0])
        self.cart_lines = [[],[],[]]
        self.at_lines   = [[] for i in range(natom)]

        idim   = 0
        count = 0
        for irep in range(len(all_QM_data)):
            for iatom in range(natom):
                rep_key = Qkeys[irep]
                Qlk_data = load_QM.find_in_Qlk(all_QM_data[rep_key][0], {
                                                                     'at_num':iatom+1, 
                                                                     'cart_dim':idim+1, 
                                                                     'lk':(1,2)      })
            
                if self.run_inp_params['NUMBER_ATOMS_PER_SITE'] > 6:
                    print("""The QM vs pos graph won't work properly for anything but 
                          Ethylene. This is because the code isn't matching the correct
                          active atoms position to the quantum momentum. This can be 
                          done using the AOM_COEFF.include file. It has been 
                          implemented in the movie maker. There will be some useful 
                          code there to help with this.""")
                    raise SystemExit("UNSUPPORTED MOLECULE TYPE")
                else:
                    act_atoms = self.run_inp_params['NUMBER_ATOMS_PER_SITE'] * \
                                self.run_inp_params['NUMBER_DIABATIC_STATES']
                    
                single_rep_pos = all_pos[Pkeys[irep]][0][0]
                single_rep_pos = single_rep_pos[:,:act_atoms]
                
                QM_R.lines[count].set_ydata(Qlk_data[istep,0])
                QM_R.lines[count].set_xdata(single_rep_pos[istep,iatom, idim])
                count += 1
        plt.draw()
        QM_R.put_time_annotation(self)
        
    @staticmethod
    def _set_Qlk_control(self):
        """
        Will set the control panel for the Qlk graphs
        """
        self.MD_dt = self.run_inp_params['NUCLEAR_TIMESTEP']

        QM_R.Qlk_dt_slider = Slider(QM_R.widg_ax, 'Time step', 0, self.Qlk_ntimesteps-1, valinit=0, valstep=1)
        QM_R.SELF = self #Need this for the on_changed function... there is probably a better way though.
        QM_R.Qlk_dt_slider.on_changed(QM_R._set_Qlk_slider_control)
        
#        self.Qlk_vlines = [self.axes[param][1].axvline(QM_R.current_dt*self.MD_dt) for param in self.non_qlk_params]
    
#    def _Qlk_cart_dim_control(self, label):
#        """
#        Will control which cartesian dimensions are plotted. 
#        (see plt.CheckButton docs)
#        
#        Inputs:
#            *label   =>  The label of the checkbutton pressed
#        """
#        for i, dim in enumerate(["X","Y","Z"]):
#            if label == dim: QM_R.cart_dims[i] = not QM_R.cart_dims[i]
#        self._set_Qlk_slider_control(self.Qlk_dt_slider.val)
    
    @staticmethod
    def _set_Qlk_slider_control(val):
        """
        Will control what happens when the slider is interacted with
        """
        # Need to be an int as it indexes
        val = int(val)
        QM_R._set_qm_vs_pos_all_reps(QM_R.SELF, val)
        
#        # Draw the time lines
#        for line in self.Qlk_vlines: line.remove()
#        self.Qlk_vlines = [self.axes[param][1].axvline(val*self.MD_dt) for param in self.non_qlk_params]
#        
#        plt.figure(self.f.number)
#        plt.draw()



















class QM_t(object):
    """
    Will plot the Qlk against time graph.
    
    Inputs:
        axes  => a list of the axes to plot on. The first item_plot_site_ener should be the 
                 widget axis the second will be the axis to plot the data.
    """
    def __init__(self, axes):
        if self.plot:
            QM_t.widget_ax = axes[0]
            QM_t.plot_ax = axes[1]
            
            #Setting initial default values
            QM_t.X, QM_t.Y, QM_t.Z, QM_t.Mag = [False, False, False, True]
            
            #Plotting
    #        self._plot_all_rep_Qlk_t()
            QM_t._plot_avg_QM_t(self)
            
            #Connect checkboxes to plot control
            if self._use_control:    QM_t._set_Qlk_t_control()
            
            QM_t.plot_ax.set_ylabel(r"$Q_{12,\nu}^{(I)}$ [$\frac{Ha s}{l}$]")
    
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
            
#            print(np.mean(QM_mag))
#            print(np.std(QM_mag))
#            ax.axhline(np.mean(QM_mag))
            self.QM_mag = QM_mag
#            ax.axhline(np.mean(QM_mag)+np.std(QM_mag)*2, color=self.colors[iatom-1])
            ln, = ax.plot(Qlk_timesteps, 
                          QM_mag[:,0], 
                          '--', 
                          color=self.colors[iatom-1])
            ln.set_visible(QM_t.Mag)
            QM_t.Maglines.append(ln)







#    #Will plot the all replica Qlk_ts
#    def _plot_all_rep_Qlk_t(self):
#        """
#        Will plot all the Qlk vs time lines for each replica
#        """
#        self.all_Qlk_t_lines = []
#        
#        ax = self.axes['qm_t'][1]
#        for Qlk_filename in self.all_Qlk_data:
##                Qlk_filename = 'run-QM-1.xyz'
#            Qlk_data = self.all_Qlk_data[Qlk_filename]
#            Qlk_timesteps = Qlk_data[1]
#            Qlk_data = Qlk_data[0]
#            num_atoms = int(np.shape(Qlk_data[0])[1]/3)
#            for iatom in range(1,num_atoms+1):
##                        if any(iatom == j for j in (2,4,11,8,1)): continue
#                QMX = load_QM.find_in_Qlk(Qlk_data, params={'at_num':iatom, 'lk':(1,2), 'cart_dim':1})
#                QMY = load_QM.find_in_Qlk(Qlk_data, params={'at_num':iatom, 'lk':(1,2), 'cart_dim':2})
#                QMZ = load_QM.find_in_Qlk(Qlk_data, params={'at_num':iatom, 'lk':(1,2), 'cart_dim':3})
#                
#                QM_mag = np.sqrt(QMX**2 + QMY**2 + QMZ**2)
##                        QM_mag /= np.max(QM_mag)
#                ln, = ax.plot(Qlk_timesteps, 
#                              QM_mag, 
#                              label="atom %i"%iatom, 
#                              color=self.colors[iatom-1])
#                self.all_Qlk_t_lines.append(ln)
#            ax.set_ylabel(r"|Q$_{lk}^{(I)}$|$^2$")
#            
#        #Initialise the replica lines
#        for line in self.all_Qlk_t_lines:
#            line.set_visible(self.all_reps_Qlk_t)