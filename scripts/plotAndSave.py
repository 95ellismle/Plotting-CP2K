#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 11:02:31 2019

@author: mellis
"""

import os; os.chdir('..')
from multiprocessing import Pool
import matplotlib.pyplot as plt

from IO import Folders as fold
from PLOT import Plot

################
## CTMQC_low_coup_2mol
#root_folder = "/scratch/mellis/surface_hop/scripts-templates-for-aom-fssh/"
#root_folder += "GENERATOR_FSSH_OS/TANH_WIDTH_CONV2/"
#
#
#plotting_parameters = ['norm', 'energy_cons']
#replicas = 'all'
#plot = True
########################################################
#
#folders = []
#for dName, dirs, fNames in os.walk(root_folder):
#    if 'run.inp' in fNames:
#        folders.append(dName)
#
#folders = [fold.make_fold_abs(i) for i in folders if os.path.isdir(i)]
#folders = [i for i in folders if os.path.isfile(i+'run.inp')]
#
#
#def do_1_fold(folder):
#    p = Plot(plot_params=plotting_parameters,
#             folder=folder,
#             reps=replicas,
#             plot=plot)
#    return p
#
#
#all_p = []
#for i, f in enumerate(folders):
#    p = do_1_fold(f)
#    plt.close()
#    print(p.run_inp_params['TANH_WIDTH'])
#    all_p.append(p)
#

saveFolder = "/homes/mellis/Documents/Graphs/CTMQC/TestingTanh2/"
for p in all_p:
    saveFilePath = saveFolder[:]
    if saveFilePath[-1] != '/':
        saveFilePath += '/'
    saveFilePath += "Scal=%s/" % str(p.run_inp_params['SCALING_FACTOR'])
    if not os.path.isdir(saveFilePath):
        os.makedirs(saveFilePath)
    saveFilePath += "TanhWidth=%s.png" % str(p.run_inp_params['TANH_WIDTH'])
    print("Currently on %s" % saveFilePath)
    p.f.savefig(saveFilePath, dpi=300)
