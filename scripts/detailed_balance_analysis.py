#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 14:26:33 2019

@author: mellis
"""

import os; os.chdir("..") #need to change back to previous dir to link files properly
from multiprocessing import Pool

import PLOT as pl
from IO import Folders as fold


###############
#CTMQC_low_coup_2mol
folders = ['',
           '/scratch/mellis/CTMQC_Detailed_Balance/test/0.05/NREP=200/USEQM=T',
           '']
root_folders = ['',
           '/scratch/mellis/CTMQC_Detailed_Balance/test/0.05/',
           '']

#folders = []
#for root_folder in root_folders:
#    if root_folder and os.path.isdir(root_folder):
#        for folder, Dpaths, files in os.walk(root_folder):
#            if 'run.inp' in files:
#                folders.append(folder)
            

#folders              = ['/scratch/mellis/surface_hop/scripts-templates-for-aom-fssh/GENERATOR_FSSH_OS/TANH_WIDTH_CONV2_1/Scal=0.0003/TANH_WIDTH=8.5e-0',]# '/scratch/mellis/flavoured-cptk/200Rep_2mol_Ehren']
plotting_parameters = ['find_pos_crash']
replicas            = 'all'
plot                = True
#######################################################


folders = [fold.make_fold_abs(i) for i in folders if os.path.isdir(i)]
folders = [i for i in folders if os.path.isfile(i+'run.inp')]
if len(folders) > 7 and plot == True:
    raise SystemExit("It seems you are trying to plot %i images!\n\n\nThat is a lot of figures to load, these would probably be better being saved as png files."%len(folders))
if plot:
    all_p = []
    for i, f in enumerate(folders):
        all_p.append(pl.do_1_folder(f, plotting_parameters, replicas, plot))
else:
    num_proc = 4#len(folders)
    if num_proc > 21: num_proc = 21
    print("Num Proc = ", num_proc)
    if num_proc > 1:
        p = Pool(num_proc)
        if __name__ == '__main__':
            all_p = p.map(pl.do_1_folder, folders)
    else:
        all_p = [pl.do_1_folder(folders[0])]



#saveFolder = "/homes/mellis/Documents/Graphs/CTMQC/Detailed_Balance/ScalingGraphs/"
#if saveFolder[-1] != '/': saveFolder += '/'
#if not os.path.isdir(saveFolder): os.makedirs(saveFolder)
#
#for p in all_p:
#    mSaveFold = saveFolder + str(p.run_inp_params['SCALING_FACTOR']) + '/'
#    if not os.path.isdir(mSaveFold): os.makedirs(mSaveFold)
##    mSavePath = mSaveFold + 