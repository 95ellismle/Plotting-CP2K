#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 14:26:33 2019

@author: mellis
"""

import os; os.chdir("..") #need to change back to previous dir to link files properly
import numpy as np 
#from multiprocessing import Pool

import PLOT as pl
#from IO import Folders as fold


###############
#CTMQC_low_coup_2mol
folders = ['',
           '/scratch/mellis/CTMQC_Detailed_Balance/test/0.05/NREP=70/USEQM=T_PRINTALL_0.0001',
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
replicas            = 'all'
plot                = True


def findSmallestDist1Rep(rep):
    """
    Will find the smallest distance between atoms in 1 replica
    """
    rep1 = load.all_pos_data[rep][0]
    nonNeonPos[rep] = [rep1[0][i][rep1[1][0,:,0] != 'Ne'] for i in range(len(rep1[0]))]
    minRepDist = 10000
    whichAtoms = [0,0]
    stepNum    = 0
    for istep, step in enumerate(nonNeonPos[rep]):
        for iat, at in enumerate(step):
            distances = np.linalg.norm(step - at, axis=1)
            minDist = np.min(distances[distances > 0])
            if minDist < minRepDist:
                minRepDist = minDist
                whichAtoms = [iat, list(distances).index(minDist)]
                stepNum = istep
    return [minRepDist, ';'.join([str(i) for i in whichAtoms]), stepNum]


load = pl.LoadData(folders[1], 'all', plot_params=['pos'], avg_on=True)
nonNeonPos = {}
repDists = []
for rep in load.all_pos_data:
    minRepDist, whichAtoms, stepNum = findSmallestDist1Rep(rep)   
    printAts = whichAtoms.split(';')
    print("Minimum distance in replica %s is %.4g.\nThis was between atoms %s and %s and occurred at step %i\n\n"%(rep, minRepDist, printAts[0], printAts[1], stepNum))
    repDists.append([rep, minRepDist, whichAtoms, stepNum])

            
repDists = np.array(repDists)
minRep = np.argmin(repDists[:,1].astype(float))

print("\n\n##########################################\n")
print("Minimum Rep Dist     %.4g"%float(repDists[minRep][1]))
print("From Replica         %s"%repDists[minRep][0])
print("Between atoms        %s and %s"%(repDists[minRep][2].split(';')[0], repDists[minRep][2].split(';')[1]))
print("At Step Number       %i"%(int(repDists[minRep][3])))
print("\n##########################################\n\n")
