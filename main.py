#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Run this from the root dir folder (the one with the load and IO folders in)
Created on Mon Oct 29 11:41:22 2018

@author: mellis

"""
import os
from multiprocessing import Pool
import matplotlib.pyplot as plt

from IO import Folders as fold
from PLOT import Plot

###############
# Warning if root folder is set to a folder with other folders in it will crawl 
# the other folders in search of inputs to plot!
rootFolder = ["",
              #"/scratch/mellis/flavoured-cptk/200Rep_2mol", 
              #"/scratch/mellis/flavoured-cptk/PopTransfer/CTMQCForce_CTMQCCoeff",
              #"/scratch/mellis/flavoured-cptk/PopTransfer/EhrenForce_EhrenCoeff",
              #"/scratch/mellis/flavoured-cptk/PopTransfer/CTMQCForce_EhrenCoeff",
              "/scratch/mellis/flavoured-cptk/PopTransfer/EhrenForce_CTMQCCoeff",
             ]

# folders = folders[:1]
plotting_parameters = ["|C|^2", "fl"]
replicas = [1] 
plot = True
num_proc = 'auto'
#######################################################

folders = []
for rootfolder in rootFolder:
   for dpath, _, files in os.walk(rootfolder):
       if os.path.isdir(dpath):
         possFolder = os.path.abspath(dpath)
         if 'run.inp' in files:
            folders.append(possFolder)
            continue

if not folders:
   print("\t\t#####################")
   print("\nSorry I can't find any folders to plot!")
   print("Make sure the run.inp file is in the required folder")


def do_1_folder(folder, plotting_parameters, replicas, plot):
    p = Plot(plot_params=plotting_parameters,
             folder=folder,
             reps=replicas,
             plot=plot)
    return p


def do_1_fold_PL(folder):
    p = Plot(plot_params=plotting_parameters,
             folder=folder,
             reps=replicas,
             plot=plot)
    return p

if num_proc == 'auto':
    num_proc = len(folders)

if plot:
    all_p = []
    for i, f in enumerate(folders):
        all_p.append(do_1_fold_PL(f))
        plt.close()
        print("Done %s" % f)
else:
    if num_proc > 21:
        num_proc = 21
    print("Num Proc = ", num_proc)
    if num_proc > 1:
        pool = Pool(num_proc)
        if __name__ == '__main__':
            all_p = pool.map(do_1_fold_PL, folders)
    else:
        all_p = [do_1_fold_PL(folders[0])]

#for pNum, p in enumerate(all_p):
#   folder = folders[pNum]
#   newSavePath = folder[folder.rfind('/')+1:]
#   saveFolder = "/homes/mellis/Documents/Graphs/CTMQC/New_QM/BugHunt_DiffTimesteps/RawDataGraphs"
#   newSavePath = "%s/%s_norm_enerDrift.png" % (saveFolder, newSavePath)
#   p.f.savefig(newSavePath, dpi=200)

#if replicas == 'all':
#    print ("Worst Reps = ", all_p[0].worst_reps)
#    print ("Best Reps = ", all_p[0].best_reps)

#plt.figure()
#for key in p.all_tot_ener:
#    print(key)
#    allTimes = []
#    allFits = []
#    for iStep in p.all_tot_ener[key]['Time'][2:]:
#        mask = p.all_tot_ener[key]['Time'] <= iStep
#        times = p.all_tot_ener[key]['Time'][mask]
#        eners = p.all_tot_ener[key]['E_cons'][mask]
#        
#        fit = np.polyfit(list(times), list(eners),1)
#        allTimes.append(iStep)
#        allFits.append(fit[0]*1000./12) 
#
#    plt.plot(allTimes, allFits)
