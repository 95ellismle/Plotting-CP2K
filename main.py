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
import gc

from IO import Folders as fold
from PLOT import Plot

###############
# Warning if root folder is set to a folder with other folders in it will crawl 
# the other folders in search of inputs to plot!
rootFolder = ["",
              #"/home/oem/Data/CTMQC/CTMQCAll",
              #"/home/oem/Data/CTMQC/EhrenAll",
              "/home/oem/Data/CTMQC/CTMQCForceEhrenCoeff",
              #"/home/oem/Data/CTMQC/EhrenForceCTMQCCoeff",
              "",
             ]

# folders = folders[:1]
plotting_parameters = ["energy_cons"]
replicas = [15] 
plot = True
num_proc = 'auto'
min_time = 0
max_time = 'all'
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
             plot=plot,
             minTime=min_time,
             maxTime=max_time,
             )
    return p


def do_1_fold_PL(folder):
    p = Plot(plot_params=plotting_parameters,
             folder=folder,
             reps=replicas,
             plot=plot,
             minTime=min_time,
             maxTime=max_time)
    return p

if num_proc == 'auto':
    num_proc = len(folders)

if plot:
    all_p = []
    for f in folders:
        p = do_1_fold_PL(f)
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

