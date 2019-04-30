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
import numpy as np

from IO import Folders as fold
from PLOT import Plot

###############
# Warning if root folder is set to a folder with other folders in it will crawl 
# the other folders in search of inputs to plot!
rootFolder = ["",
              #"/scratch/mellis/flavoured-cptk/200Rep_2molRenorm",
              "/scratch/mellis/flavoured-cptk/200Rep_2mol",
              #"/scratch/mellis/flavoured-cptk/LongQM",
              #"/scratch/mellis/flavoured-cptk/NormCons/Ehren",
              #"/scratch/mellis/flavoured-cptk/NormCons/CTMQC",
              "",
             ]
             
plotting_parameters = ['norm', "qm_t"]
replicas = 'all'
plot = True
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


def do_1_folder(folder, plotting_parameters, replicas, plot,
                minTime, maxTime):
    p = Plot(plot_params=plotting_parameters,
             folder=folder,
             reps=replicas,
             plot=plot,
             minTime=minTime,
             maxTime=maxTime,
             )
    return p

#allStats = []
all_p = []
for f in folders:
    p = do_1_folder(folder=f,
                    plotting_parameters=plotting_parameters,
                    replicas=replicas,
                    plot=plot,
                    minTime=min_time,
                    maxTime=max_time
                    )
    all_p.append(p)
    #allStats.append({
    #  'ener_cons': np.mean(p.ener_drift_per_rep['Total']),
    #  'norm': p.norm_drift,
    #  'scaling': p.run_inp_params['SCALING_FACTOR'],
    #  'NS': p.run_inp_params['NUCLEAR_TIMESTEP'],
    #  'ES': p.run_inp_params['NUCLEAR_TIMESTEP']//p.run_inp_params['ELECTRONIC_PARTIAL_STEP']
    #                })
    print("Done %s" % f)

plt.show()

