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
              #"/scratch/mellis/flavoured-cptk/Timesteps/0.01NS_0.002ES",
              #"/scratch/mellis/flavoured-cptk/Timesteps/0.1NS_0.02ES",
              "/scratch/mellis/flavoured-cptk/NormCons/CTMQC",
              #"/scratch/mellis/flavoured-cptk/200Rep_2mol",
              #"/scratch/mellis/surface_hop/scripts-templates-for-aom-fssh/GENERATOR_FSSH_OS",
              "",
             ]


# folders = folders[:1]
plotting_parameters = ["ener_cons", "norm"]
replicas = 'all'
min_time = 0 
max_time = 1000
step = 10
savePath = "/homes/mellis/Documents/Graphs/CTMQC/New_QM/tmpImg"
iter_min_time = 0  # Where to start iterating if some pics already rendered.
#######################################################

folders = []
for rootfolder in rootFolder:
   for dpath, _, files in os.walk(rootfolder):
       if os.path.isdir(dpath):
         possFolder = os.path.abspath(dpath)
         if 'run.inp' in files:
            folders.append(possFolder)
            continue


if len(folders) != 1:
   print("folders = ", folders)
   raise SystemExit("This only works with 1 folder at the moment")

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


# Init params
p = do_1_folder(folder=folders[0],
                plotting_parameters=plotting_parameters,
                replicas=replicas,
                plot=True,
                minTime=min_time,
                maxTime=max_time
                )
ylims = []
for ax in p.axes:
   AX = p.axes[ax][1]
   if AX:
      ylims.append(AX.get_ylim())
gc.collect()

numDigits = len(str(max_time))

#all_p = []
count = 0
for maxT in range(iter_min_time, max_time, step):
    if maxT <= min_time:
       print("Skipping step %i" % maxT)
       continue
    p = do_1_folder(folder=folders[0],
                    plotting_parameters=plotting_parameters,
                    replicas=replicas,
                    plot=True,
                    minTime=min_time,
                    maxTime=maxT
                    )
    for axNum, ax in enumerate(p.axes):
      AX = p.axes[ax][1]
      if AX:
         AX.set_xlim([min_time, max_time])
         AX.set_ylim(ylims[axNum])

    fileName = str(count)
    fileName = "/" + "0" * (numDigits - len(fileName)) +fileName + ".png"
    p.f.savefig(savePath + fileName)
    gc.collect()
    count += 1
    print("Done %s" % fileName)

