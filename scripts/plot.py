#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Run this from the root dir folder (the one with the load and IO folders in)
Created on Mon Oct 29 11:41:22 2018

@author: mellis

"""
import os
from multiprocessing import Pool

from IO import Folders as fold
from PLOT import Plot

###############
folders = ["/scratch/mellis/flavoured-cptk/200Rep_2mol"]

# folders = folders[:1]
plotting_parameters = ['fl']
replicas = 'all'
plot = False
num_proc = 'auto'
#######################################################

folders = [fold.make_fold_abs(i) for i in folders if os.path.isdir(i)]
folders = [i for i in folders if os.path.isfile(i+'run.inp')]


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
#        plt.close()
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
