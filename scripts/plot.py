#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 11:41:22 2018

@author: mellis
"""
import os; os.chdir('..')
from multiprocessing import Pool

from IO import Folders as fold
from PLOT import Plot

###############
# CTMQC_low_coup_2mol
folders = ['',
#            '/scratch/mellis/flavoured-cptk/DiffScalings/Coup=10meV',
           "/scratch/mellis/flavoured-cptk/200Rep_2mol/",
#           "/scratch/mellis/flavoured-cptk/Spare/NOTANH",
           '',
           ]

root_folders = ['',
                '/scratch/mellis/CTMQC_Detailed_Balance/test/0.05/',
                '']


plotting_parameters = ['norm', 'qm_t', 'rlk']
replicas = 'all'
plot = True
#######################################################

folders = [fold.make_fold_abs(i) for i in folders if os.path.isdir(i)]
folders = [i for i in folders if os.path.isfile(i+'run.inp')]

if len(folders) > 7 and plot:
    msg = "It seems you are trying to plot %i images!"
    msg += "\n\n\nThat is a lot of figures to load."
    msg += "These would probably be better being saved as "
    msg += "png files." % len(folders)
    raise SystemExit(msg)


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


if plot:
    all_p = []
    for i, f in enumerate(folders):
        all_p.append(do_1_fold_PL(f))
else:
    num_proc = len(folders)
    if num_proc > 21:
        num_proc = 21
    print("Num Proc = ", num_proc)
    if num_proc > 1:
        p = Pool(num_proc)
        if __name__ == '__main__':
            all_p = p.map(do_1_fold_PL, folders)
    else:
        all_p = [do_1_fold_PL(folders[0])]
