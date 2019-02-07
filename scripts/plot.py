#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 11:41:22 2018

@author: mellis
"""
import os; os.chdir('..')
from multiprocessing import Pool
import numpy as np
import matplotlib.pyplot as plt

from IO import Folders as fold
from PLOT import Plot

###############
# CTMQC_low_coup_2mol
folders = ['',
#            '/scratch/mellis/flavoured-cptk/DiffScalings/Coup=10meV',
           "/scratch/mellis/flavoured-cptk/200Rep_2mol",
           '',
           ]

#root_fold = '/scratch/mellis/surface_hop/scripts-templates-for-aom-fssh/' + \
#            "GENERATOR_FSSH_OS/TANH_WIDTH_CONV300K_Proper"
#folders = []
#for dname, dirs, files in os.walk(root_fold):
#    if 'run.inp' in files:
#        folders.append(dname)

plotting_parameters = ['norm', 'energy_cons', 'Rlk']
replicas = 'all'
plot = True
num_proc = 5
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


if num_proc == 'auto':
    num_proc = len(folders)

if plot:
    all_p = []
    for i, f in enumerate(folders):
        all_p.append(do_1_fold_PL(f))
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

## Now get all scalings, couplings, numSteps etc...
#allScal = list(set([p.run_inp_params['SCALING_FACTOR'] for p in all_p]))
#coup = {i: [] for i in allScal}
#Edrifts = {i: [] for i in allScal}
#norms = {i: [] for i in allScal}
#tanh_widths = {i: [] for i in allScal}
#numSteps = {i: [] for i in allScal}
#for p in all_p:
#    scal = p.run_inp_params['SCALING_FACTOR']
#    tanh_width = p.run_inp_params['TANH_WIDTH']
#    norm_drift = p.norm_drift
#    ener_drift = np.mean(p.ener_drift_per_rep['tot'])
#    tanh_widths[scal].append(tanh_width)
#    norms[scal].append(norm_drift)
#    Edrifts[scal].append(ener_drift)
#    coup[scal].append(p.coupling)
#    numSteps[scal].append(p.num_di_coeff_steps)
#
#xlim = 2
#f, a = plt.subplots(3)
#rgb = ['r', 'g', 'b']
#for i, key in enumerate([0.0003, 0.003, 0.03]):
#    window = 3
#    x = tanh_widths[key]
#    y = np.absolute(norms[key])
#
#    xy = np.array(sorted(zip(x, y)))
#    x = xy[:, 0]
#    y = xy[:, 1]
#
#    coup[key] = np.array(coup[key])
#    coup[key] = coup[key].astype(float)
#    currCoup = np.mean(coup[key])
#    a[i].plot(x, y, 'o', label=r"Coupling = %.2f meV" % (currCoup),
#              color=rgb[i])
#    a[i].set_title("")
#
#    ytxt = np.max(y) - (np.max(y) - np.min(y))*0.9
#    xtxt = xlim*0.7
#    a[i].annotate(r"Coupling = %.2f meV" % (np.mean(coup[key])),
#                  (xtxt, ytxt),
#                  fontsize=24)
#
##    a[i].set_xlim([-0.001,xlim])
#    if i < 2:
#        a[i].tick_params(axis='x', labelbottom=False)
#
#a[2].set_xlabel(r"Tanh Width [au$_f$]")
#a[1].set_ylabel(r"Energy Drift [$\frac{Ha}{atom \ ps}$]")
#
#plt.subplots_adjust(top=0.945,
#                    bottom=0.095,
#                    left=0.075,
#                    right=0.95,
#                    hspace=0.134,
#                    wspace=1.0)
