#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 12 18:28:54 2018

@author: mellis
"""
import os; os.chdir('..')
from multiprocessing import Pool
import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl

from IO import Folders as fold
from PLOT import Plot

# Data often too large to store as a pickle so this won't be used very much!
# Use the pickled object instead of expensive reading!
dumpFile = "AllTanhWidthData.pkl"
if os.path.isfile(dumpFile):
    print("\n\nUsing dump file!\n\n")
    print("If you don't want this change the dump file's name!\n\n\n")
    with open(dumpFile, 'r') as f:
        all_p = pkl.load(f)

# If there is no pickled object you'll have to read data
else:
    root_fold = '/scratch/mellis/surface_hop/' + \
                "scripts-templates-for-aom-fssh/GENERATOR_FSSH_OS/" + \
                "TANH_WIDTH_CONV300K_Proper"
    folders = []
    for dname, dirs, files in os.walk(root_fold):
        if 'run.inp' in files:
            folders.append(dname)

    plotting_parameters = ['energy_cons', '|C|^2', 'norm']
    replicas = 'all'
    plot = False
    num_proc = 'auto'
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
        if num_proc > 11:
            num_proc = 11
        print("Num Proc = ", num_proc)
        if num_proc > 1:
            pool = Pool(num_proc)
            if __name__ == '__main__':
                all_p = pool.map(do_1_fold_PL, folders)
        else:
            all_p = [do_1_fold_PL(folders[0])]

# Now get all scalings, couplings, numSteps etc...
allScal = list(set([p.run_inp_params['SCALING_FACTOR'] for p in all_p]))
couplings = {i: [] for i in allScal}
Edrifts = {i: [] for i in allScal}
norms = {i: [] for i in allScal}
tanh_widths = {i: [] for i in allScal}
numSteps = {i: [] for i in allScal}
for p in all_p:
    scal = p.run_inp_params['SCALING_FACTOR']
    tanh_width = p.run_inp_params['TANH_WIDTH']
    norm_drift = p.norm_drift
    ener_drift = np.mean(p.ener_drift_per_rep['tot'])
    tanh_widths[scal].append(tanh_width)
    norms[scal].append(norm_drift)
    Edrifts[scal].append(ener_drift)
    couplings[scal].append(p.coupling)
    numSteps[scal].append(p.num_di_coeff_steps)

toPlot = 'norm'

xlim = np.max([tanh_widths[i] for i in tanh_widths])
f, a = plt.subplots(3)
rgb = ['r', 'g', 'b']
for i, key in enumerate([0.0003, 0.003, 0.03]):
    window = 3
    x = tanh_widths[key]
    if toPlot == 'ener':
        y = np.absolute(Edrifts[key])
        ylab = r"Energy Drift [$\frac{Ha}{atom \ ps}$]"
    elif toPlot == 'norm':
        y = np.absolute(norms[key])
        ylab = r"Norm  Drift [ps$^{-1}$]"
    elif toPlot == 'steps':
        y = numSteps[key]
        ylab = r"Num Steps"

    xy = np.array(sorted(zip(x, y)))
    x = xy[:, 0]
    y = xy[:, 1]

    couplings[key] = np.array(couplings[key])
    couplings[key] = couplings[key].astype(float)
    currCoup = np.mean(couplings[key])
    a[i].plot(x, y, 'o', label=r"Coupling = %.2f meV" % (currCoup),
              color=rgb[i])
    a[i].set_title("")
    a[i].set_yscale("log")

    ytxt = np.max(y) - (np.max(y) - np.min(y))*0.2
    xtxt = xlim*0.7
    a[i].annotate(r"Coupling = %.2f meV" % (np.mean(couplings[key])),
                  (xtxt, ytxt), fontsize=24)

    if i < 2:
        a[i].tick_params(axis='x', labelbottom=False)

a[2].set_xlabel(r"Tanh Width [au$_f$]")
a[1].set_ylabel(ylab)

plt.subplots_adjust(top=0.945,
                    bottom=0.095,
                    left=0.075,
                    right=0.95,
                    hspace=0.134,
                    wspace=1.0)

save_p = all_p[:]
# Data too large to store as a pickle!
## Save the classes as a pickle
#with open(dumpFile, 'w') as f:
#    pkl.dump(all_p, f)
