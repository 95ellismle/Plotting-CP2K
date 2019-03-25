#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 11:41:22 2018

@author: mellis



TODO:
    * make the code work for new quantum momentum (without lk)
      This will involve translating the input parameter (qm_t) to
      either qlk or qm0 in the PLOT class. Once this is done I will
      just need to write a new plotting class/loading function for
      the qm0 case.
"""
import os; os.chdir('..')
from multiprocessing import Pool

from IO import Folders as fold
from PLOT import Plot

###############
# CTMQC_low_coup_2mol
# folders = ['',
#           '/scratch/mellis/flavoured-cptk/DiffScalings/Coup=10meV',
#           "/scratch/mellis/flavoured-cptk/200Rep_2mol",
#           '',
#           ]

root_fold = '/scratch/mellis/surface_hop/scripts-templates-for-aom-fssh/' + \
            "GENERATOR_FSSH_OS/TANH_WIDTH_CONV300K_Proper_1/" + \
            "Scal=0.0003/TANH_WIDTH=0.0001"
root_fold = "/scratch/mellis/flavoured-cptk/200Rep_2mol/"
folders = []
for dname, dirs, files in os.walk(root_fold):
    if 'run.inp' in files:
        folders.append(dname)
# folders = folders[:1]
plotting_parameters = ['norm', 'qm_t']
replicas = 'all'
plot = True
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

if replicas == 'all':
    print "Worst Reps = ", all_p[0].worst_reps
    print "Best Reps = ", all_p[0].best_reps

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