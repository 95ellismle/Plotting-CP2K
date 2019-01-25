#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 13:31:22 2019

@author: mellis
"""

import os; os.chdir('..')
import matplotlib.pyplot as plt
import numpy as np

import PLOT

numReps = 200
folder = '/scratch/mellis/flavoured-cptk/DiffScalings/Coup=1meV/'

p = PLOT.Plot(plot_params=['energy_cons', 'norm', '|C|^2'],
              folder=folder,
              reps='all',
              plot=False)

allNormDrifts = p.norm_drift_per_rep
allEnerDrifts = p.ener_drift_per_rep['tot']

# Sort the drifts and find the best N ones
allSortedNormDrifts = sorted(allNormDrifts)[:numReps]
allSortedEnerDrifts = sorted(allEnerDrifts)[:numReps]

# Find the replica index of these
bestNormReps = [list(allNormDrifts).index(i) for i in allSortedNormDrifts]
bestEnerReps = [list(allEnerDrifts).index(i) for i in allSortedEnerDrifts]

# Find the coefficient data corresponding to those replicas
numStates = p.run_inp_params['NUMBER_DIABATIC_STATES']
timesteps = p.all_Acoeff_data['run-coeff_ad_1-1.xyz'][2]
numSteps = len(timesteps)
coeffData = np.zeros((numReps, numSteps, numStates))

for repi, rep in enumerate(bestEnerReps):
    coeffFileName = 'run-coeff_ad_%i-1.xyz' % (rep + 1)
    coeffData[repi] = p.all_Acoeff_data[coeffFileName][3]  # Populations
avgCoeffData = np.mean(coeffData, axis=0)

# Plot the best N reps
lastNumSteps = 1000000
f, a = plt.subplots()
for repData in coeffData:
    a.plot(timesteps[-lastNumSteps:],
           repData[:, 0][-lastNumSteps:],
           'b-',
           alpha=0.1,
           lw=0.5)
    a.plot(timesteps[-lastNumSteps:],
           repData[:, 1][-lastNumSteps:],
           'r-',
           alpha=0.1,
           lw=0.5)

# Plot Averages
a.plot(timesteps[-lastNumSteps:],
       avgCoeffData[:, 0][-lastNumSteps:],
       'b-',
       lw=1)
a.plot(timesteps[-lastNumSteps:],
       avgCoeffData[:, 1][-lastNumSteps:],
       'r-',
       lw=1)

totReps = p.run_inp_params['NREP']
a.set_title("Best %i/%i Replicas (norm cons)" % (numReps, totReps))
a.set_ylabel(r"|C|$^2$")
a.set_xlabel(r"Timestep [fs]")

plt.tight_layout()
