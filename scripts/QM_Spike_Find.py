#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 15:02:02 2019

@author: oem
"""

import os; os.chdir('..')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

import PLOT 
from load import load_QM

folder = '/home/oem/Data/ToPlot/VarSigma/'

load = PLOT.LoadData(folder, [1], ['qm_t'])

qlk_data = load_QM.find_in_Qlk(load.all_Qlk_data['run-QM-1.xyz'][0], 
                               {'at_num':2, 'cart_dim':3})
qlk_data = qlk_data[:,0]
qlk_timesteps = load.all_Qlk_data['run-QM-1.xyz'][1]


dt = np.mean(np.diff(qlk_timesteps))
grad_data = np.gradient(qlk_data, dt)
mean_val = np.abs(np.mean(qlk_data))
std_var  = np.abs(np.std(qlk_data))


bad_times = []
for qlki in range(len(qlk_data)):
    if np.abs(qlk_data[qlki]) > mean_val + 2*std_var:
        bad_times.append(qlk_timesteps[qlki])


rects = []
for x in bad_times:
    rects.append(Rectangle((0+x,-1), 0.1, 2))

pc = PatchCollection(rects, facecolor='k', edgecolor=None, alpha=0.1)



# Create figure and axes
fig, ax = plt.subplots(1)
ax.plot(qlk_timesteps, qlk_data)
ax.plot(qlk_timesteps, qlk_data*grad_data)
ylims = ax.get_ylim()

ax.add_collection(pc)
ax.set_ylim(ylims)

plt.show()