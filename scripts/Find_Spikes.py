#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 13:02:29 2019

@author: mellis
"""

import os; os.chdir('..')
import matplotlib.pyplot as plt
import numpy as np

import PLOT

from load import load_QM

folder = '/scratch/mellis/flavoured-cptk/200Rep_2mol'

AllData = PLOT.LoadData(folder,
                        'all',
                        plot_params=['qm_t', 'rlk'],
                        avg_on=True)

appendTime = True
allAtRepData = []
allRepTimesteps = []
for iat in range(1, 12):

    allRepData = []
    for repNum in range(1, 201):
        if appendTime:
            allRepTimesteps.append(
                    AllData.all_Qlk_data['run-QM-%i.xyz' % repNum][1]
                                  )

        qlk_data = AllData.all_Qlk_data['run-QM-%i.xyz' % repNum][0]

        X = load_QM.find_in_Qlk(qlk_data, {'at_num': iat,
                                           'cart_dim': 1,
                                           'lk': (1, 2)
                                           }
                                )
        Y = load_QM.find_in_Qlk(qlk_data, {'at_num': iat,
                                           'cart_dim': 2,
                                           'lk': (1, 2)
                                           }
                                )
        Z = load_QM.find_in_Qlk(qlk_data, {'at_num': iat,
                                           'cart_dim': 3,
                                           'lk': (1, 2)
                                           }
                                )

        mag = np.sqrt(X**2 + Y**2 + Z**2)
        allRepData.append(mag)

    appendTime = False
    allAtRepData.append(np.array(allRepData))

allAtRepData = np.array(allAtRepData)

allAtRlkData = []
for iat in range(1, 12):

    rlk_data = AllData.Rlk_data[0]
    X = load_QM.find_in_Qlk(rlk_data, {'at_num': iat,
                                       'cart_dim': 1,
                                       'lk': (1, 2)
                                       }
                            )
    Y = load_QM.find_in_Qlk(rlk_data, {'at_num': iat,
                                       'cart_dim': 2,
                                       'lk': (1, 2)
                                       }
                            )
    Z = load_QM.find_in_Qlk(rlk_data, {'at_num': iat,
                                       'cart_dim': 3,
                                       'lk': (1, 2)
                                       }
                            )

    mag = np.sqrt(X**2 + Y**2 + Z**2)
    allAtRlkData.append(mag)

allAtRlkData = np.array(allAtRlkData)

f, a = plt.subplots()
iat = 9
for timesteps, data in zip(allRepTimesteps, allAtRepData[iat]):

    # Plot all reps (normalised)
    a.plot(timesteps,
           data / np.max(data),
           lw=0.5,
           color='k',
           alpha=0.1)

# Plot Rlk Data
data = allAtRlkData[iat]
a.plot(timesteps,
       data / np.max(data),
       lw=1,
       color='y',
       alpha=1)



a.set_xlabel("Timestep [fs]")
a.set_ylabel(r"|Qlk|$^2$ (atom %i) [fs]" % iat)
