#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 12:44:32 2019

@author: mellis
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

rootfolder = "/homes/mellis/Documents/Graphs/Tully_Models"
modelNum = 1
initialMom = 25

things =  (rootfolder, modelNum, initialMom)
popsFilePath = "%s/Model%i/CTMQC_%iK/Pops.csv" % things
decoFilePath = "%s/Model%i/CTMQC_%iK/Deco.csv" % things
names = ['exact', 'SH', 'Eh', 'MQC', 'CTMQC']

for i, fPath in enumerate((popsFilePath, decoFilePath)):
    with open(fPath, 'r') as f:
        lines = f.read().split('\n')[2:]
        headers = [name + j for name in names for j in ['_x', '_y']]
        txt = [','.join(headers)] + lines
    
    with open("tmp.csv", 'w') as f:
        f.write('\n'.join(txt))
        
    data = pd.read_csv("tmp.csv")

    # Sort by x data
    for name in names:
        xname, yname = name+'_x', name+'_y'
        data[yname][data[yname] > 1] = 1
        data[yname][data[yname] < 0] = 0
        
    colors = ['r', 'b', (0, 1, 1), (0.7, 0.7, 0), 'g']
    
    plt.figure()
    for col, name in zip(colors, names):
        xname, yname = name+'_x', name+'_y'
        xdata, ydata = data[xname], data[yname]
    
        sortedData = sorted(zip(xdata, ydata))
        xdata = np.array([i[0] for i in sortedData])
        ydata = np.array([i[1] for i in sortedData])
        
        plt.plot(xdata, ydata, color=col, ls='-', lw=3)
#        plt.plot(xdata, 1-ydata, color=col, ls='-')
    
    plt.show()
    
    
    os.remove("tmp.csv")