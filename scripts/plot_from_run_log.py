#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 10:45:48 2018

@author: mellis
"""

import numpy as np
import matplotlib.pyplot as plt


folder = "/scratch/mellis/flavoured-cptk/200Rep_2mol/"
filepath = folder + 'run.log'
find_str = "Check2"


with open(filepath, 'r') as f:
    found_lines = [line.split('=')[1].strip('\n').strip(' ') for line in f if find_str in line]    
    found_lines = np.array(found_lines).astype(float)
    print(found_lines)

with open(filepath.replace(".log", ".inp"), 'r') as f:
    for line in f:
        if 'TANH_WIDTH' in line:
            tanh_width = line.strip(" ").strip("TANH_WIDTH").strip()
            tanh_width = float(tanh_width)
            

plt.plot(np.arange(len(found_lines))*0.1, found_lines)
plt.xlabel("Time [fs]")
plt.ylabel("Equation S26")
plt.ylim([-0.001,0.0004])

plt.title("Check of Equation S26 (it should stay at 0)")
