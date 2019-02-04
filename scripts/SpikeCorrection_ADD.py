#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 13:25:52 2019

@author: mellis
"""

import matplotlib.pyplot as plt
import numpy as np

width = 0.001
height = 1e-5

x1 = np.arange(0, 0.01, 0.00001)
y1 = height * np.exp(-(x1 / width)**2)

x2 = np.arange(-0.01, 0, 0.00001)
y2 = -height * np.exp(-(x2 / width)**2)

plt.plot(x1, y1)
plt.plot(x2, y2)

plt.xlabel(r"$\sum_{J} Y^{(I)}_{lk, \nu}(t)$")
plt.ylabel("Add-on Correction")
