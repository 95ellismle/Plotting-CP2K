#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  7 11:32:36 2018

@author: mellis
"""
import numpy as np
import os

from load import load_xyz as XYZ
from load import load_utils as Utils

# Will load a single coefficient file
def load_coeff(filepath, min_step=0, max_step='all', stride=1, ignore_steps=[]):
    data, cols, timesteps = XYZ.read_xyz_file(filename=filepath, 
                                              num_data_cols=2, 
                                              min_step=min_step, 
                                              max_step=max_step, 
                                              stride=stride, 
                                              ignore_steps=ignore_steps)
    pops = np.linalg.norm(data, axis=2)**2
    data = np.array([np.array([complex(*j) for j in i]) for i in data])
    return np.array(data), np.array(cols), np.array(timesteps), np.array(pops)

# Reads all the Qlk files from a given folder
def load_all_coeff_in_folder(folder, min_step=0, max_step='all', stride=1, ignore_steps=[], filename_must_not_contain=[], filename_must_contain=[], reps='all'):
    return Utils.load_all_in_folder(folder=folder, 
                                    func=load_coeff, 
                                    args=[min_step, max_step, stride, ignore_steps], 
                                    filename_must_not_contain=filename_must_not_contain, 
                                    filename_must_contain=filename_must_contain, 
                                    reps=reps)

