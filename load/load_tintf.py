#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 19 16:02:39 2018

@author: mellis
"""

from load import load_xyz as XYZ
from load import load_utils as Utils

#import numpy as np

# Reads 1 QM file
def load_tintf(filepath, min_step=0, max_step='all', stride=1, ignore_steps=[]):
    data, cols, timesteps = XYZ.read_xyz_file(filename=filepath, 
                                              num_data_cols =3,
                                              min_step=min_step, 
                                              max_step=max_step, 
                                              stride=stride, 
                                              ignore_steps=ignore_steps)
    cols[:,:,1] = cols[:,:,1].astype(int)
    return (data, cols), timesteps

# Reads all the Qlk files from a given folder
def load_all_tintf_in_folder(folder, min_step=0, max_step='all', stride=1, ignore_steps=[], reps='all'):
    return Utils.load_all_in_folder(folder=folder, 
                                    func=load_tintf, 
                                    args=[min_step, max_step, stride, ignore_steps], 
                                    filename_must_contain=['t_int_frc','xyz'], 
                                    filename_must_not_contain=[], 
                                    reps=reps)

