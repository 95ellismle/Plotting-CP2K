#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 14:20:52 2018

@author: mellis
"""


from load import load_xyz as XYZ
from load import load_utils as Utils

#import numpy as np


# Reads 1 frc file
def load_pos(filepath, min_step=0, max_step='all', stride=1, ignore_steps=[]):
    """
    Load forces
    """
    data, cols, timesteps = XYZ.read_xyz_file(filename=filepath,
                                              num_data_cols=False,
                                              min_step=min_step,
                                              max_step=max_step,
                                              stride=stride,
                                              ignore_steps=ignore_steps)
    return [data, cols, timesteps]

# Reads all the Qlk files from a given folder
def load_all_vel_in_folder(folder, min_step=0, max_step='all', stride=1, ignore_steps=[], reps='all'):
    """
    Loads all forces in a folder
    """
    return Utils.load_all_in_folder(folder=folder,
                                    func=load_pos,
                                    args=[min_step, max_step, stride, ignore_steps],
                                    filename_must_contain=['vel','xyz'],
                                    filename_must_not_contain=['init'],
                                    reps=reps)

# Reads all the Qlk files from a given folder
def load_all_pos_in_folder(folder, min_step=0, max_step='all', stride=1, ignore_steps=[], reps='all'):
    """
    Loads all forces in a folder
    """
    return Utils.load_all_in_folder(folder=folder,
                                    func=load_pos,
                                    args=[min_step, max_step, stride, ignore_steps],
                                    filename_must_contain=['pos','xyz'],
                                    filename_must_not_contain=['init'],
                                    reps=reps)
