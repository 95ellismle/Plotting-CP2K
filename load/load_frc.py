#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 19 16:27:47 2018

@author: mellis
"""


from load import load_xyz as XYZ
from load import load_utils as Utils

#import numpy as np


# Reads all the Qlk files from a given folder
def load_all_qm_frc_in_folder(folder, min_step=0,
                              max_step='all', stride=1,
                              ignore_steps=[], reps='all'):
    """
    Loads all diabatic forces in a folder
    """
    return Utils.load_all_in_folder(folder=folder, 
                                    func=load_frc, 
                                    args=[min_step, max_step, stride, ignore_steps], 
                                    filename_must_contain=['run-QM_frc','xyz'], 
                                    filename_must_not_contain=['ad'], 
                                    reps=reps)



# Reads 1 frc file
def load_frc(filepath, min_step=0, max_step='all', stride=1, ignore_steps=[]):
    """
    Load diabatic forces
    """
    data, cols, timesteps = XYZ.read_xyz_file(filename=filepath, 
                                              num_data_cols =3,
                                              min_step=min_step, 
                                              max_step=max_step, 
                                              stride=stride, 
                                              ignore_steps=ignore_steps)
    return data, cols, timesteps


# Reads all the Qlk files from a given folder
def load_all_frc_in_folder(folder, min_step=0, max_step='all', stride=1, ignore_steps=[], reps='all'):
    """
    Loads all diabatic forces in a folder
    """
    return Utils.load_all_in_folder(folder=folder, 
                                    func=load_frc, 
                                    args=[min_step, max_step, stride, ignore_steps], 
                                    filename_must_contain=['run-frc_','xyz'], 
                                    filename_must_not_contain=['ad', 'QM'], 
                                    reps=reps)


# Reads 1 frc file
def load_ad_frc(filepath, min_step=0, max_step='all', stride=1, ignore_steps=[]):
    """
    Load adiabatic forces
    """
    data, cols, timesteps = XYZ.read_xyz_file(filename=filepath, 
                                              num_data_cols =3,
                                              min_step=min_step, 
                                              max_step=max_step, 
                                              stride=stride, 
                                              ignore_steps=ignore_steps)
    cols[:,:,1] = cols[:,:,1].astype(int)
    return data, cols, timesteps


# Reads all the Qlk files from a given folder
def load_all_ad_frc_in_folder(folder, min_step=0, max_step='all', stride=1, ignore_steps=[], reps='all'):
    """
    Loads all adiabatic forces in a folder
    """
    return Utils.load_all_in_folder(folder=folder, 
                                    func=load_ad_frc, 
                                    args=[min_step, max_step, stride, ignore_steps], 
                                    filename_must_contain=['frc','xyz', 'ad'], 
                                    filename_must_not_contain=[], 
                                    reps=reps)
