#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 19 16:02:39 2018

@author: mellis
"""

from load import load_xyz as XYZ
from load import load_utils as Utils

import numpy as np

# Reads 1 QM file
def load_tintf(filepath, min_step=0, max_step='all', stride=1, ignore_steps=[]):
    data, cols, timesteps = XYZ.read_xyz_file(filename=filepath,
                                              num_data_cols =3,
                                              min_step=min_step,
                                              max_step=max_step,
                                              stride=stride,
                                              ignore_steps=ignore_steps)
    cols[:,:,1] = cols[:,:,1].astype(int)
    return data, cols, timesteps

# Reads all the Qlk files from a given folder
def load_all_tintf_in_folder(folder, min_step=0, max_step='all', stride=1, ignore_steps=[], reps='all'):
    return Utils.load_all_in_folder(folder=folder,
                                    func=load_tintf,
                                    args=[min_step, max_step, stride, ignore_steps],
                                    filename_must_contain=['t_int_frc','xyz'],
                                    filename_must_not_contain=[],
                                    reps=reps)


def find_in_histF(data, cols, params):
    """
    Will find the data corresponding to the params in the dictionary.

    N.B. Integer index starts at 1

    Inputs:
        * data = just the data from load_tintf
        * cols = just the cols from load_tintf
        * params = Dictionary of paramters:
            - 'step_num': int,
            - 'at_num': int,
            - 'state': int
    """
    params['state'] += 1
    # Get vital metadata
    num_atoms = sum(cols[0, :, 1] == '1')
    num_states = len(set(cols[0, :, 1]))

    # Get the data from atom i
    if 'at_num' in params:
        at_nums = [params['at_num'] + (num_atoms*i) for i in range(num_states)]
        data = data[:, at_nums]
        cols = cols[:, at_nums, :]

    # Get data from state i
    if 'state' in params:
        mask = cols[:, :, 1] == str(params['state'])
        data = np.array([data[i][mask[i]] for i in range(len(mask))])
        cols = np.array([cols[i][mask[i]] for i in range(len(mask))])
    
    # Get data from step i
    if 'step_num' in params:
        data = data[params['step_num']]
        cols = cols[params['step_num']]

    return data, cols
