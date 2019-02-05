#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 11:54:01 2019

@author: oem
"""
from load import load_xyz as XYZ
from load import load_utils as Utils

import numpy as np


def load_K(filename, min_step=0, max_step='all',
           stride=1, ignore_steps=[]):
    """
    Will load a sigma.list file. This contains info on the nuclear width.

    Inputs:
        * filename     = filename to open and parse
        * min_step     = first step to read
        * max_step     = last step to read
        * stride       = stride to take when reading
        * ignore_steps = which steps to ignore
    """
    with open(filename, 'r') as f:
        txt = f.read()
    ltxt = txt.split('\n')

    # Find the timesteps
    time_delim = XYZ.find_time_delimeter(ltxt[:1], filename)
    timesteps = [XYZ.string_between(ltxt[i], 'e =', time_delim[0])
                 for i in range(0, len(ltxt), 2)]
    timesteps = np.array([float(i) for i in timesteps if i.strip()])

    widths = np.array([ltxt[i].split(' ') for i in range(1, len(ltxt), 2)])
    widths = np.array([[float(j.strip(',')) for j in i if j] for i in widths])

    return widths, timesteps


def load_all_K_in_folder(folder, min_step=0, max_step='all',
                         stride=1, ignore_steps=[],
                         filename_must_not_contain=[],
                         filename_must_contain=[], reps='all'):
    """
    Will load all the K files in a specified folder. This will output an
    Ordered Dictionary with filenames as keys and data as values.

    Inputs:
        * folder = the folder to look for sigma files in
        * min_step     = first step to read
        * max_step     = last step to read
        * stride       = stride to take when reading
        * ignore_steps = which steps to ignore
        * filename_must_contain = any strings the filename must contain (list)
        * filename_must_not_contain = any strings the filename must not contain
                                      (list)
        * reps         = which replicas to get data for
    """
    filename_must_contain.append("K")
    filename_must_contain.append("list")
    return Utils.load_all_in_folder(folder=folder,
                                    func=load_K,
                                    args=[min_step,
                                          max_step,
                                          stride,
                                          ignore_steps],
                                    filename_must_not_contain=filename_must_not_contain,
                                    filename_must_contain=filename_must_contain, 
                                    reps=reps)
