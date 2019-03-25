#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 15 13:52:01 2018

@author: mellis
"""

from load import load_xyz as XYZ
from load import load_utils as Utils

import numpy as np

params_to_ind_conv = {'at_num': 0,
                      'l': 2,
                      'k': 3,
                      'cart_dim': 1}


def load_QM_0(filepath,
              min_step=0,
              max_step='all',
              stride=1,
              ignore_steps=[]):
    """
    Will load a single traj QM_0 (the quantum momentum without the lk indices)
    """
    data, cols, timesteps = XYZ.read_xyz_file(filename=filepath,
                                              num_data_cols=3,
                                              min_step=min_step,
                                              max_step=max_step,
                                              stride=stride,
                                              ignore_steps=ignore_steps)
    cols = np.array(cols).astype(int)
    return (data, cols), timesteps


# Reads all the Qlk files from a given folder
def load_all_QM_0_in_folder(folder,
                            min_step=0,
                            max_step='all',
                            stride=1,
                            ignore_steps=[],
                            reps='all'):
    """
    Will load all the QM_0 files in 1 folder.
    """
    return Utils.load_all_in_folder(folder=folder,
                                    func=load_QM_0,
                                    args=[min_step,
                                          max_step,
                                          stride,
                                          ignore_steps],
                                    filename_must_contain=['QM_0', 'xyz'],
                                    filename_must_not_contain=[],
                                    reps=reps)


# Reads 1 QM file
def load_Qlk(filepath,
             min_step=0,
             max_step='all',
             stride=1,
             ignore_steps=[]):

    data, cols, timesteps = XYZ.read_xyz_file(filename=filepath,
                                              num_data_cols=1,
                                              min_step=min_step,
                                              max_step=max_step,
                                              stride=stride,
                                              ignore_steps=ignore_steps)
    cols = np.array(cols).astype(int)
    return (data, cols), timesteps


# Reads all the Qlk files from a given folder
def load_all_Qlk_in_folder(folder,
                           min_step=0,
                           max_step='all',
                           stride=1,
                           ignore_steps=[],
                           reps='all'):

    return Utils.load_all_in_folder(folder=folder,
                                    func=load_Qlk,
                                    args=[min_step,
                                          max_step,
                                          stride,
                                          ignore_steps],
                                    filename_must_contain=['QM', 'xyz'],
                                    filename_must_not_contain=[],
                                    reps=reps)


# Reads all the Qlk files from a given folder
def load_all_Rlk_in_folder(folder,
                           min_step=0,
                           max_step='all',
                           stride=1,
                           ignore_steps=[],
                           reps='all'):

    rlk_data = Utils.load_all_in_folder(folder=folder,
                                        func=load_Qlk,
                                        args=[min_step,
                                              max_step,
                                              stride,
                                              ignore_steps],
                                        filename_must_contain=['rlk', 'xyz'],
                                        filename_must_not_contain=[],
                                        reps=reps)
    keys = list(rlk_data.keys())
    if len(keys) > 1:
        files = '\n\t*'.join(keys)
        raise SystemExit("Found more than 1 Rlk file!\n\n\t*%s" % files)

    rlk_data = rlk_data[keys[0]]
    return rlk_data


# Will get a value given the index in the Qlk data structure
def _get_in_Qlk_int_(data, cols, ind, val):
    if type(val) == int or type(val) == np.int64:
        if len(np.shape(data)) == 3:
            mask = cols[:, :, ind] == val
        elif len(np.shape(data)) == 2:
            mask = cols[:, ind] == val
        if len(mask) != len(data) or len(mask) != len(cols):
            minLen = np.min([len(mask), len(data), len(cols)])
            data = data[:minLen]
            cols = cols[:minLen]
            mask = mask[:minLen]

        data = data[mask]
        cols = cols[mask]
    else:
        msg = """
The `find_in_Qlk' function requires all values in the params dict to be a
 single integer.

You passed %s which is a %s
""" % (val, type(val))
        raise TypeError(msg)
    return data, cols


# Will return the list of all possible values of the particular index of Qlk
def _get_all_Qlk_col_vals_(cols, ind):
    if type(ind) == str:
        ind = params_to_ind_conv[ind]

    return list(set(cols[:, :, ind][0]))


# Will help in picking out certain parts of the Qlk data
def find_in_Qlk(Qlk_data, params):
    """
    Will find data corresponding to the params entered.

    Inputs:
        Qlk_data  =>  the qlk to search in (must just be a tuple of data and
                                            cols)
        params    =>  A dictionary with the parameters to search for.
                      The optional key labels are:
                          * at_num    [int]
                          * cart_dim  [int]
                          * lk        [list or tuple <int> dimension 2]
                          * step_num  [int or list of ints]
    """
    data, cols = Qlk_data
    if params.get("step_num") is not None:
        data = data[params['step_num']]
        cols = cols[params['step_num']]
    if params.get('at_num') is not None:
        data, cols = _get_in_Qlk_int_(data, cols, 0, params['at_num'])
    if params.get('cart_dim') is not None:
        data, cols = _get_in_Qlk_int_(data, cols, 1, params['cart_dim'])
    if params.get('lk') is not None:
        max_val, min_val = np.max(params['lk']), np.min(params['lk'])
        data, cols = _get_in_Qlk_int_(data, cols, 2, max_val)  # get l
        data, cols = _get_in_Qlk_int_(data, cols, 3, min_val)  # then get k
    return data
