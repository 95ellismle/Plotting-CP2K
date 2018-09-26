#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 15 13:52:01 2018

@author: mellis
"""

from load import load_xyz as XYZ
from load import load_utils as Utils

import numpy as np

params_to_ind_conv = {'at_num':0,
                      'l':2,
                      'k':3,
                      'cart_dim':1}

# Reads 1 QM file
def load_Qlk(filepath, min_step=0, max_step='all', stride=1, ignore_steps=[]):
    data, cols, timesteps = XYZ.read_xyz_file(filename=filepath, 
                                              num_data_cols =1,
                                              min_step=min_step, 
                                              max_step=max_step, 
                                              stride=stride, 
                                              ignore_steps=ignore_steps)
    cols = np.array(cols).astype(int)
    return (data, cols), timesteps

# Reads all the Qlk files from a given folder
def load_all_Qlk_in_folder(folder, min_step=0, max_step='all', stride=1, ignore_steps=[], reps='all'):
    return Utils.load_all_in_folder(folder=folder, 
                                    func=load_Qlk, 
                                    args=[min_step, max_step, stride, ignore_steps], 
                                    filename_must_contain=['QM','xyz'], 
                                    filename_must_not_contain=[], 
                                    reps=reps)

# Will get a value given the index in the Qlk data structure
def _get_in_Qlk_int_(data, cols, ind, val):
    if type(val) == int or type(val) == np.int64:
        if len(np.shape(data)) == 3:
            mask = cols[:,:,ind] == val
        elif len(np.shape(data)) == 2:
            mask = cols[:,ind] == val
        data = data[mask]
        cols = cols[mask]
    else:
        msg = """
The `find_in_Qlk' function requires all values in the params dict to be a single integer.

You passed %s which is a %s
"""%(val , type(val))
        raise TypeError(msg)
    return data, cols

# Will return the list of all possible values of the particular index of Qlk
def _get_all_Qlk_col_vals_(cols, ind):
    if type(ind) == str:
        ind = params_to_ind_conv[ind]
    return list(set(cols[:,:,ind][0]))



# Will help in picking out certain parts of the Qlk data
def find_in_Qlk(Qlk_data, params):
    data, cols = Qlk_data
    if params.get('at_num'):
        data, cols = _get_in_Qlk_int_(data, cols, 0, params['at_num'])
    if params.get('cart_dim'):
        data, cols = _get_in_Qlk_int_(data, cols, 1, params['cart_dim'])
    if params.get('lk'):
        max_val, min_val = np.max(params['lk']), np.min(params['lk'])
        data, cols = _get_in_Qlk_int_(data, cols, 2, max_val) #get l
        data, cols = _get_in_Qlk_int_(data, cols, 3, min_val) #then get k
    return data
