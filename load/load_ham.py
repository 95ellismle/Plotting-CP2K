#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  7 15:24:40 2018

@author: mellis
"""

from load import load_xyz as XYZ
from load import load_utils as Utils

import os

# Will try and find the dimension of the hamiltonian
def find_num_basis(filepath):
    with open(filepath, 'r') as f:
        for line in f:
            all_line = [i for i in line.split(' ') if i.strip() and '\n' != i and '\n' != i.strip()]
            num_line = [i for i in all_line if XYZ.is_num(i)]
            if float(len(num_line))/len(all_line) > 0.8:
                floats_in_num_line = [i for i in num_line if '.' in i]
                num_basis = len(floats_in_num_line)
                break
        else:
            raise SystemExit("Sorry I can't find the dimension of the hamiltonian!\n\nPlease look at the load_ham.py file in load/")
    return num_basis

# Will load a single hamiltonian file
def load_ham(filepath,
             num_basis,
             metadata,
             min_step=0,
             max_step='all',
             stride=1,
             ignore_steps=[]):
    data, cols, timesteps = XYZ.read_xyz_file(filepath, 
                                              num_basis,
                                              min_step=min_step, 
                                              max_step=max_step, 
                                              stride=stride, 
                                              ignore_steps=ignore_steps,
                                              metadata=metadata)
    return [data*27000, cols, timesteps]


# Reads all the Qlk files from a given folder
def load_all_ham_in_folder(folder,
                           min_step=0,
                           max_step='all',
                           stride=1,
                           ignore_steps=[],
                           reps='all'):
    allHamFiles = [i for i in os.listdir(folder) if 'run-ham' in i]
    filepath = folder + allHamFiles[0]

    num_basis = find_num_basis(filepath)
    metadata = XYZ.get_xyz_step_metadata2(filepath)
    return Utils.load_all_in_folder(folder=folder,
                                    func=load_ham,
                                    args=[num_basis,
                                          metadata,
                                          min_step,
                                          max_step,
                                          stride,
                                          ignore_steps],
                                    filename_must_contain=['ham','xyz'],
                                    filename_must_not_contain=[],
                                    reps=reps)


## Will find and load all hamiltonian files
#def load_all_ham_in_folder(folder, min_step=0, max_step='all', stride=1, ignore_steps=[], filename_must_not_contain='', filename_must_contain='', reps='all'):
#    if not os.path.isdir(folder):
#        raise SystemError("Sorry the folder given can't be found...")
#    else:
#        if folder[-1] != '/':
#            folder = folder + '/'
#    if filename_must_not_contain:
#        ham_files = [folder+i for i in os.listdir(folder) if '.xyz' in i and 'ham' in i and filename_must_not_contain not in i and filename_must_contain in i]
#    else:
#        ham_files = [folder+i for i in os.listdir(folder) if '.xyz' in i and 'ham' in i and filename_must_contain in i]
#    
#    ham_files = Utils.files_with_correct_reps(ham_files, reps) # Only read files with the correct rep num
#
#    all_ham_data = {f[f.rfind('/')+1:]: load_ham(f, 
#                                    min_step=min_step, 
#                                    max_step=max_step, 
#                                    stride=stride, 
#                                    ignore_steps=ignore_steps) for f in ham_files}
#    return all_ham_data
#
#    
