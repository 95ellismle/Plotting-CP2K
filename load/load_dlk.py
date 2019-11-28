#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 15 13:52:01 2018

@author: mellis
"""

from load import load_xyz as XYZ
from load import load_utils as Utils

import pandas as pd
import os


# Reads 1 NACV file
def load_dlk(filepath,
             min_step=0,
             max_step='all',
             stride=1,
             ignore_steps=[]):

    if not os.path.isfile(filepath):
       print("Can't find filepath %s" % filepath)
       raise SystemExit("BREAK")

    # Read the data
    df = pd.read_csv(filepath)

    return df


# Reads all the dlk files from a given folder
def load_all_dlk_in_folder(folder,
                           min_step=0,
                           max_step='all',
                           stride=1,
                           ignore_steps=[],
                           reps='all'):

    return Utils.load_all_in_folder(folder=folder,
                                    func=load_dlk,
                                    args=[min_step,
                                          max_step,
                                          stride,
                                          ignore_steps],
                                    filename_must_contain=['NACV', 'csv'],
                                    filename_must_not_contain=[],
                                    reps=reps)

