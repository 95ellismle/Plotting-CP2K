#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 17 11:38:11 2018

@author: mellis
"""

import pandas as pd

from load import load_utils as Utils


# Will load the total energy file
def load_ener_dat(filepath):
    headers = ['Step','Time','Kin','Temp','Pot','E_cons','CPU']
    data = pd.read_csv(filepath, names=headers, delim_whitespace=True, skiprows=[0])
    return data

# Load all adiabatic energy files in a folder
def load_all_ener_dat(folder, reps):
    return Utils.load_all_in_folder(folder, load_ener_dat, filename_must_contain=['dat', 'ener'], reps=reps)


# Will load the adiabtic energy csv
def load_ener_ad(filepath):  return pd.read_csv(filepath)

# Load all adiabatic energy files in a folder
def load_all_ener_ad(folder, reps):
    return Utils.load_all_in_folder(folder, load_ener_ad, filename_must_contain=['csv', 'ad_ener'], reps=reps)

