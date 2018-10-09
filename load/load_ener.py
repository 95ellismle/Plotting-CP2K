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


def load_ener_time_data(filepath):
    """
    Will load the run-ener file which contains the kinetic, potential consvered
    quantity and time taken etc... 
    
    Inputs:
        * filepath => the filepath of the file that should be loaded [str]
    
    Outputs:
        * dataframe containing parsed ener data. [pd.DataFrame]
    """
    cols = ['step',
            'curr_time',     
            'E_kin',     
            'temp',     
            'E_pot',     
            'E_tot', 
            'time_per_step']
    df = pd.read_csv(filepath, delim_whitespace=True, names=cols, skiprows=[0])
    
    return df

def load_all_ener_time_in_folder(folder, reps):
    """
    Will load all ener files in a folder. This means the kinetic, potential and 
    conserved energies and the time taken.
    
    Inputs:
        * folder => the folder to look in [str]
        * reps   => Which reps to load (can be 'all') [str] or [list of int]
    
    Outpus:
        * dictionary of each file with parse energy data. Filename is key and 
          data is value.
    """
    return Utils.load_all_in_folder(folder, 
                                    load_ener_time_data, 
                                    filename_must_contain=['csv', 'ener'], 
                                    filename_must_not_contain['ad'], 
                                    reps=reps)
