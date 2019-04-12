#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 17 11:38:11 2018

@author: mellis
"""

import pandas as pd

from load import load_utils as Utils


# Will load the total energy file
def load_ener_dat(filepath, max_time):
    headerReplacers = {'StepNr': 'Step','Time_[fs]':'Time',
                       'Kinetic_[a.u.]': 'Kin',
                       'Temperature_[K]': 'Temp',
                       'Potential_[a.u.]': 'Pot',
                       'ConsQty_[a.u.]': 'E_cons',
                       'CPU_[s]': 'CPU'}

    headers = []
    with open(filepath, 'r') as f:
      for line in f:
         if line:
            for elem in line.split():
               if elem in headerReplacers:
                   headers.append(headerReplacers[elem])
               else:
                   headers.append(elem)
            break
    data = pd.read_csv(filepath, names=headers, delim_whitespace=True, skiprows=[0])
    if max_time == 'all':
        return data
    else:
        return data[data['Time'] < max_time]
    
# Load all adiabatic energy files in a folder
def load_all_ener_dat(folder, reps, max_time):
    return Utils.load_all_in_folder(folder, 
                                    load_ener_dat, 
                                    args=[max_time], 
                                    filename_must_contain=['dat', 'ener'], 
                                    reps=reps)


# Will load the adiabtic energy csv
def load_ener_ad(filepath, max_time, min_time):
    data = pd.read_csv(filepath)
    if min_time > 0:
        data = data[data['Time'] > min_time]
    if max_time == 'all':
        return data
    else:
        return data[data['Time'] < max_time]
    
# Load all adiabatic energy files in a folder
def load_all_ener_ad(folder, reps, max_time='all', min_time=0):
    return Utils.load_all_in_folder(folder, load_ener_ad, args=[max_time, min_time], filename_must_contain=['csv', 'ad_ener'], reps=reps)


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
    return Utils.load_all_in_folder(folder = folder, 
                                    func = load_ener_time_data, 
                                    filename_must_contain = ['csv', 'ener'], 
                                    filename_must_not_contain = ['ad'], 
                                    reps = reps)
