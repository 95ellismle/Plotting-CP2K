#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 11:00:42 2018

@author: mellis
"""

import re
import os

from IO import Folders as fold


def get_all_run_inp_variables(filepath):
    """
    Will parse all the variables in a run.inp file and return them in a 
    dictionary. This dictionary will have the setting name as a key and the
    setting value as the value.
    
    Inputs:
        * filepath => filepath of the run.inp file [str]
    """
    with open(filepath, 'r') as f:
        ltxt = [i for i in f.read().split('\n') if i]
        variables = {}
        for line in ltxt:
            line = line.strip()
            
            #Don't want comment lines or include lines
            if any(line[0] == j for j in ('#','&','@','!')): continue            
            line = line.split('!')[0]
            line = line.split('#')[0]
            line = [i for i in line.split(' ') if i]
            
            # Remove unit info
            line = [i for i in line if '[' not in i and ']' not in i]
            
            if len(line) == 2:
                variables[line[0].upper()] = line[1]
            elif len(line) == 1:
                variables[line[0].upper()] = True
            else:
                print("""The line
                      
        ''' %s '''

seems to be causing some trouble.

It's length after processing and being split by ' ' is not 2 or 1."""%line)

    # Remove unwanted settings
    for i in ['MD','FORMAT','UNIT']:
        if variables.get(i): variables.pop(i)
        
    return variables



def get_run_inp_in_folder(folder):
    """
    Will parse the variables from the run.inp files in a folder. These are 
    files which match the 'run*.\.inp' syntax. Will return a list of 
    run.inp files with a list of filenames.
    
    Inputs:
        * folder  =>  the folder in which to for the run.inp files
    """
    folder = fold.make_fold_abs(folder)
    all_files = os.listdir(folder)
    run_inp_files = [re.findall('run*.\.inp',i) for i in all_files]
    run_inp_files = [folder+i[0] for i in run_inp_files if i]
    run_inp_vars = [get_all_run_inp_variables(i) for i in run_inp_files]
    
    return run_inp_vars, run_inp_files

    