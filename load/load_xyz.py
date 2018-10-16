#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  7 11:34:43 2018

@author: mellis
"""

import os 
import difflib as dfl
from collections import OrderedDict
import numpy as np

# Returns the substring between 2 substrings within a string
def string_between(Str, substr1, substr2):
    Str = Str[Str.find(substr1)+len(substr1):]
    Str = Str[:Str.find(substr2)]
    return Str

# Checks if a line of text is an atom line in a xyz file
def is_atom_line(line):
    line = [i for i in line.split(' ') if i]
    if len(line) < 3:
        return False
    fline = [is_num(i) for i in line]
    if not sum(fline[-3:]) == 3:
        return False
    else:
        return True

# Checks whether a string can be a number
def is_num(Str):
    try:
        float(Str)
        return True
    except:
        return False

# Will determine the number of atoms in an xyz file
def num_atoms_find(ltxt):
    start_atoms, finish_atoms = 0,0
    for i,line in enumerate(ltxt):
        if (is_atom_line(line)) == True:
            start_atoms = i
            break
    for i,line in enumerate(ltxt[start_atoms:],start=start_atoms):
        if (is_atom_line(line) == False):
            finish_atoms=i
            break
    return start_atoms, finish_atoms-start_atoms

# Finds the number of lines in one step of the xyz file data
def find_num_lines_in_xyz_file_step(ltxt, filename):
    first_line = ltxt[0]
    num_lines_in_step = 1
    for i in ltxt[2:]: # Loops over all the line of text
        num_lines_in_step += 1
        #If any lines are very similar to the first line then assume the step is repeating
        if dfl.SequenceMatcher(None, first_line, i).ratio() > 0.8:
            return num_lines_in_step
    raise SystemExit("Unable to determine number of steps in:\n\n%s"%filename)

# Finds the number of title lines and number of atoms with a step
def find_num_title_lines(step): # should be the text in a step split by line
    num_title_lines = 0
    for line in step:
        if is_atom_line(line):
            break
        num_title_lines += 1
    return num_title_lines

# Finds the delimeter for the time-step in the xyz_file title
def find_time_delimeter(step, filename):
    for linenum,txt in enumerate(step):
        txt = txt.lower()
        if 'time' in txt:
            break
    else:
        raise SystemExit ("Can't find the word 'time' in this data:\n\n%s\n\n\tFilename:%s"%(str(step), filename) )
    prev_char, count = False, 0
    txt = txt[txt.find("time"):]
    for char in txt.replace(" ",""):
        isnum = (char.isdigit() or char == '.')
        if isnum != prev_char:
            count += 1
        prev_char = isnum
        if count == 2:
            break
    if char.isdigit(): return '\n', linenum
    else: return char, linenum
    raise SystemExit("Cannot find the delimeter for the time-step info in the following xyz_file:\n\n%s\n\nstep = %s"%(filename,step))

# Will get necessary metadata from an xyz file such as time step_delim, lines_in_step etc...
# This will also create the step_data dictionary with the data of each step in
def get_xyz_step_metadata(ltxt, filename):
    most_stable = False
    if any('*' in i for i in ltxt[:300]):
        most_stable = True
    if not most_stable:
        num_title_lines, num_atoms = num_atoms_find(ltxt)
        lines_in_step = num_title_lines + num_atoms
        if len(ltxt) > lines_in_step+1: # take lines from the second step instead of first as it is more reliable
           step_data = {i: ltxt[i*lines_in_step:(i+1)*lines_in_step] for i in range(1,2)}
        else: #If there is no second step take lines from the first
           step_data = {1:ltxt[:lines_in_step]}
    else:
        lines_in_step = find_num_title_lines(ltxt)
        step_data = {i: ltxt[i*lines_in_step:(i+1)*lines_in_step] for i in range(1,2)}
        num_title_lines = find_num_title_lines(step_data[1])
    time_delim, time_ind = find_time_delimeter(step_data[1][:num_title_lines], filename)
    return time_delim, time_ind, lines_in_step, num_title_lines

# Reads an xyz_file
# Would like to create a mask here to avoid reading the atoms to ignore.
# This function is quite obscure and terse as this is a bottle neck in the code and has been optimised.
# It relies heavily on numpy arrays and list\dictionary comphrensions to give speed things up.
def read_xyz_file(filename, num_data_cols, min_step=0, max_step='all', stride=1, ignore_steps=[]):
    num_data_cols = -num_data_cols
    ltxt = open_read(filename).split('\n')
    time_delim, time_ind, lines_in_step, num_title_lines = get_xyz_step_metadata(ltxt, filename)
    abs_max_step = int(len(ltxt)/lines_in_step)
    if max_step == 'all' or max_step > abs_max_step:
        max_step = abs_max_step
    # The OrderedDict doesn't seem to have major overheads as dictionary access aren't the main bottleneck here.
    # It is also much easier to use!
    step_data = OrderedDict() # The OrderedDict keeps the order of the frames for saving etc...
    all_steps = [i for i in range(min_step, max_step, stride) if i not in ignore_steps]
    for i in all_steps:
        step_data[i] = ltxt[i*lines_in_step:(i+1)*lines_in_step]
        step_data[i] = (step_data[i][:num_title_lines],step_data[i][num_title_lines:])
    time_steps = np.array([string_between(step_data[i][0][time_ind], "time = ", time_delim) for i in step_data]).astype(float)
    for i in all_steps:
        step_data[i] = [j.split(' ') for j in step_data[i][1]]
        step_data[i] = np.array([[k for k in j if k] for j in step_data[i]])
    data = np.array([step_data[i][:,num_data_cols:] for i in step_data]).astype(float)
    spare_info = np.array([step_data[i][:,:num_data_cols] for i in step_data])
    return data, spare_info, time_steps


# Reads a file and closes it
def open_read(filename, throw_error=True):
    if os.path.isfile(filename):
        f = open(filename, 'r')
        txt = f.read()
        f.close()
        return txt
    else:
        if throw_error:
            raise SystemExit("The %s file doesn't exist!"%filename)
        return False