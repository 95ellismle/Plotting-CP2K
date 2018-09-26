#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  7 15:39:12 2018

@author: mellis
"""

import os
import re
import collections


# Will find the replica number in the filename
def find_rep_num_in_name(filename):
#    labs = ['pos','ham','ener','coeff','vel']
    regex_matches = {'n-pos':"-\d*-\d*\.xyz", 
                     'ham':"-\d*-\d*\.xyz",
                     'n-ener_':"_\d*-\d*\.dat",
                     'coeff-':"f-\d*-\d*\.xyz",
                     'n-vel':"-\d*-\d*\.xyz",
                     'n-frc':"_\d*-\d*\.xyz",
                     't_frc':"_\d*\.xyz",
                     'd_frc':"_\d*\.xyz",
                     "coeff_a":"_\d*-\d*\.xyz",
                     "QM":"-\d*\.xyz",
                     "d_ener":"-\d*-\d*\.csv",}
    regex_matches2 = {'n-pos':"\d*-\d*\.", 
                     'ham':"\d*-\d*\.",
                     'n-ener_':"\d*-\d*\.",
                     'coeff-':"\d*-\d*\.",
                     'n-vel':"\d*-\d*\.",
                     'n-frc':"\d*-\d*\.",
                     't_frc':"\d*\.",
                     'd_frc':"\d*\.",
                     "coeff_a":"\d*-\d*\.",
                     "QM":"\d*\.",
                     "d_ener":"\d*-\d*\.",}
    final_delim = {'n-pos':"-", 
                   'ham':"-",
                   'n-ener_':"-",
                   'coeff-':"-",
                   'n-vel':"-",
                   'n-frc':"-",
                   't_frc':".",
                   'd_frc':".",
                   "coeff_a":"-",
                   "QM":".",
                   "d_ener":"-",}
    
    for lab in regex_matches:
        if lab in filename:
            reduced_fname = re.findall(regex_matches[lab], filename)
            if len(reduced_fname) == 1:
                rep = re.findall(regex_matches2[lab], reduced_fname[0])
                if len(rep) == 1:
                    for ichar, char in enumerate(rep[0]):
                        if char == final_delim[lab]:
                            try:    
                                rep = int(rep[0][:ichar])
                                return rep
                            except TypeError:
                                print("Sorry something went wrong extracting the replica number from a file")
                                print("Cannot convert %s to an integer"%str(rep[0][:ichar]))
                                print("\nFilename = %s"%filename)
                                print("Please make sure the replica number is followed by a hyphen i.e. 'run_pos-(1-)1.xyz'")
                                raise SystemExit("Type Conversion Error")
                            break
                    else:
                        print("Cannot find a hyphen in the filename.")
                        print("Please make sure the replica number is followed by a hyphen i.e. 'run_pos-(1-)1.xyz'")
                        print("\nFilename = %s"%filename)
                        raise SystemExit("Missing Final Delimeter")
                    break
                else:
                    raise SystemExit("""Sorry I couldn't find the replica number with the regex.
                                 
Filename = '%s'            regex = '\d*-\d*\.'

Filename after regex = '%s' """%(reduced_fname, str(rep)))
            else:
                raise SystemExit("""Sorry I the pre-programmed regex doesn't work for this file.
                                 
Filename = '%s'            regex = '%s'

Filename after regex = '%s' """%(filename, regex_matches[lab], str(reduced_fname)))
    else:
        raise SystemExit("Sorry I couldn't find the file type (pos, vel, frc, coeff etc...), something went wrong!\n\nFilename = %s"%(filename))     


# Will find how many/which replicas are available and return them based on which ones are wanted
def sort_reps(folder, reps):
    all_files = os.listdir(folder)
    labs = ['n-frc', 'd_frc', 'coeff_a', 'n-vel', 'n-pos', 't_frc', 'ham', 'coeff-', 'd_ener', 'QM', 'n-ener_',]
    files_sorted = {lab:
                   [f for f in all_files if lab in f] for lab in labs }
    for lab in labs:
        avail_reps = [find_rep_num_in_name(f) for f in files_sorted[lab]]
        break
            
    if type(reps) == str:
        reps = avail_reps
    elif type(reps) == list:
        reps = [i for i in reps if i in avail_reps]
        
    return reps


# Given a list of files and a list of replica numbers find which files match the list of reps
def files_with_correct_reps(files, reps):
    rep_nums = {f:find_rep_num_in_name(f) for f in files}
    files = [i[1] for i in sorted(zip(rep_nums.values(), rep_nums.keys()))]
    if type(reps) == int:
        reps = [reps]
    elif 'all' in reps:
        return files
    files = [f for f in files if rep_nums[f] in reps]
    return files    

# Will apply a function to load all relevant files in a folder
def load_all_in_folder(folder, func, args=[], filename_must_not_contain=[], filename_must_contain=[], reps='all'):
    if not os.path.isdir(folder):
            raise SystemError("Sorry the folder given can't be found...")
    else:
        if folder[-1] != '/':
            folder = folder + '/'

    filename_must_not_contain.append(".sw")
    all_files = [folder+i for i in os.listdir(folder) if all(j in i for j in filename_must_contain) and all(k not in i for k in filename_must_not_contain)]
    all_files = files_with_correct_reps(all_files, reps) # Only read files with the correct rep num
    args = [[f]+list(args) for f in all_files]
    all_data = collections.OrderedDict()
    for arg in args:
        all_data[arg[0][arg[0].rfind('/')+1:]] = func(*arg)
    return all_data






