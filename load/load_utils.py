#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  7 15:39:12 2018

@author: mellis
"""

import os
import re
import collections
import multiprocessing as mp
import numpy as np


# Will find the replica number in the filename
def find_rep_num_in_name(filename):
    """
    Will find the replica number in the output filename. This should be only
    the filename and not the folder.
    """
#    labs = ['pos','ham','ener','coeff','vel']
    regex_matches = {'n-pos': r"-\d*-\d*\.xyz",
                     'ham': r"-\d*-\d*\.xyz",
                     'n-ener_': r"_\d*-\d*\.dat",
                     'coeff-': r"f-\d*-\d*\.xyz",
                     #'coeff_': r"f_\d*-\d*\.xyz",
                     'n-vel': r"-\d*-\d*\.xyz",
                     'n-frc': r"_\d*-\d*\.xyz",
                     't_frc': r"_\d*\.xyz",
                     'd_frc': r"_\d*\.xyz",
                     "coeff_a": r"_\d*-\d*\.xyz",
                     "QM-": r"-\d*\.xyz",
                     "QM_": r"-\d*\.xyz",
                     "sigma": r"-\d*\.list",
                     "d_ener": r"-\d*-\d*\.csv",
                     "n-K": r"-K-\d*.list",
                     }

    regex_matches2 = {'n-pos': r"\d*-\d*\.",
                      'ham': r"\d*-\d*\.",
                      'n-ener_': r"\d*-\d*\.",
                      'coeff-': r"\d*-\d*\.",
                      #'coeff_': r"\d*-\d*\.",
                      'n-vel': r"\d*-\d*\.",
                      'n-frc': r"\d*-\d*\.",
                      't_frc': r"\d*\.",
                      'd_frc': r"\d*\.",
                      "coeff_a": r"\d*-\d*\.",
                      "QM-": r"\d*\.",
                      "QM_": r"\d*\.",
                      "sigma": r"\d*\.",
                      "d_ener": r"\d*-\d*\.",
                      "n-K": r"\d*\.",
                      }
    final_delim = {'n-pos': r"-",
                   'ham': r"-",
                   'n-ener_': r"-",
                   'coeff-': r"-",
                   #'coeff_': r"-",
                   'n-vel': r"-",
                   'n-frc': r"-",
                   't_frc': r".",
                   'd_frc': r".",
                   "coeff_a": r"-",
                   "QM-": r".",
                   "QM_": r".",
                   "sigma": r".",
                   "d_ener": r"-",
                   "n-K": r".",
                   }
    if 'run-rlk-' in filename:
        return 1

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
                            except ValueError:
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

# Given a list of files and a list of replica numbers find which files match the list of reps
def files_with_correct_reps(files, reps):
#    try:  # The try except is for SH compatibility
    rep_nums = {f:find_rep_num_in_name(f[f.rfind('/'):]) for f in files}
#    except:
#        print("""Replica sorting isn't working for this data, switching it off.
#
#This means all data files will be read.
#
#If this is Surface Hopping this is OK as it does not output different trajectory
#info in the same format as CTMQC.
#
#
#
#""")
#        return files
    files = [i[1] for i in sorted(zip(rep_nums.values(), rep_nums.keys()))]
    if type(reps) == int:
        reps = [reps]
    elif 'all' in reps:
        return files
    files = [f for f in files if rep_nums[f] in reps]
    return files


# A helper function for parallelising the loading of the files
def helperLoad(args):
    func, args = args
    return func(*args)


# Will apply a function to load all relevant files in a folder
def load_all_in_folder(folder, func, args=[], filename_must_not_contain=[],
                       filename_must_contain=[], reps='all'):
    """
    Will repeat a file loading function for all files in a folder depending on
    what rules for finding the files are specified in the arguments.

    Inputs:
        * folder  =>  folder to look in for files [str]
        * func    =>  func to apply to each file  [func]
        * args    =>  arguments to be applied to the file loading function in
                      order [list]
        * filename_must_not_contain => any strs the filenames must not contain
        * filename_must_contain     => any strs the filenames must contain
        * reps    => Which reps to load (can be 'all' or list of ints)

    Outputs:
        * A dictionary with filenames as keys and the return of the func passed
          as values.
    """
    if not os.path.isdir(folder):
            raise SystemError("Sorry the folder given can't be found...")
    else:
        if folder[-1] != '/':
            folder = folder + '/'
    filename_must_not_contain.append(".sw")
    all_files1 = [folder+i for i in os.listdir(folder)
                  if all(j in i for j in filename_must_contain) and
                  all(k not in i for k in filename_must_not_contain)]
    # Error checking
    if not all_files1:
        mustntContMsg = ", ".join(filename_must_not_contain)
        msg = "\n\tSorry I can't find any files with the correct filenames!\n"
        msg += "\n\n\t\tmust_contain = %s" % (", ".join(filename_must_contain))
        msg += "\n\t\tmustn't contain = %s" % (mustntContMsg)
        msg += "\n\n\tFiles found = %s\n\n" % (', '.join(all_files1))
        print(msg)

    # Only read files with the correct rep num
    all_files2 = files_with_correct_reps(all_files1, reps)
    if not all_files2:
        valid_rep_nums = [find_rep_num_in_name(f) for f in all_files1]
        corrReps = ",\t".join([str(i) for i in valid_rep_nums])
        msg = "\n\n\n\tSorry I can't find any files with the correct replica "
        msg += "number!\n\n\n\t"
        msg += "Valid replica numbers are:\n\t\t* %s" % corrReps
        print(msg)

    # The dictionary stores the data
    all_data = collections.OrderedDict()
    if len(all_files2) > 1:
      doParallel = True
    else:
      doParallel = False
    if doParallel:  # Read files in parallel
        args = [(func, [f]+list(args)) for f in all_files2]
        numProc = mp.cpu_count()
        if (numProc > 20):
            numProc = 20
        if (numProc > len(args)):
            numProc = len(args)
        p = mp.Pool(numProc)
        res = p.map(helperLoad, args)
        p.close()
        p.join()

        for arg, result in zip(args, res):
            fName = arg[1][0][arg[1][0].rfind('/')+1:]
            all_data[fName] = result
    else:  # Read files in serial
        args = [[f]+list(args) for f in all_files2]
        for arg in args:
            fName = arg[0][arg[0].rfind('/')+1:]
            print(fName)
            all_data[fName] = func(*arg)

    if not all_data:
        return False
    return all_data


def reshape_by_state(data, cols):
   """
   Will reshape the ad_frc and ad_mom arrays so that the arrays are
   in a (Nsteps, Nstates, Natoms, 3) order.

   Inputs:
      * dataDict => the dictionary with all the data.

   Outputs:
      * the dataDictionary that contains the reshaped data.
   """
   # Get vital metadata
   oldShape = data.shape
   numStates = max(cols[0, :, 1].astype(int))
   numAtoms = oldShape[1]/numStates
   if int(numAtoms) != numAtoms:
      msg = "Something went wrong calculating the new shape for the"
      msg += " ad frc array"
      raise SystemExit(msg)
   numAtoms = int(numAtoms)
   
   # Create new array
   newShape = (oldShape[0], numStates, numAtoms, oldShape[2])
   #newShape = (1, numStates, numAtoms, oldShape[2])
   newData = np.zeros(newShape)
   
   # Populate new array
   for istate in range(1, numStates+1):
      stateMask = cols[0, :, 1] == str(istate)
      for istep in range(oldShape[0]):
         newData[istep][istate-1] = data[istep][stateMask]

   return newData

