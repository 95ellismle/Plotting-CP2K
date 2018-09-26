#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  7 15:42:51 2018

@author: mellis
"""
import numpy as np
import pandas as pd

from load import load_coeff


# Will transform the diabatic coefficients with the hamiltonian
def transform_di_2_ad(udata, ucols, utimesteps, hdata, htimesteps):
    hmask = [h in utimesteps for h in htimesteps]
    umask = [u in htimesteps for u in utimesteps]
    
    udata = udata[umask]
    hdata = hdata[hmask]
    ucols = ucols[umask]
    utimesteps = utimesteps[umask]
    
    Us = np.linalg.eigh(hdata)[1]
    C = np.array([np.dot(U, u) for U, u in zip(Us, udata)])
    pops = np.array([np.array([np.linalg.norm(j)**2 for j in i]) for i in C])
    return C, ucols, utimesteps, pops

# Transforms all the diabatic coefficients to adiabatic and returns the dict
def trans_all_diab_to_adiab(all_Dcoeff_data, all_ham_data, reps):
    print("Sorry I can't find any adiabatic coefficient data, I'll try and transform the diabatic ones")
    all_Acoeff_data = {fC.replace("coeff", "coeff_ad"):
                        transform_di_2_ad(all_Dcoeff_data[fC][0],
                                     all_Dcoeff_data[fC][1],
                                     all_Dcoeff_data[fC][2],
                                     all_ham_data[fH][0],
                                     all_ham_data[fH][2])
                        for fC, fH in zip(all_Dcoeff_data, all_ham_data)}
    return all_Acoeff_data

# Will average the coeff data dictionnary
def avg_C_data_dict(all_coeff_data):
    min_len = np.min([[len(all_coeff_data[f][i]) for f in all_coeff_data] for i in range(4)])
    all_ad =  [[all_coeff_data[f][0][:min_len] for f in all_coeff_data if 'ad' in f],
                    '',
                    [all_coeff_data[f][2][:min_len] for f in all_coeff_data if 'ad' in f],
                    [all_coeff_data[f][3][:min_len] for f in all_coeff_data if 'ad' in f]]
    ad_avg = {}
    if all_ad[0]:
        ad_avg = {'avg_ad':[np.mean(all_ad[0], axis=0),
                       '', 
                       np.mean(all_ad[2], axis=0), 
                       np.mean(all_ad[3], axis=0)]}    
    
    all_di = [[all_coeff_data[f][0][:min_len] for f in all_coeff_data if 'ad' not in f],
                    '',
                    [all_coeff_data[f][2][:min_len] for f in all_coeff_data if 'ad' not in f],
                    [all_coeff_data[f][3][:min_len] for f in all_coeff_data if 'ad' not in f]]
    di_avg = {}
    if all_di[0]:
        di_avg = {'avg_di':[np.mean(all_di[0], axis=0),
                       '', 
                       np.mean(all_di[2], axis=0), 
                       np.mean(all_di[3], axis=0)]}
    return ad_avg, di_avg

# Will average the coeff data dictionnary
def avg_coeff_data(all_coeff_data):
    #Find the max length of data that all arrays can give
    max_len = np.min([[len(all_coeff_data[f][i]) for f in all_coeff_data] for i in range(4)])
    #Splice to the max length
    all_coeff =  [[all_coeff_data[f][0][:max_len] for f in all_coeff_data],
                    '',
                    [all_coeff_data[f][2][:max_len] for f in all_coeff_data],
                    [all_coeff_data[f][3][:max_len] for f in all_coeff_data]]
    
    #Now average
    if all_coeff[0]:
        coeff_avg = [np.mean(all_coeff[0], axis=0),
                       '', 
                       np.mean(all_coeff[2], axis=0), 
                       np.mean(all_coeff[3], axis=0)]  

    return coeff_avg

# Will average the hamiltonians
def avg_H_data_dict(all_ham_data):
    min_len = np.min([[len(all_ham_data[f][i]) for f in all_ham_data] for i in range(3)])
    all_hams = [[all_ham_data[f][0][:min_len] for f in all_ham_data ],
                '',
                [all_ham_data[f][2][:min_len] for f in all_ham_data ]]
    avg_hams = {}
    if all_hams[0]:
        avg_hams = {'avg_ham':[np.mean(all_hams[0], axis=0),
                       '', 
                       np.mean(all_hams[2], axis=0)]}   
    return avg_hams

# Will average the energy data
def avg_E_data_dict(all_ener_data):
    cols_dict = {col:[] for col in all_ener_data[list(all_ener_data.keys())[0]]}
    for filename in all_ener_data:
        for col in cols_dict:
            cols_dict[col].append(all_ener_data[filename][col])
    return pd.DataFrame({col:np.mean(cols_dict[col], axis=0) for col in cols_dict})

# Will get the number of replicas and check if the user is asking to plot loads of graphs
def get_num_reps(all_data, avg_on, found_num_reps):
        num_reps = len(all_data)
        if num_reps > 5 and not avg_on and not found_num_reps:
            Q = raw_input("Are you sure you want to plot each one of %i replicas seperately?\n\nI can average them if you set that in the settings file.\n\t(yes|no):\t"%(num_reps))
            if 'y' in Q.lower():
                return num_reps, True
            else:
                raise SystemExit("OK, I am exiting so you can average the replica data.\n\nTo average the replicas set the variable 'avg_reps = True'")
        return num_reps, True
#Will get the avg_coupling, site_ener_diff 
def get_coup_data(H_data, key):
    site_ener = H_data[key][0][:,0,0]-H_data[key][0][:,1,1]
    couplings = H_data[key][0][:,0,1]
    avg_couplings = np.mean(couplings)
    timesteps = H_data[key][2]
    
    return site_ener, couplings, avg_couplings, timesteps


# Will load the Adiab coeff data or transform it from the diabatic coefficents.
def load_Acoeff_data(folder, reps, all_ham_data):
    all_coeff_data = load_coeff.load_all_coeff_in_folder(folder, filename_must_contain=['ad','xyz','coeff'], reps=reps)
    # Transform Diabatic
    if not all_coeff_data:
        all_Dcoeff_data = load_coeff.load_all_coeff_in_folder(folder, filename_must_contain=['coeff','xyz'], filename_must_not_contain=['ad'], reps=reps)
        all_Acoeff_data = trans_all_diab_to_adiab(all_Dcoeff_data=all_Dcoeff_data,
                                                             all_ham_data = all_ham_data, 
                                                             reps=reps)
        if all_Acoeff_data:
            print("Sucessfully transformed the diab coefficients!")
        else:
            raise SystemExit("Sorry I couldn't transform the diabatic coefficients. Something went wrong")
    else:
        all_Acoeff_data = all_coeff_data
    return all_Acoeff_data