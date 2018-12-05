#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  7 15:42:51 2018

@author: mellis
"""
import numpy as np
import pandas as pd
import collections
import itertools as IT

from load import load_coeff


def calc_U_matrix(H):
    """
    Calculates the U matrix for each replica according to the hamiltonian.
    
    Inputs:
        * all_ham_data  => all the hamiltonian data in a dictionary with 
                           filenames as keys
    """
    H = -H
    
    E, U = np.linalg.eigh(H) #get eigenvectors

#    # Sort eigenvectors by eigenvalues        
#    idx = E.argsort() #Sort by absolute energy
#    E = E[idx]
#    U = U[:,idx]
#    
    
#    U_T_1 = np.matrix(np.linalg.inv(U).T)
#    # To check the unitary transformation matrix
#    if np.max(U_T_1 - U).real > 1e-14:
#        raise SystemExit("""Something is wrong with the transformation matrix U
#
#Max deviation from U^{\dagger} - U = %.2g"""%(np.max(np.transpose(np.linalg.inv(U)) - U).real))
    return E, U

# Will transform the diabatic coefficients with the hamiltonian
def transform_di_2_ad(udata, ucols, utimesteps, hdata, htimesteps):
    hmask = [h in utimesteps for h in htimesteps]
    umask = [u in htimesteps for u in utimesteps]
    
    udata = udata[umask]
    hdata = hdata[hmask]
    ucols = ucols[umask]
    utimesteps = utimesteps[umask]
    
    Us = np.array([calc_U_matrix(H)[1] for H in hdata])
    
    C = np.array([np.dot(U.T, u) for U, u in zip(Us, udata)])
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

def avg_pos_data(all_pos_data):
    """
    Will average the postition data.
    
    Inputs:
        * all_pos_data => dictionary containing all the position data (filenames as keys)
    
    Outputs:
        * avg position data (same format dictionary as input)
    """
    min_len = np.min([[len(all_pos_data[f][0][0]), len(all_pos_data[f][0][1]), len(all_pos_data[f][1])] for f in all_pos_data])
    all_pos = [[all_pos_data[f][0][0][:min_len] for f in all_pos_data ],
                '',
                [all_pos_data[f][1][:min_len] for f in all_pos_data ]]
    avg_pos = {}
    if all_pos[0]:
        # Should return in same format as input
        avg_pos = {'avg_pos':[(np.mean(all_pos[0], axis=0),
                               np.zeros(len(all_pos[0][0]))), 
                       np.mean(all_pos[2], axis=0)]}   
    return avg_pos

def avg_Qlk_data(all_Qlk_data):
    """
    Will average the Qlk data and return a dictionary with the data.
    
    Inputs:
        * all_Qlk_data  =>  A dictionary containing all the data to be averaged
    
    Outputs:
        The averaged Qlk_data in the same format as the inputted dictionary.
    """
    avg_Qlk_data = np.mean([all_Qlk_data[rep][0][0] for rep in all_Qlk_data], axis=0)
    avg_Qlk_timesteps = np.mean([all_Qlk_data[rep][1] for rep in all_Qlk_data], axis=0)
    
    avg_Qlk_data = np.array(avg_Qlk_data)
    avg_Qlk_timesteps = np.array(avg_Qlk_timesteps)
    
    Qkeys = list(all_Qlk_data.keys())
    avg_Qlk_data = {'avg_Qlk':
                              [ [avg_Qlk_data,
                                 all_Qlk_data[Qkeys[0]][0][1]],
                               avg_Qlk_timesteps]}
    return avg_Qlk_data

def sum_hist_f_data(all_tintf_data):
    """
    Will sum all the time-integrated adiabatic forces.
    
    Inputs:
        * all_tintf_data    =>  A dictionary containing all the history force 
                                data.
    
    Ouputs:
        * The averaged tintf data.
    """
    Tcols = all_tintf_data[list(all_tintf_data.keys())[0]][0][1]
    nstates = max(Tcols[0,:,1].astype(int))
    all_state_combs = list(IT.combinations(range(nstates), 2))
    all_tintf_lk = {"%i%i"%(i,j):'' for i,j in all_state_combs}

    Tdata = all_tintf_data[list(all_tintf_data.keys())[0]][0][0]
    nsteps = len(Tdata)
    natom = int(len(Tdata[0])/nstates)
    for i,j in all_state_combs:
        all_tintf_lk["%i%i"%(i,j)] = np.zeros((nsteps, natom, 3))
    
    #Actually do the calc
    for Tkey in all_tintf_data:
        (Tdata, Tcols), _ = all_tintf_data[Tkey]

        for l,k in all_state_combs:
            fk = np.array([i[Tcols[0,:,1] == str(k+1)] for i in Tdata])
            fl = np.array([i[Tcols[0,:,1] == str(l+1)] for i in Tdata])
            
            tmp = [F for i, F in enumerate(fl-fk)]            
            all_tintf_lk['%i%i'%(l,k)] += tmp
        
    Tsteps = all_tintf_data[list(all_tintf_data.keys())[0]][1]
    return [all_tintf_lk, Tsteps]


def sum_hist_f_CC_data(all_tintf_data, all_A_coeff_data):
    """
    Will sum all the time-integrated adiabatic forces with the adiabatic 
    populations.
    
    Inputs:
        * all_tintf_data    =>  A dictionary containing all the history force 
                                data.
    
    Ouputs:
        * The averaged tintf data.
    """
    Tcols = all_tintf_data[list(all_tintf_data.keys())[0]][0][1]
    nstates = max(Tcols[0,:,1].astype(int))
    all_state_combs = list(IT.combinations(range(nstates), 2))
    all_tintf_lk = {"%i%i"%(i,j):'' for i,j in all_state_combs}

    #Match the timesteps for the 2 sets of data
    for Akey, Tkey in zip(all_A_coeff_data, all_tintf_data):
        all_A_coeff_data[Akey], all_tintf_data[Tkey] = match_timesteps_4col_2col( all_A_coeff_data[Akey],
                                                                            all_tintf_data[Tkey])

    Tdata = all_tintf_data[list(all_tintf_data.keys())[0]][0][0]
    nsteps = len(Tdata)
    natom = int(len(Tdata[0])/nstates)
    for i,j in all_state_combs:
        all_tintf_lk["%i%i"%(i,j)] = np.zeros((nsteps, natom, 3))
    
    #Actually do the calc
    for Akey, Tkey in zip(all_A_coeff_data, all_tintf_data):
        (Tdata, Tcols), _ = all_tintf_data[Tkey]
        _, _, _, pops  = all_A_coeff_data[Akey]

        for l,k in all_state_combs:
            fk = np.array([i[Tcols[0,:,1] == str(k+1)] for i in Tdata])
            fl = np.array([i[Tcols[0,:,1] == str(l+1)] for i in Tdata])
            Clk = pops[:,l] * pops[:,k]
            
            tmp = [Clk[i] * F for i, F in enumerate(fl-fk)]            
            all_tintf_lk['%i%i'%(l,k)] += tmp
        
    Tsteps = all_tintf_data[list(all_tintf_data.keys())[0]][1]
    return [all_tintf_lk, Tsteps]

def match_timesteps(DCT1, DCT2):
    """
    Will return data with the same timesteps.
    
    Inputs:
        * DCT1  =>  (data1, cols1, timesteps1)  [list]
        * DCT2  =>  (data2, cols2, timesteps2)  [list]
    
    Outputs:
        * DCT1 and DCT2 (same as inputs but spliced)
    """
    # Make sure that all data in data1 is in data2
    mask1 = [i in DCT2[2] for i in DCT1[2]]
    DCT1[0] = DCT1[0][mask1]
    DCT1[1] = DCT1[1][mask1]
    DCT1[2] = DCT1[2][mask1]
    
    #Make sure any data in data2 is in data1    
    mask2 = [i in DCT1[2] for i in DCT2[2]]
    DCT2[0] = DCT2[0][mask2]
    DCT2[1] = DCT2[1][mask2]
    DCT2[2] = DCT2[2][mask2]
    
    return DCT1, DCT2, mask1, mask2

# There is a method in plot_QM.QM_R that also matches data. It probs belongs 
# here.
def match_timesteps_4col_2col(data_4col, data_2col):
    """
    Will return 2 sets of data which share the same timesteps
    """
    (data2, cols2), timesteps2 = data_2col
    data4, cols4, timesteps4, populations4   = data_4col
    DCT1, DCT2, mask1, mask2 = match_timesteps([data2, cols2, timesteps2], 
                                               [data4, cols4, timesteps4])
    populations4 = populations4[mask2]
    
    return [data4, cols4, timesteps4, populations4], [(data2, cols2), timesteps2]
    

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
def load_Acoeff_data(folder, reps, all_ham_data, max_step='all', min_step=0, stride=1):
    all_coeff_data = load_coeff.load_all_coeff_in_folder(folder, 
                                                         filename_must_contain=['ad','xyz','coeff'], 
                                                         reps=reps,
                                                         max_step=max_step, 
                                                         min_step=min_step, 
                                                         stride=stride)
    # Transform Diabatic
    if not all_coeff_data:
        all_Dcoeff_data = load_coeff.load_all_coeff_in_folder(folder, 
                                                              filename_must_contain=['coeff','xyz'], 
                                                              filename_must_not_contain=['ad'], 
                                                              reps=reps,
                                                              max_step=max_step, 
                                                              min_step=min_step, 
                                                              stride=stride)
        
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


def print_timings(timings_dict, ntabs=0, max_len=50):
    """
    Will print timing data in a pretty way.
    
    Inputs:
        * timings_dict  =>  A dictionary containing all the timings data
        * ntabs         =>  OPTIONAL (not recommended to change). Number of 
                            tabs in the printout.
    Ouputs:
        None
    """
    def print_line(line):
        line = "|"+line+" "*(max_len-len(line))+"|"
        print(line)
    
    tab = "    "
    for Tkey in timings_dict:
        if type(timings_dict[Tkey]) == dict or type(timings_dict[Tkey]) == collections.OrderedDict:
            print("|" + " "*max_len + "|")
            line = "%s"%(Tkey)
            print_line(line)
            print_timings(timings_dict[Tkey], ntabs+2)
        else:
            line = "%s* %s:"%(tab*ntabs, Tkey)
            str_num = "%.0e s"%timings_dict[Tkey]
            line = line + " "*(max_len-26 - (len(line) + len(str_num)) + ntabs*5) + str_num
            print_line(line)
    return 
