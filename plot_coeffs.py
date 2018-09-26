#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  7 12:48:21 2018

@author: mellis
"""

from load import load_coeff
from load import load_ham

from IO import Folders as fold

from Plot import Plot as plt
from Plot import plot_utils

#import numpy as np

"""
Specify in the plot_params which quantities you would like plotting. 

Plot params should be a list of strings.

Quantities to plot are:
    * 'site_ener'     = site_energies
    * 'norm'          = norm (of diabatic coefficients)
    * '|u|^2'         = population of diabatic coeffficients
    * '|C|^2'         = population of adiabatic coefficients
"""


plot_params    = ['norm']  
reps = 'all'                              #Which replicas to plot
avg_reps = True                           #Whether or not to average the values between replicas
folder = fold.make_fold_abs('/scratch/mellis/flavoured-cptk/200Rep_3mol') #The folder to look in for the data





"""Remove some global variables"""
for var in ["all_Dcoeff_data", "all_Acoeff_data", 'all_ham_data', 'all_coeff_data']:
    if var in globals():
        del globals()[var]



"""Setting some vars"""
found_num_reps = False


"""Load Coupling Data"""
all_ham_data = load_ham.load_all_ham_in_folder(folder, reps=reps)
avg_ham_data = plot_utils.avg_H_data_dict(all_ham_data)
avg_site_ener, avg_couplings, avg_avg_couplings, timesteps = plot_utils.get_coup_data(avg_ham_data, 'avg_ham')
all_site_ener = [plot_utils.get_coup_data(all_ham_data, ham_key) for ham_key in all_ham_data]
all_site_ener = [[i[0], i[3]] for i in all_site_ener]



######## SHOULD MOVE THIS TO THE LOAD_COEFF.py module #########
def load_all_di_coeffs():
    """ 
    Loads all the diabatic coefficients, no input. Returns a dict of coeffs with filenames as keys.
    """
    if 'all_Dcoeff_data' not in globals():
        all_Dcoeff_data = load_coeff.load_all_coeff_in_folder(folder, filename_must_contain=['xyz','coeff'], filename_must_not_contain=['ad'], reps=reps)
    else:
        all_Dcoeff_data = globals()['all_Dcoeff_data']
        if 'avg_di' in all_Dcoeff_data.keys():
            all_Dcoeff_data = load_coeff.load_all_coeff_in_folder(folder, filename_must_contain=['xyz','coeff'], filename_must_not_contain=['ad'], reps=reps)
    return all_Dcoeff_data    
    
def load_all_ad_coeffs():     
    """ 
    Loads all the adiabatic coefficients, no input. Returns a dict of coeffs with filenames as keys.
    """
    if "all_Acoeff_data" not in globals(): 
        all_Acoeff_data = plot_utils.load_Acoeff_data(folder, reps, all_ham_data)
    else:
        all_Acoeff_data = globals()['all_Acoeff_data']
        if 'avg_ad' in all_Acoeff_data.keys():
            all_Acoeff_data = plot_utils.load_Acoeff_data(folder, reps, all_ham_data)
    return all_Acoeff_data
        


"""
Plotting 3 specified parameters.
"""
if plot_params:
    plot_params = [i.lower().strip() for i in plot_params]
    diab_watchwords = ['|u|^2','norm']
    if any(i in j for i in diab_watchwords for j in plot_params): all_Dcoeff_data = load_all_di_coeffs()
    else:                                            all_Dcoeff_data = {i:'' for i in range(3000)}
    if '|c|^2' in plot_params:        all_Acoeff_data = load_all_ad_coeffs()
    else:                          all_Acoeff_data = {i:'' for i in range(3000)}
    
    
    FA = plt.axes_arrangement(plot_params)
    style_dict = {}
    if avg_reps:
        
        #Plot the individual replicas faded
        for Dfilename, Afilename, site_ener in zip(all_Dcoeff_data, all_Acoeff_data, all_site_ener):
            style_dict['alpha'] = 0.1
            style_dict['label'] = ''
            style_dict['color'] = 'g'
            plt.plot_params(all_Dcoeff_data[Dfilename], all_Acoeff_data[Afilename], site_ener, plot_params, style_dict=style_dict, FA=FA)
        style_dict.pop('color')
        style_dict['lw']=2
        style_dict.pop('alpha')
        
        
        #Average the data here and put it in a list for plot_params func
        if all_Dcoeff_data.keys()[0] != 0:
            _, all_Dcoeff_data = plot_utils.avg_C_data_dict(all_Dcoeff_data)
        if all_Acoeff_data.keys()[0] != 0:
            all_Dcoeff_data, _ = plot_utils.avg_C_data_dict(all_Acoeff_data)
            
    
    #Plot the averages (or replicas in different graphs)
    for Dfilename, Afilename, site_ener in zip(all_Dcoeff_data, all_Acoeff_data, all_site_ener):
        plt.plot_params(all_Dcoeff_data[Dfilename], all_Acoeff_data[Afilename], site_ener, plot_params, style_dict=style_dict, FA=FA)
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
#'''
#Plotting the adiabatic coefficients (if they haven't been printed off 
#then transform them using the hamiltonian)
#'''
#if plot_adiabatic:
#    
#    all_Acoeff_data = plot_utils.load_Acoeff_data(folder, reps, all_ham_data)
#    num_reps, found_num_reps = plot_utils.get_num_reps(all_Acoeff_data, avg_reps, found_num_reps)
#    
#    if avg_reps:
#        all_Acoeff_data, _ = plot_utils.avg_C_data_dict(all_Acoeff_data)
#    
#    if all_Acoeff_data:
#        print("Plotting Adiabatic Coeffs")
#        Afigs, Aaxs = plt.plot_pops(all_Acoeff_data, num_reps=num_reps)
#
#
#
#
#'''
#Plotting the diabatic coefficients
#'''
#if plot_diabatic:
#    all_Dcoeff_data = load_all_di_coeffs()
#    num_reps, found_num_reps = plot_utils.get_num_reps(all_Dcoeff_data, avg_reps, found_num_reps)
#
#    if avg_reps:
#        _, all_Dcoeff_data = plot_utils.avg_C_data_dict(all_Dcoeff_data)
#    
#    if all_Dcoeff_data:
#        print("Plotting Diabatic Coeffs")
#        Dfigs, Daxs = plt.plot_pops(all_Dcoeff_data, num_reps=num_reps)
#
#
#
#
#
#'''
#Plotting the site energy, adiabatic and diabatic coefficient graph
#'''
#if plot_site_ener:
#    
#    all_Dcoeff_data = load_all_di_coeffs()
#    all_Acoeff_data = load_all_ad_coeffs()
#    
#    num_reps, found_num_reps = plot_utils.get_num_reps(all_Acoeff_data, avg_reps, found_num_reps)
#    
#    if avg_reps:
#        all_coeff_data = all_Dcoeff_data.copy() #Merging 2 dicts
#        all_coeff_data.update(all_Acoeff_data)
#        all_Acoeff_data, all_Dcoeff_data = plot_utils.avg_C_data_dict(all_coeff_data) #Averaging
#        
#        plt.plot_site_ener_diab_pops_adiab_pops(all_Dcoeff_data['avg_di'][2], all_Dcoeff_data['avg_di'][3],
#                                     all_Acoeff_data['avg_ad'][2], all_Acoeff_data['avg_ad'][3],
#                                     avg_ham_data['avg_ham'][2], avg_site_ener,
#                                     num_reps = num_reps) 
#    
#    else:
#        for Dfile, Afile, site_ener, Hfile in zip(all_Dcoeff_data, all_Acoeff_data, all_site_ener, all_ham_data):
#            plt.plot_site_ener_diab_pops_adiab_pops(all_Dcoeff_data[Dfile][2], all_Dcoeff_data[Dfile][3],
#                             all_Acoeff_data[Afile][2], all_Acoeff_data[Afile][3],
#                             all_ham_data[Hfile][2], site_ener,
#                             num_reps = 1) 
#            
#    
#    
#
#
#'''
#Plotting the norm conservation for the adiabatic and diabatic coefficients
#'''
#if plot_norm_cons:
#    
#    all_Dcoeff_data = load_all_di_coeffs()
##    all_Acoeff_data = load_all_ad_coeffs()
#    Dfigs, Daxs = plt.plot_norms(all_Dcoeff_data, ' for CTMQC')
##    Afigs, Aaxs = plt.plot_norms(all_Acoeff_data, ' for CTMQC')
#    