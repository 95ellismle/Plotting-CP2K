#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 17 11:38:10 2018

@author: mellis
"""

from load import load_ener

from IO import Folders as fold

#from Plot import Plot as plt
from Plot import plot_utils


#Global settings
plot_ad_ener = True                        #Plot only the adiabatic coefficients
plot_tot_ener  = False                     #Plot only the diabatic coefficients
reps = [1]                               #Which replicas to plot
avg_reps = False                            #Whether or not to average the values between replicas


#Settings local to particular plots
#   ad_ener
states = 'all'



#The folder to look in for the data
folder = fold.make_fold_abs('/scratch/mellis/flavoured-cptk/200Rep_2mol') 

if plot_ad_ener:
    all_ad_ener_data = load_ener.load_all_ener_ad(folder, reps=reps)
    if not all_ad_ener_data:
        raise IOError("Can't find any data, please check folder.")
    if avg_reps:
        all_ad_ener_data = plot_utils.avg_E_data_dict(all_ad_ener_data)
    
    import matplotlib.pyplot as plt
    
    ener_data = all_ad_ener_data
    
    if states == 'all':
        states = [int(i.lower().replace('state_','')) for i in ener_data.columns if 'state' in i.lower()]
    
    
    # Plots a single graph
    f,a = plt.subplots()
    for l in states:
        a.plot(ener_data['Time'], ener_data["State_%i"%l], label="State %i"%l)
    a.legend()