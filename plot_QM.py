#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 15 13:51:30 2018

@author: mellis
"""

from load import load_QM

from IO import Folders as fold

#from Plot import Plot as plt
#from Plot import plot_utils



"""
Loading works well, plotting the values still needs work. Need to put labels on graph, sort out averaging the data etc...
"""




folder = fold.make_fold_abs('/scratch/mellis/flavoured-cptk/200Rep_3mol') 
params = {
    'at_num'    : 'all',        # Which atom number to plot
    'lk'        : (2,1),            # Which state, l, to plot
    'cart_dim'  : 2,                # Which Cartesian Dimension to plot
}
reps      = [1]             # Which replicas to plot
avg_reps  = False          # Whether or not to average the values between replicas


import matplotlib.pyplot as plt



all_Qlk_data = load_QM.load_all_Qlk_in_folder(folder, reps=reps)
num_reps = len(all_Qlk_data)
for irep, filename in enumerate(all_Qlk_data):  
    Qlk_data = all_Qlk_data[filename][0]
    
    # Convert any 'all's into a list
    for P in params:
        if type(params[P]) == str and 'all' in params[P].lower():
            if P == 'lk':
                all_states = load_QM._get_all_Qlk_col_vals_(Qlk_data[1], 'l')
                params[P] = []
                for l in all_states:
                    for k in range(1,l):
                        params[P].append((l,k))
            else:
                params[P] = load_QM._get_all_Qlk_col_vals_(Qlk_data[1], P)
            
    # Plot any QMs
    if any(type(params[p]) == list for p in params):
        for P in params:
            if type(params[P]) == list:
                f,a = plt.subplots()
                for param in params[P]:
                    new_params = params.copy()
                    new_params[P] = param
                    data = load_QM.find_in_Qlk(Qlk_data, new_params) 
                    a.plot(all_Qlk_data[filename][1], data, label="%s = %s"%(P, param))
                a.legend(fontsize=19)
    else:
        f,a = plt.subplots()
        data = load_QM.find_in_Qlk(Qlk_data, params)
        a.plot(all_Qlk_data[filename][1], data, avg_reps, num_reps)
        a.set_ylim([min(data),max(data)])


    
        
        
    