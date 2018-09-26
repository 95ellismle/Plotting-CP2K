import matplotlib.pyplot as plt
import numpy as np


#Will plot with stylings determined by the style dict
def plot_style(axis, X, Y, style):
    """
    Will plot the X and Y data dependent on the style in the dicitonary style.
    """
    #Set defaults
#    if style.get('color') == None: style['color'] = 'k'
    if style.get('label') == None: style['label'] = ''
    if style.get('ls')    == None: style['ls'] = '-'
    if style.get('alpha') == None: style['alpha'] = 1.0
    if style.get('lw')    == None: style['lw'] = 1
    if style.get('xlab')  == None: style['xlab'] = ''
    if style.get('ylab')  == None: style['ylab'] = ''
    if style.get('title') == None: style['title'] = ''
    
    # Plot with styles
    l, = axis.plot(X,Y, ls=style['ls'], lw=style['lw'], alpha=style['alpha'], label=style['label'])
    
    # Add labels
    axis.set_ylabel(style['ylab'])
    axis.set_xlabel(style['xlab'])
    axis.set_title(style['title'])
    
    if style.get("color"): l.set_color(style['color'])
    if style.get('xlims'): axis.set_xlim(style['xlims'])

# Will plot the populations of the expansion coefficients
def plot_pops(all_coeff_data, axis, ad_or_di, style_dict={}):
    
    data, cols, timesteps, pops = all_coeff_data
    
    for ibasis in range(len(pops[0])):
        if style_dict.get('label') == None: style_dict['label'] = "State %i"%(ibasis+1)
        plot_style(axis, timesteps, pops[:,ibasis], style=style_dict)
    axis.legend()
    style_dict.pop('ylab')


# Will plot the site_energy, diabatic populations and adiabatic populations in one graph stacked
def plot_site_ener(site_ener, axis, style_dict={}):
    
    Stimesteps = site_ener[1]
    site_ener  = site_ener[0]
    #Plot 1: Site Ener
    if style_dict.get('title') == None: style_dict['title'] = r"Site Energy Difference vs Time"
    if not style_dict.get('ylab'): style_dict['ylab'] = r"$\Delta$E (Ha)"
    plot_style(axis, Stimesteps, site_ener, style_dict)
    style_dict.pop('ylab')



# Gets the norms of the pops
def plot_norms(all_coeff_data, axis, style_dict={}):
    coeffs, cols, timesteps, pops = all_coeff_data
    norms = np.sum(pops, axis=1)
    
    plot_style(axis, timesteps, norms, style_dict)
    
    axis.set_ylabel(r"$\sum_{l} |$u$_{l}|^2$ ")
    
#    min_val, max_val = np.min([0,min([min(line.get_ydata()) for line in axis.lines])]), np.max([1,max([max(line.get_ydata()) for line in axis.lines])])
#    axis.set_ylim([min_val, max_val])
    
    axis.set_title("Diabatic Norm Conservation")
    
#    fit = np.polyfit(timesteps, norms, 1)
#    axis.annotate(r"Norm Drift: $\frac{d}{dt}$[$\sum_{l} |$u$_{l}|^2$] = %.2gt + %.2g"%(fit[0], fit[1]), ( (np.max(timesteps)-np.min(timesteps))/999., 1-(np.max(norms)-np.min(norms))/1.5), fontsize=24)
    axis.axhline(1, color='k', ls='--')
    axis.get_yaxis().get_major_formatter().set_useOffset(False)

#Decides what arrangement of axes to use
def axes_arrangement(params):
    if len(params)   == 1:        
        f,a = plt.subplots(1)
        f,a = [f], [a]
    elif len(params) <= 3:        f,a = plt.subplots(len(params))
    elif len(params) == 4:        f,a = plt.subplots(2,2)
    else: raise SystemExit("Sorry no rule has been set for a arrangement of more than 4 axes!")
    return f, a

# Will plot specified parameters
def plot_params(all_Dcoeff_data, all_Acoeff_data, site_ener, params, style_dict={}, FA=False):
    """
    Will plot the graphs of the parameters specified in params. 
    """
    if FA == False:    f, a = axes_arrangement(params)
    else: f, a = FA
    
    max_ax = len(params)
    if max_ax == 4: max_ax = 2    

    for i,param in enumerate(params):
        if param == 'norm':
            plot_norms(all_Dcoeff_data, a[i], style_dict=style_dict)
        elif param == "|u|^2":
            if not style_dict.get('title') and i ==0 : style_dict['title'] = 'Diabatic Population vs Time'
            if not style_dict.get('ylab'): style_dict['ylab'] = r'$|u_l|^2$'
            plot_pops(all_Dcoeff_data, a[i], 'Diabatic', style_dict=style_dict)
        elif param == "|c|^2":
            if not style_dict.get('title') and i ==0 : style_dict['title'] = 'Adiabatic Population vs Time'
            if not style_dict.get('ylab'): style_dict['ylab'] = r'$|C_l|^2$'
            plot_pops(all_Acoeff_data, a[i], 'Aiabatic', style_dict=style_dict)      
        elif param == "site_ener":
            plot_site_ener(site_ener, a[i], style_dict=style_dict)
    a[i].set_xlabel("Time (fs)")
#    plt.tight_layout()



# Will close all windows
def close_all():
    for i in range(100):
        plt.close()