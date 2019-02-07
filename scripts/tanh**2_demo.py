#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 10:53:29 2018

@author: mellis
"""

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import numpy as np

tanh_width = 0.001


fig, ax = plt.subplots()
#plt.subplots_adjust(bottom=0.25)
ax_slider = plt.axes([0.25, 0.1, 0.65, 0.03])

x = np.array(list(np.arange(-0.1,-0.001,0.0001)) + list(np.arange(0.001,0.1,0.0001)))
y = 1/(x)
ax.plot(x,y/200, label=r'$\frac{1}{x}$')


x1 = np.arange(-0.1,0.1,0.000012)
y1 = np.tanh(x1/tanh_width)**2
ln1, = ax.plot(x1,y1)



x2 = np.arange(-0.1,0.1,0.000012)
y2 = (1/x2) * np.tanh(x2/tanh_width)**2
ln2,  = ax.plot(x2,y2/200)


#[np.tanh(1.08866)**2/(1.08866*tanh_width)]
maxX = 1.08866*tanh_width
ln4, = ax.plot([maxX],
               (1/200) * (0.5826/tanh_width), 'ko')

ann = ax.annotate("Width = %.2g"%tanh_width, (-0.11,3), fontsize=22)


twidth = Slider(ax_slider, 'Tanh Width', 0.0005, 0.01, valinit=0.001)
def update(num):
    global ln3    
    global ln4
    
    ln2.set_ydata((1/(200*x2)) * np.tanh(x2/twidth.val)**2)
    ln1.set_ydata(np.tanh(x1/twidth.val)**2)
    
    ln4.remove()
    maxX = 1.08866*twidth.val
    ln4, = ax.plot([maxX],
                   (1/200) * (0.5826/twidth.val), 'ko')
    
    ann.set_text("Width = %.2g"%(twidth.val))
    fig.canvas.draw_idle()

twidth.on_changed(update)

ax.set_yticks([])
ax.set_xticks([])
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)

ax.annotate(r"$\frac{1}{x}$", (0.001,5), (-0.03,5), arrowprops={'arrowstyle':"->",'connectionstyle':"arc3"}, fontsize=26)
ax.annotate(r"tanh($\frac{x}{w}$)$^2$", (0.075,1), (0.06,2), arrowprops={'arrowstyle':"->",'connectionstyle':"arc3"}, fontsize=22)
