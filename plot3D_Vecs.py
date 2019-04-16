#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 14:05:06 2019

@author: mellis
"""
from load import LOAD as ld

import mayavi.mlab as mlab
import numpy as np
import os


folder = "/home/oem/Data/CTMQC/QM_Vis"
sizes = {'C': 1.2, 'H': 0.6, 'Ne': 0.5}
plotting_parameters = ["qm_t"]
replicas = 'all'
plot = True
min_time = 0
max_time = 'all'

saveOn = True
saveStride = 20


imgFolder = "./img"
draw = {i: i in plotting_parameters for i in ld.dependencies}
pts = {}
plotting_parameters.append("pos")

CTMQCdata = ld.LoadData(folder=folder,
                        reps=replicas,
                        plot_params=plotting_parameters,
                        avg_on=True,
                        minTime=min_time,
                        maxTime=max_time)

mlab.figure(bgcolor=(0.2, 0.2, 0.2), size=(1000, 1000))
nsteps = CTMQCdata.num_pos_steps


def getActiveAtoms(pos, cols, atomName=False):
    """
    Will return an array of positions that are active depending on the cols
    variable. This will obviously have to change with the AOM_COEFF file though
    the cols are just for proof of concept.
    """
    if atomName:
        mask = cols[0] == atomName
    else:
        mask = cols[0] != 'Ne'
    return np.array([i[mask] for i in pos]), np.array([i[mask] for i in cols])


def reshapeData(arrData):
   """
   Will flatten the 2 columns that are to be plotted (i.e. num atoms, num states etc...)
   Will keep the num steps and num dim columns the same.
   Will return an array of shape (Nstep, 3, Natom*Nrep)
   """
   arrData = np.swapaxes(arrData, 2, 3)
   arrData = np.swapaxes(arrData, 1, 3)
   shape = (arrData.shape[0]*arrData.shape[1], arrData.shape[2], arrData.shape[3])
   arrData = np.reshape(arrData, shape)
   return arrData


if draw['pos']: 
   Cpos = [getActiveAtoms(CTMQCdata.all_pos_data[i][0],
                            CTMQCdata.all_pos_data[i][1], "C")[0]
             for i in CTMQCdata.all_pos_data]
   Cpos = reshapeData(Cpos)

   Hpos = np.array([getActiveAtoms(CTMQCdata.all_pos_data[i][0],
                                     CTMQCdata.all_pos_data[i][1], "H")[0]
                      for i in CTMQCdata.all_pos_data])
   Hpos = reshapeData(Hpos)

# Get vital information for plotting any vectors
if any((i != 'pos' for i in plotting_parameters)):
   tmp = [getActiveAtoms(CTMQCdata.all_pos_data[i][0], 
                         CTMQCdata.all_pos_data[i][1]) for i in CTMQCdata.all_pos_data]
   allActPos = np.array([i[0] for i in tmp])
   allActPos = reshapeData(allActPos)

# Get the qm info
if draw['qm_t']:
   qm = np.array([CTMQCdata.all_QM_0_data[i][0] for i in CTMQCdata.all_QM_0_data])
   qm = reshapeData(qm)
     

# Plot atoms
if draw['pos']:
   pts['Cpos'] = mlab.points3d(Cpos[:, 0, 0],
                               Cpos[:, 1, 0],
                               Cpos[:, 2, 0],
                               scale_factor=sizes['C'], color=(0, 0, 0))
   pts['Hpos'] = mlab.points3d(Hpos[:, 0, 0],
                               Hpos[:, 1, 0],
                               Hpos[:, 2, 0],
                               scale_factor=sizes['H'], color=(1, 1, 0))

if draw['qm_t']:
   pts['qm'] = mlab.quiver3d(allActPos[:, 0, 0],
                         allActPos[:, 1, 0],
                         allActPos[:, 2, 0], 
                         qm[:, 0, 0],
                         qm[:, 1, 0], 
                         qm[:, 2, 0])


if not os.path.isdir(imgFolder):
   os.makedirs(imgFolder)
numZeros = len(str(nsteps - 1))
filePaths = [(numZeros - len(str(i)))*"0" + str(i) for i in range(nsteps)]
filePaths = ["./img/%s.jpg" % i for i in filePaths]
#mlab.savefig(filePaths[0], (1000, 1000), mlab.gcf())

@mlab.animate(delay=10)
def anim(draw, pts):
    for step in range(1, nsteps):
        # Update Positions
        if draw['pos']:
           pts['Cpos'].mlab_source.set(x=Cpos[:, 0, step],
                                       y=Cpos[:, 1, step],
                                       z=Cpos[:, 2, step])
           pts['Hpos'].mlab_source.set(x=Hpos[:, 0, step],
                                       y=Hpos[:, 1, step],
                                       z=Hpos[:, 2, step])
        
        if draw['qm_t']:
           pts['qm'].mlab_source.set(x=allActPos[:, 0, step],
                                     y=allActPos[:, 1, step],
                                     z=allActPos[:, 2, step],
                                     u=qm[:, 0, step],
                                     v=qm[:, 1, step],
                                     w=qm[:, 2, step])

        if saveOn and step % saveStride:
           mlab.savefig(filePaths[step], figure=mlab.gcf())
        yield

anim(draw, pts)
mlab.show()
