#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 14:05:06 2019

@author: mellis
"""

from load import load_pos
from load import load_QM
from load import load_tintf
from load import load_frc

import mayavi.mlab as mlab
import numpy as np


folder = "/home/oem/Data/CTMQC/CTMQCForceEhrenCoeff/"
sizes = {'C': 1.8, 'H': 0.7, 'Ne': 1}
replica = 15

allPosData = load_pos.load_all_pos_in_folder(folder, reps=[replica])
allQMData = load_QM.load_all_QM_0_in_folder(folder, reps=[replica])
allFrcData = load_frc.load_all_frc_in_folder(folder, reps=[replica])
allAMData = load_tintf.load_all_tintf_in_folder(folder, reps=[replica])

pKeys = list(allPosData.keys())
qmKeys = list(allQMData.keys())
frcKeys = list(allFrcData.keys())
amKeys = list(allAMData.keys())


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

pos, rcols, timesteps = allPosData[pKeys[0]]
qm, qmcols, timesteps = allQMData[qmKeys[0]]
frc, fcols, timesteps = allFrcData[frcKeys[0]]
am, amcols, timesteps = allAMData[amKeys[0]]

nsteps = len(qm)
if (nsteps != len(pos) or nsteps != len(frc)):
    print("Len(pos) = %i" % len(pos))
    print("Len(frc) = %i" % len(frc))
    print("Len(qm) = %i" % len(qm))
    raise SystemExit("Incorrect number of steps")

# Get vital information about the quantum momentum
qmPos, tmpCols = getActiveAtoms(pos, rcols)
qmx, qmy, qmz = qm[0].T
posx, posy, posz = qmPos[0].T
# Separate positions into carbon and hydrogen
carbonPos, _ = getActiveAtoms(qmPos, tmpCols, "C")
hydrogPos, _ = getActiveAtoms(qmPos, tmpCols, "H")
neonPos, _ = getActiveAtoms(pos, rcols, "Ne")
# Init forces
activeFrc, _ = getActiveAtoms(frc, fcols)
frcx, frcy, frcz = activeFrc[0].T
#Init ad Momentum
stateMom, _ = load_tintf.find_in_histF(am, amcols, {"state": 0,
                                                    "step_num": 0})
amx, amy, amz = stateMom.T


# Plot atoms
cPts = mlab.points3d(*carbonPos[0].T, scale_factor=sizes['C'], color=(0, 0, 0))
hPts = mlab.points3d(*hydrogPos[0].T, scale_factor=sizes['H'], color=(1, 1, 0))
nePts = mlab.points3d(*hydrogPos[0].T, scale_factor=sizes['Ne'], color=(1, 1, 1))
# Plot vectors
#qmPts = mlab.quiver3d(posx, posy, posz, qmx, qmy, qmz)
frcPts = mlab.quiver3d(posx, posy, posz, frcx, frcy, frcz)
#amPts = mlab.quiver3d(posx, posy, posz, amx, amy, amz)


@mlab.animate(delay=10)
def anim():
    for i in range(nsteps):
        cPts.mlab_source.points = carbonPos[i]
        hPts.mlab_source.points = hydrogPos[i]
        nePts.mlab_source.points = neonPos[i]

        
        stateMom, _ = load_tintf.find_in_histF(am,
                                               amcols,
                                               {"state": 1,
                                                "step_num": i})
        #amPts.mlab_source.points = qmPos[i]
        #amPts.mlab_source.u = stateMom.T[0]
        #amPts.mlab_source.v = stateMom.T[1]
        #amPts.mlab_source.w = stateMom.T[2]
        
        frcPts.mlab_source.points = qmPos[i]
        frcPts.mlab_source.u = activeFrc[i].T[0]
        frcPts.mlab_source.v = activeFrc[i].T[1]
        frcPts.mlab_source.w = activeFrc[i].T[2]

        #qmPts.mlab_source.points = qmPos[i]
        #qmPts.mlab_source.u = qm[i].T[0]
        #qmPts.mlab_source.v = qm[i].T[1]
        #qmPts.mlab_source.w = qm[i].T[2]
        yield

anim()
mlab.show()
