#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May  1 12:49:10 2019

@author: mellis
"""

import numpy as np
import matplotlib.pyplot as plt

dx = 0.005
x = np.arange(-10, 10, dx)
tullyModel = 3

def create_H1(x, A=0.01, B=1.6, C=0.005, D=1.0):
    """
    Will create the Hamiltonian in Tully Model 1
    """
    if x > 0:
        V11 = A*(1 - np.exp(-B*x))
    elif x < 0:
        V11 = -A*(1 - np.exp(B*x))
    else:
        V11 = 0

    V22 = -V11

    V12 = C * np.exp(-D*(x**2))

    return np.matrix([[V11, V12], [V12, V22]])


def create_H2(x, A=0.1, B=0.28, C=0.015, D=0.06, E0=0.05):
    """
    Will create the Hamiltonian in Tully Model 2
    """
    V11 = 0
    V22 = -A * np.exp(-B*(x**2)) + E0
    V12 = C*np.exp(-D*(x**2))

    return np.matrix([[V11, V12], [V12, V22]])


def create_H3(x, A=6e-4, B=0.1, C=0.9):
    """
    Will create the Hamiltonian in Tully Model 3
    """
    V11 = A
    if x < 0:
        V12 = B*np.exp(C*x)
    elif x > 0:
        V12 = B*(2-np.exp(-C*x))
    else:
        V12 = B
    V22 = -A

    return np.matrix([[V11, V12], [V12, V22]])


if tullyModel == 1:
    allH = [create_H1(i) for i in x]
if tullyModel == 2:
    allH = [create_H2(i) for i in x]
if tullyModel == 3:
    allH = [create_H3(i) for i in x]

def getEigProps(H):
    """
    Will return eigen properties (values and vectors) that are usable
    (corrected) minus signs in the code.
    """
    E, U = np.linalg.eigh(H)
    if tullyModel == 2:
        E1, _ = np.linalg.eig(H)
        if E1[0] > E1[1]:
            U[0, 1] = -U[0, 1]
            U[1, 1] = -U[1, 1]
    return E, U

adData = [getEigProps(H) for H in allH]
adEner = np.array([i[0] for i in adData])
adStates = np.array([i[1] for i in adData])


# Calculate the NACV
gradStates = np.array(np.gradient(adStates, dx, axis=0))
NACV = [np.matmul(phi, grad) for phi, grad in zip(adStates, gradStates)]
NACV = np.array(NACV)

plt.plot(x, adEner[:, 0], 'g-', label=r"E$^{ad}_{1}$")
plt.plot(x, adEner[:, 1], 'b-', label=r"E$^{ad}_{2}$")

if tullyModel == 1:
    NACVx = np.concatenate((x[x > 0+dx], x[x < 0-dx]))
    NACV = np.concatenate((NACV[x > 0+dx], NACV[x < 0-dx]))
    plt.plot(NACVx, NACV[:, 1, 0]/50, 'r.', label=r"$\frac{\mathbf{d}_{10}}{50}$")
    plt.title("Tully Model 1 -Single Avoided Crossing")
elif tullyModel == 2:
    plt.plot(x, NACV[:, 1, 0]/12, 'r.', label=r"$\frac{\mathbf{d}_{10}}{12}$")
    plt.title("Tully Model 2 -Dual Avoided Crossing")
elif tullyModel == 3:
    plt.plot(x, NACV[:, 0, 1], 'r.', label=r"$\mathbf{d}_{10}$")
    plt.title("Tully Model 3 -Extended Coupling")

plt.xlabel(r"$\mathbf{R}$")
#plt.xlabel(r"E")
plt.legend()
