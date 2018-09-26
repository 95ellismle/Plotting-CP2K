#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  7 13:01:57 2018

@author: mellis
"""

import os

# Adds the slash to the end of a folder string
def folder_correct(folder):
    if os.path.isdir(folder):
        if folder[-1] != '/':
            folder = folder + '/'
        return folder
    else:
        return None

# Makes a relative folder path absolute
def make_fold_abs(folder):
    if os.path.isdir(folder):
        return folder_correct(os.path.abspath(folder))
    else:
        raise SystemExit("Sorry I can't find the folder:\n\n\t* %s"%folder)