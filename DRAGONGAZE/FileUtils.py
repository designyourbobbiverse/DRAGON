'''
GAZE1D.py
DRAGONGAZE
Created by Bobbie Markwick on 07/07/2026.
'''

from Config import *
import os
import h5py
import numpy as np
from HDF5_keys import *

'''
File Utilities
'''
def h5FileName(n):
    if n < 0: return hdf_dir + "/" + h5_base_filename
    return h5FileName(-1) + "_" + str(n).zfill(5) + ".h5"
    
def imgFileName1D(n):
    return img_dir + "/" + img_base_filename + "_" + str(n).zfill(5) + ".png"

def imgFileName2D(n, key):
    return img_dir + "/" + fileprefix[key] + "_" + img_base_filename + "_" + str(n).zfill(5) + ".png"

def imgFileName3D(n, key, axis):
    return img_dir + "/" + fileprefix[key] + "_" + axis + "_" + img_base_filename + "_" + str(n).zfill(5) + ".png"

def fileExists(n):
    return os.path.exists(h5FileName(n))

'''
Reading Arrays
'''
def readField(n, field):
    with h5py.File(h5FileName(n), "r") as f:
        if field == "Bmag":
            Bx = f[keys["Bx"]][:]
            By = f[keys["By"]][:]
            Bz = f[keys["Bz"]][:]
            return np.sqrt(Bx**2 + By**2 + Bz**2)
        return f[keys[field]][:]
    return  []


'''
Ranges
'''
def pad_range(ranges, frac=0.05):
    lo, hi = ranges[0], ranges[1]
    width = hi - lo
    # If range is flat or nearly flat, invent a small symmetric margin
    if width == 0:
        if lo != 0: margin = abs(lo) * frac
        else: margin = frac
    else: margin = width * frac

    ranges[0] -= margin
    ranges[1] += margin
    
def range(key):
    rng = [np.inf, -np.inf]

    n = 0
    while fileExists(n):
        val = readField(n,key)
        n += 1

        rng[0] = min(rng[0], np.min(val))
        rng[1] = max(rng[1], np.max(val))

    return rng
