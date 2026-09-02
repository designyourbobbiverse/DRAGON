'''
FileUtils.py
DRAGONGAZE
Created by Bobbie Markwick on 07/07/2026.

Shared low-level helpers imported (via `from .FileUtils import *`) across GAZE*.py plotting scripts:
building input/output file paths, checking which frames exist, reading fields out of a frame's HDF5 file,
and computing padded value ranges for consistent axis/colour limits across an animation.

Part of the DRAGONGAZE package (see __init__.py) -- imports here are package-relative,
so this module (and every GAZE*.py script) must be run as part of the package rather
than as a bare script; see __init__.py for the invocation this requires.
'''

from .Config import *
import os
import h5py
import numpy as np
from .HDF5_keys import *

'''
File Utilities
'''
def h5FileName(n):
    #n<0 is an internal sentinel used to get the base path without a frame number/extension (see the recursive call below)
    if n < 0: return hdf_dir + "/" + h5_base_filename
    return h5FileName(-1) + "_" + str(n).zfill(5) + ".h5"

def imgFileName1D(n):
    #1D output isn't per-field (GAZE1D.py/GAZE1DMHD.py always plot every field in one multi-panel figure), so no prefix is needed
    return img_dir + "/" + img_base_filename + "_" + str(n).zfill(5) + ".png"

def imgFileName2D(n, key):
    prefix = fileprefix[key] if key in fileprefix else key
    return img_dir + "/" + prefix + "_" + img_base_filename + "_" + str(n).zfill(5) + ".png"

def imgFileName3D(n, key, axis):
    #axis is the two-letter slice label ("xy"/"xz"/"yz") appended by GAZE3D.py, not necessarily matching command-line order
    prefix = fileprefix[key] if key in fileprefix else key
    return img_dir + "/" + prefix + "_" + axis + "_" + img_base_filename + "_" + str(n).zfill(5) + ".png"

def fileExists(n):
    #Used by every script to discover how many frames exist: loop n = 0, 1, 2, ... until this returns False
    return os.path.exists(h5FileName(n))

'''
Reading Arrays
'''
def readField(n, field):
    #Reads one field's array for frame n.
    #"Bmag" isn't stored in the HDF5 file directly (see HDF5_keys.py), so it's computed here from the stored Bx/By/Bz components instead of looked up like everything else.
    #To add another derived (rather than stored) field, add another special case here following this pattern.
    with h5py.File(h5FileName(n), "r") as f:
        if field == "Bmag":
            Bx = f[keys["Bx"]][:]
            By = f[keys["By"]][:]
            Bz = f[keys["Bz"]][:]
            return np.sqrt(Bx**2 + By**2 + Bz**2)
        if field in keys:  return f[keys[field]][:]
        return f[key_passive + field][:]
    return  []


'''
Ranges
'''
def pad_range(ranges, frac=0.05):
    #Pads a [lo, hi] pair (in place) by a fractional margin so plotted values don't sit flush against the axis/colour limits
    lo, hi = ranges[0], ranges[1]
    width = hi - lo
    # If range is flat or nearly flat, invent a small symmetric margin
    if width == 0:
        if lo != 0: margin = abs(lo) * frac
        else: margin = frac
    else: margin = width * frac

    ranges[0] -= margin
    ranges[1] += margin

def field_range(key):
    #Scans every existing frame to find a field's global [min, max], so an animation's axis/colour limits can be fixed once up front (via pad_range) instead of jumping around from frame to frame.
    rng = [np.inf, -np.inf]

    n = 0
    while fileExists(n):
        val = readField(n,key)
        n += 1

        rng[0] = min(rng[0], np.min(val))
        rng[1] = max(rng[1], np.max(val))

    return rng
