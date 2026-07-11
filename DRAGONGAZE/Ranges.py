'''
Ranges.py
DRAGONGAZE
Created by Bobbie Markwick on 07/07/2026.
'''
import numpy as np
import h5py
from FileUtils import *

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
        val = []
        with h5py.File(h5FileName(n), "r") as f:
            val = f[key][:] # read dataset into NumPy array
        n += 1

        rng[0] = min(rng[0], np.min(val))
        rng[1] = max(rng[1], np.max(val))

    return rng
