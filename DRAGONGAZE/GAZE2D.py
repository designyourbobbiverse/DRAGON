'''
GAZE1D.py
DRAGONGAZE
Created by Bobbie Markwick on 07/07/2026.
'''

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize

import h5py
from HDF5_keys import *
import numpy as np

from Config import *

from FileUtils import *
from Ranges import range

import sys

def readMetadata(n):
    t,dx,dy = 0,0,0
    with h5py.File(h5FileName(n), "r") as f:
        dx = f.attrs[key_dx]
        dy = f.attrs[key_dy]
        t = f.attrs[key_time]
    return  t,dx,dy
    
def readFile(n, field):
    with h5py.File(h5FileName(n), "r") as f:
        return f[keys[field]][:]
    return  []
    
def plotField(data, n, dx, dy, rng, key):
    ny, nx = data.shape

    x = np.arange(nx) * dx
    y = np.arange(ny) * dy
    
    if x_mode < 0: x.flip()
    elif x_mode == 0: x -= max(x)/2
    
    if y_mode < 0: y.flip()
    elif y_mode == 0: y -= max(y)/2
    
    if log_plots[key]:
        arr = np.maximum(data, rng[0])
        norm=LogNorm(vmin = rng[0], vmax = rng[1])
    else:
        arr = data
        norm=Normalize(vmin = rng[0], vmax = rng[1])
    
    plt.imshow(
        arr,
        origin="lower",
        interpolation="nearest",
        aspect="equal",
        cmap = cmaps[key],
        norm=norm
    )
    plt.colorbar(label=labels[key])
    png_name = imgFileName(n,key)
    plt.savefig(png_name, dpi=256)
    print(f"Wrote {png_name}")
    plt.close()
    #png_name = imgFileName(n)
    ##plt.savefig(png_name, dpi=64)
    #print(f"Wrote {png_name}")
    #plt.close(fig)


rngs = {}
for field in sys.argv[1:]:
    rngs[field] = range(keys[field])
    
n = 0
while fileExists(n):
    t, dx, dy = readMetadata(n)
    for field in sys.argv[1:]:
        data = readFile(n, field)
        plotField(data,n,dx,dy,rngs[field],field)
    n+=1

