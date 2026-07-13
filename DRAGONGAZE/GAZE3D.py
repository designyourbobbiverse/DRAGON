'''
GAZE1D.py
DRAGONGAZE
Created by Bobbie Markwick on 07/07/2026.
'''

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize
import math

import h5py
from HDF5_keys import *
import numpy as np

from Config import *

from FileUtils import *
from Ranges import range

import sys

def readMetadata(n):
    with h5py.File(h5FileName(n), "r") as f:
        dx = f.attrs[key_dx]
        dy = f.attrs[key_dy]
        dz = f.attrs[key_dz]
        t = f.attrs[key_time]
    return  t,dx,dy,dz
    
def readFile(n, field):
    with h5py.File(h5FileName(n), "r") as f:
        return f[keys[field]][:]
    return  []
    
def plotField(data, n, dx, dy, dz, rng, key, dim):
    nz, ny, nx = data.shape
    #Calculate the axes
    if x_mode < 0: xR = 0; xL = nx*dx
    elif x_mode == 0: xR = nx*dx*0.5; xL = -xR
    else: xL = 0; xR = nx*dx
    
    if y_mode < 0: yR = 0; yL = ny*dy
    elif y_mode == 0: yR = ny*dy*0.5; yL = -yR
    else: yL = 0; yR = ny*dy
    
    if z_mode < 0: zR = 0; zL = nz*dz
    elif z_mode == 0: zR = nz*dz*0.5; zL = -zR
    else: zL = 0; zR = nz*dz
    
    if dim == 0:
        ax_w_px, ax_h_px = ny,nz
        extent = [yL,yR,zL,zR]
        axis_label = "yz"
        _data = data[:,:,nx//2]
    elif dim == 1:
        ax_w_px, ax_h_px = nx,nz
        extent = [xL,xR,zL,zR]
        axis_label = "xz"
        _data = data[:,ny//2,:]
    elif dim == 2:
        ax_w_px, ax_h_px = nx,ny
        extent = [xL,xR,yL,yR]
        axis_label = "xy"
        _data = data[nz//2,:,:]
   
   #Calculate image size
    if ax_w_px < 400 or ax_h_px < 400:
        sx, sy = math.ceil(400/ax_w_px), math.ceil(400/ax_h_px)
        a = max([sx,sy,1])
        ax_w_px *= a
        ax_h_px *= a
    left_px   = 80
    right_px  = 100
    bottom_px = 70
    top_px    = 90
    fig_w_px = left_px + ax_w_px + right_px
    fig_h_px = bottom_px + ax_h_px + top_px

    fig = plt.figure(figsize=(fig_w_px / dpi, fig_h_px / dpi), dpi=dpi)
    ax = fig.add_axes([left_px / fig_w_px, bottom_px / fig_h_px, ax_w_px / fig_w_px, ax_h_px / fig_h_px])
    
    #Scale logarithmically or linearly
    if log_plots[key]:
        arr = np.maximum(_data, rng[0])
        norm=LogNorm(vmin = rng[0], vmax = rng[1])
    else:
        arr = _data
        norm=Normalize(vmin = rng[0], vmax = rng[1])
    #plot the data
    img = ax.imshow(
        arr,
        origin="lower",
        interpolation="nearest",
        cmap = cmaps[key],
        extent = extent,
        norm=norm
    )
    #labels and titles
    if dim == 0:
        ax.set_xlabel(y_label)
        ax.set_ylabel(z_label)
    elif dim == 1:
        ax.set_xlabel(x_label)
        ax.set_ylabel(z_label)
    elif dim == 2:
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
    ax.set_title(f"{titles[key]} [t = {t:.4f}{time_unit}]")
    fig.suptitle(plot_title, fontsize=22, y=0.96)
    #Colour Bar
    cbar_left_px = left_px + ax_w_px + 20
    cbar_w_px = 25
    cax = fig.add_axes([
        cbar_left_px / fig_w_px,
        bottom_px / fig_h_px,
        cbar_w_px / fig_w_px,
        ax_h_px / fig_h_px,
    ])
    fig.colorbar(img, cax=cax)
    #save
    png_name = imgFileName3D(n,key, axis_label)
    fig.savefig(png_name, dpi=100)
    print(f"Wrote {png_name}")
    plt.close(fig)


rngs = {}
for field in sys.argv[1:]:
    n = field.find("-")
    if n >= 0 and field[:n] not in rngs: rngs[field[:n]] = range(keys[field[:n]])
    else: rngs[field] = range(keys[field])
    
    
    
n = 0
while fileExists(n):
    t, dx, dy, dz = readMetadata(n)
    for field in sys.argv[1:]:
        idx = field.find("-")
        if idx >= 0:
            axis_code = field[idx+1:]
            field = field[:idx]
            data = readFile(n, field)
            if axis_code == "yz" or axis_code == "zy":   plotField(data,n,dx,dy,dz,rngs[field],field,0)
            elif axis_code == "xz" or axis_code == "zx": plotField(data,n,dx,dy,dz,rngs[field],field,1)
            elif axis_code == "xy" or axis_code == "yx": plotField(data,n,dx,dy,dz,rngs[field],field,2)
        else:
            data = readFile(n, field)
            plotField(data,n,dx,dy,dz,rngs[field],field,0)
            plotField(data,n,dx,dy,dz,rngs[field],field,1)
            plotField(data,n,dx,dy,dz,rngs[field],field,2)
    n+=1

