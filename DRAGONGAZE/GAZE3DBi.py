'''
GAZE3DBi.py
DRAGONGAZE
Created by Bobbie Markwick on 01/09/2026.

Field-generic plotting script for 3D runs, plotting two variables.
Part of the DRAGONGAZE package -- run as a module, from the directory that
CONTAINS DRAGONGAZE/ (see package __init__.py for why):
    python -m DRAGONGAZE.GAZE3DBi <red_field> <blue_field> <axisp_pair>
Produces a mid-plane-slice bivariate heatmap, with red_field in red and blue_field in blue
See GAZE2DBi.py for the 2D equivalent.
'''

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize
import math

import h5py
from .HDF5_keys import *
import numpy as np

from .Config import *

from .FileUtils import *

import sys

def readMetadata(n):
    #Grid spacing and simulation time for frame n, read from the file's own HDF5 attributes
    with h5py.File(h5FileName(n), "r") as f:
        dx = f.attrs[key_dx]
        dy = f.attrs[key_dy]
        dz = f.attrs[key_dz]
        t = f.attrs[key_time]
    return  t,dx,dy,dz

def bivar(field1, field2, rng1, rng2):
    A = (np.asarray(field1, dtype=float) - rng1[0]) / (rng1[1] - rng1[0])
    B = (np.asarray(field2, dtype=float) - rng2[0]) / (rng2[1] - rng2[0])
    w00 = np.array([1, 1, 1])   # white
    w10 = np.array([1, 0, 0])   # red
    w01 = np.array([0, 0, 1])   # blue
    w11 = np.array([1, 0, 1])   # magenta
    return ((1 - A)[..., None] * (1 - B)[..., None] * w00 +
            A[..., None] * (1 - B)[..., None] * w10 +
            (1 - A)[..., None] * B[..., None] * w01 +
            A[..., None] * B[..., None] * w11)
            
def plotBivar(data, n, dx, dy, dim):
    #Renders a single mid-plane slice of one field's 3D array as a heatmap image and saves it to img_dir.
    #dim selects which slice: 0="yz" (fixed x, at the x mid-plane), 1="xz" (fixed y), 2="xy" (fixed z).
    nz, ny, nx, _ = data.shape
    #Calculate the axes (x_mode/y_mode/z_mode from Config.py choose where each origin is labelled -- see LORE Section 5.3)
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
        _data = data[:,:,nx//2]   #slice at the x mid-plane
    elif dim == 1:
        ax_w_px, ax_h_px = nx,nz
        extent = [xL,xR,zL,zR]
        axis_label = "xz"
        _data = data[:,ny//2,:]   #slice at the y mid-plane
    elif dim == 2:
        ax_w_px, ax_h_px = nx,ny
        extent = [xL,xR,yL,yR]
        axis_label = "xy"
        _data = data[nz//2,:,:]   #slice at the z mid-plane

   #Calculate image size: upscale small grids to a ~400px floor on each axis so a coarse simulation doesn't produce a postage-stamp image
    if ax_w_px < 400 or ax_h_px < 400:
        sx, sy = math.ceil(400/ax_w_px), math.ceil(400/ax_h_px)
        a = max([sx,sy,1])
        ax_w_px *= a
        ax_h_px *= a
    #Fixed pixel margins around the plot axes (converted to inches via dpi below), plus room on the right for a separate colorbar axis (see cax below)
    left_px   = 80
    right_px  = 80
    bottom_px = 70
    top_px    = 100
    fig_w_px = left_px + ax_w_px + right_px
    fig_h_px = bottom_px + ax_h_px + top_px

    fig = plt.figure(figsize=(fig_w_px / dpi, fig_h_px / dpi), dpi=dpi)
    ax = fig.add_axes([left_px / fig_w_px, bottom_px / fig_h_px, ax_w_px / fig_w_px, ax_h_px / fig_h_px])

    ax.imshow(_data, origin="lower", interpolation="bilinear")
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
    ax.set_title(f"[t = {t:.4f}{time_unit}]")
    fig.suptitle(plot_title, fontsize=22, y=0.96)
    #save
    png_name = imgFileName3D(n,plot_title, axis_label)
    fig.savefig(png_name, dpi=100)
    print(f"Wrote {png_name}")
    plt.close(fig)

#Compute each requested field's global [min, max] across every existing frame up front (FileUtils.field_range), so colour limits stay fixed for the whole run rather than jumping frame to frame.
rngs = {}
for field in sys.argv[1:]:
    rngs[field] = field_range(field)

#Plot every frame that currently exists (see FileUtils.fileExists); one image per field per frame
n = 0
while fileExists(n):
    t, dx, dy = readMetadata(n)
    f1, f2 = sys.argv[1], sys.argv[2]
    A,B = readField(n, f1), readField(n, f2)
    data = bivar(A,B, rngs[f1], rngs[f2])
    
    axis_code = sys.argv[3]
    if axis_code == "yz" or axis_code == "zy":   plotBivar(data,n,dx,dy,0)
    elif axis_code == "xz" or axis_code == "zx": plotBivar(data,n,dx,dy,1)
    elif axis_code == "xy" or axis_code == "yx": plotBivar(data,n,dx,dy,2)
    
    n+=1


