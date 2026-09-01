'''
GAZE3D.py
DRAGONGAZE
Created by Bobbie Markwick on 07/07/2026.

Field-generic plotting script for 3D runs.
Part of the DRAGONGAZE package -- run as a module, from the directory that
CONTAINS DRAGONGAZE/ (see package __init__.py for why):
    python -m DRAGONGAZE.GAZE3D <field1>[-<axis pair>] <field2>[-<axis pair>] ...
(e.g. `python -m DRAGONGAZE.GAZE3D rho Bmag-xy`);
produces a mid-plane-slice heatmap PNG per named field (and per requested slice) per existing output frame.
A field with no "-xy"/"-xz"/"-yz" suffix (order within the pair doesn't matter, e.g. "yx" == "xy") produces all three mid-plane slices.
See GAZE2D.py for the 2D equivalent.
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

def plotField(data, n, dx, dy, dz, rng, key, dim):
    #Renders a single mid-plane slice of one field's 3D array as a heatmap image and saves it to img_dir.
    #dim selects which slice: 0="yz" (fixed x, at the x mid-plane), 1="xz" (fixed y), 2="xy" (fixed z).
    #rng is this field's precomputed [min, max] (see rngs below), shared across every frame so an animation's colour scale stays fixed rather than rescaling frame to frame.
    nz, ny, nx = data.shape
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
    right_px  = 100
    bottom_px = 70
    top_px    = 90
    fig_w_px = left_px + ax_w_px + right_px
    fig_h_px = bottom_px + ax_h_px + top_px

    fig = plt.figure(figsize=(fig_w_px / dpi, fig_h_px / dpi), dpi=dpi)
    ax = fig.add_axes([left_px / fig_w_px, bottom_px / fig_h_px, ax_w_px / fig_w_px, ax_h_px / fig_h_px])

    #Scale logarithmically or linearly
    if key in log_plots and log_plots[key]:
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
        cmap = cmaps[key] if key in cmaps else cmaps["default"],
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
    if key in titles: ax.set_title(f"{titles[key]} [t = {t:.4f}{time_unit}]")
    else: ax.set_title(f"{key} [t = {t:.4f}{time_unit}]")
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


#Compute each requested field's global [min, max] across every existing frame up front (FileUtils.field_range), so colour limits stay fixed for the whole run rather than jumping frame to frame.
#Keyed by the bare field name (the "-xy" axis suffix, if any, is stripped first) since the range doesn't depend on which slice is being plotted.
#rngs[field[:n]] is only computed once even if the same field is requested twice with two different suffixes (e.g. "rho-xy rho-xz").
rngs = {}
for field in sys.argv[1:]:
    n = field.find("-")
    if n >= 0 and field[:n] not in rngs: rngs[field[:n]] = field_range(field[:n])
    else: rngs[field] = field_range(field)



#Plot every frame that currently exists (see FileUtils.fileExists)
n = 0
while fileExists(n):
    t, dx, dy, dz = readMetadata(n)
    for field in sys.argv[1:]:
        #Parse an optional "-<axis pair>" suffix off this command-line argument.
        #Order within the pair doesn't matter (xy/yx, xz/zx, yz/zy are all accepted) since they select the same slice.
        idx = field.find("-")
        if idx >= 0:
            axis_code = field[idx+1:]
            field = field[:idx]
            data = readField(n, field)
            if axis_code == "yz" or axis_code == "zy":   plotField(data,n,dx,dy,dz,rngs[field],field,0)
            elif axis_code == "xz" or axis_code == "zx": plotField(data,n,dx,dy,dz,rngs[field],field,1)
            elif axis_code == "xy" or axis_code == "yx": plotField(data,n,dx,dy,dz,rngs[field],field,2)
        else:
            #No suffix: render all three mid-plane slices for this field
            data = readField(n, field)
            plotField(data,n,dx,dy,dz,rngs[field],field,0)
            plotField(data,n,dx,dy,dz,rngs[field],field,1)
            plotField(data,n,dx,dy,dz,rngs[field],field,2)
    n+=1

