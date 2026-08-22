'''
GAZE2D.py
DRAGONGAZE
Created by Bobbie Markwick on 07/07/2026.

Field-generic plotting script for 2D runs.
Part of the DRAGONGAZE package -- run as a module, from the directory that
CONTAINS DRAGONGAZE/ (see package __init__.py for why):
    python -m DRAGONGAZE.GAZE2D <field1> <field2> ...
(e.g. `python -m DRAGONGAZE.GAZE2D rho p Bmag`).
Produces one heatmap PNG per named field per existing output frame.
Unlike GAZE1D.py/GAZE1DMHD.py this script has no fixed layout. See FileUtils.readField for the list of valid field names, and Config.py to control how each looks.
For 3D runs, see GAZE3D.py.
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
    t,dx,dy = 0,0,0
    with h5py.File(h5FileName(n), "r") as f:
        dx = f.attrs[key_dx]
        dy = f.attrs[key_dy]
        t = f.attrs[key_time]
    return  t,dx,dy

def plotField(data, n, dx, dy, rng, key):
    #Renders a single field's array as one heatmap image and saves it to img_dir.
    #rng is this field's precomputed [min, max] (see rngs below), shared across every frame so an animation's colour scale stays fixed rather than rescaling frame to frame.
    ny, nx = data.shape
    #Calculate image size: upscale small grids to a ~400px floor on each axis so a coarse simulation doesn't produce a postage-stamp image
    ax_w_px, ax_h_px = nx, ny
    if ax_w_px < 400 or ax_h_px < 400:
        sx, sy = math.ceil(400/nx), math.ceil(400/ny)
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

    #Calculate the axes (x_mode/y_mode from Config.py choose where the origin is labelled -- see LORE Section 5.3)
    if x_mode < 0: xR = 0; xL = nx*dx
    elif x_mode == 0: xR = nx*dx*0.5; xL = -xR
    else: xL = 0; xR = nx*dx
    if y_mode < 0: yR = 0; yL = ny*dy
    elif y_mode == 0: yR = ny*dy*0.5; yL = -yR
    else: yL = 0; yR = ny*dy
    extent = [xL,xR,yL,yR]
    #Scale logarithmically or linearly
    if log_plots[key]:
        arr = np.maximum(data, rng[0])
        norm=LogNorm(vmin = rng[0], vmax = rng[1])
    else:
        arr = data
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
    png_name = imgFileName2D(n,key)
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
    for field in sys.argv[1:]:
        data = readField(n, field)
        plotField(data,n,dx,dy,rngs[field],field)
    n+=1

