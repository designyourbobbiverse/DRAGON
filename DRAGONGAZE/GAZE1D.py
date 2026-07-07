'''
HDF5_keys.py
DRAGONGAZE
Created by Bobbie Markwick on 07/07/2026.
'''

import matplotlib.pyplot as plt
import h5py
import numpy
from HDF5_keys import *
from Config import *
import os

def h5FileName(n):
    return hdf_dir + "/" + h5_base_filename + "_" + str(n).zfill(5) + ".h5"
    
def imgFileName(n):
    return img_dir + "/" + img_base_filename + "_" + str(n).zfill(5) + ".png"

def fileExists(n):
    return os.path.exists(h5FileName(n))
    
def readFile(n):
    rho, vx, p, E, dx, t = [],[],[],[], 0,0
    with h5py.File(h5FileName(n), "r") as f:
        rho = f[key_rho][:]              # read dataset into NumPy array
        vx  = f[key_vx][:]
        p   = f[key_p][:]
        E   = f[key_E][:]

        dx = f.attrs[key_dx]
        t = f.attrs[key_time]
    return rho, vx, p, E, t, dx
    
def plotFile(n, rho,vx,p,E,t,dx):
    x = numpy.arange(len(rho)) * dx
    
    fig, axs = plt.subplots(2,2, figsize=(10, 9), sharex=True)
    axs = axs.ravel()
    
    if log_plots: axs[0].semilogy(x, rho)
    else: axs[0].plot(x, rho)
    axs[0].set_ylabel(rho_label)
    axs[0].set_title(rho_title)

    axs[1].plot(x, vx)
    axs[1].set_ylabel(v_label)
    axs[1].set_title(v_title)

    if log_plots: axs[2].semilogy(x, p)
    else: axs[2].plot(x, p)
    axs[2].set_ylabel(p_label)
    axs[2].set_xlabel(x_label)
    axs[2].set_title(p_title)


    if log_plots: axs[3].semilogy(x, E)
    else: axs[3].plot(x, E)
    axs[3].set_ylabel(E_label)
    axs[3].set_xlabel(x_label)
    axs[3].set_title(E_title)

    title = f"{plot_title} [t = {t:.4f}{time_unit}]"
    fig.suptitle(title)
    fig.tight_layout()
    
    png_name = imgFileName(n)
    fig.savefig(png_name, dpi=64)
    print(f"Wrote {png_name}")
    plt.close(fig)


n = 0
while fileExists(n):
    rho, vx, p, E, t, dx = readFile(n)
    plotFile(n,rho,vx,p,E,t,dx)
    n+=1


