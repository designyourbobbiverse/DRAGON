'''
GAZE1D.py
DRAGONGAZE
Created by Bobbie Markwick on 07/07/2026.
'''

import matplotlib.pyplot as plt
import h5py
import numpy as np
from HDF5_keys import *
from Config import *
from FileUtils import *
from Ranges import range


    
def readFile(n):
    rho, vx, vy,By,p, E, dx, t = [],[],[],[],[],[], 0,0
    with h5py.File(h5FileName(n), "r") as f:
        rho = f[key_rho][:]              # read dataset into NumPy array
        vx  = f[key_vx][:]
        vy  = f[key_vy][:]
        By  = f[key_By][:]
        p   = f[key_p][:]
        E   = f[key_E][:]

        dx = f.attrs[key_dx]
        t = f.attrs[key_time]
    return rho, vx, vy, By, p, E, t, dx
    
def plotFile(n, rho,vx,vy,By,p,E,t,dx, rho_rng, v_rng, vy_rng, B_rng, p_rng, E_rng):
    x = np.arange(len(rho)) * dx
    
    if x_mode < 0: x.flip()
    elif x_mode == 0: x -= max(x)/2
    
    fig, axs = plt.subplots(3,2, figsize=(10, 14), sharex=True)
    axs = axs.ravel()
    
    if log_plots: axs[0].semilogy(x, rho)
    else: axs[0].plot(x, rho)
    axs[0].set_ylim(rho_rng[0], rho_rng[1])
    axs[0].set_ylabel(rho_label)
    axs[0].set_title(rho_title)

    axs[1].plot(x, vx)
    axs[1].set_ylim(v_rng[0], v_rng[1])
    axs[1].set_ylabel(v_label)
    axs[1].set_title(v_title)

    axs[2].plot(x, vy)
    axs[2].set_ylim(vy_rng[0], vy_rng[1])
    axs[2].set_ylabel(vy_label)
    axs[2].set_title(vy_title)

    if log_plots: axs[3].semilogy(x, p)
    else: axs[3].plot(x, p)
    axs[3].set_ylim(p_rng[0], p_rng[1])
    axs[3].set_ylabel(p_label)
    axs[3].set_title(p_title)


    axs[4].plot(x, By)
    axs[4].set_ylim(B_rng[0], B_rng[1])
    axs[4].set_ylabel(B_label)
    axs[4].set_xlabel(x_label)
    axs[4].set_title(B_title)
    
    if log_plots: axs[5].semilogy(x, E)
    else: axs[5].plot(x, E)
    axs[5].set_ylim(E_rng[0], E_rng[1])
    axs[5].set_ylabel(E_label)
    axs[5].set_xlabel(x_label)
    axs[5].set_title(E_title)
    
   


    title = f"{plot_title} [t = {t:.4f}{time_unit}]"
    fig.suptitle(title)
    fig.tight_layout()
    
    png_name = imgFileName(n)
    fig.savefig(png_name, dpi=64)
    print(f"Wrote {png_name}")
    plt.close(fig)


rho_rng = range(key_rho)
v_rng = range(key_vx)
vy_rng = range(key_vy)
B_rng = range(key_By)
p_rng = range(key_p)
E_rng = range(key_E)

n = 0
while fileExists(n):
    rho, vx,vy,By, p, E, t, dx = readFile(n)
    plotFile(n,rho,vx,vy,By,p,E,t,dx, rho_rng, v_rng, vy_rng, B_rng, p_rng, E_rng)
    n+=1


