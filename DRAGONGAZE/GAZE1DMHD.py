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


    
def readFile(n):
    rho, vx, vy,By,p, E, dx, t = [],[],[],[],[],[], 0,0
    with h5py.File(h5FileName(n), "r") as f:
        rho = f[key_rho][:]              # read dataset into NumPy array]
        p   = f[key_p][:]
        E   = f[key_E][:]

        vx  = f[key_vx][:]
        vy  = f[key_vy][:]
        vz  = f[key_vz][:]

        By  = f[key_By][:]
        Bz  = f[key_Bz][:]

        dx = f.attrs[key_dx]
        t = f.attrs[key_time]
    return  rho, vx,vy,vz,By,Bz,p, E, t, dx
    
def plotFile(n, rho,vx,vy,By,p,E,t,dx, rho_rng, v_rng, p_rng, E_rng, vy_rng, By_rng, vz_rng, Bz_rng):
    x = np.arange(len(rho)) * dx
    
    if x_mode < 0: np.flip(x)
    elif x_mode == 0: x -= max(x)/2
    
    fig, axs = plt.subplots(4,2, figsize=(10, 20), sharex=True)
    axs = axs.ravel()
    
    if log_plots["rho"]: axs[0].semilogy(x, rho)
    else: axs[0].plot(x, rho)
    axs[0].set_ylim(rho_rng[0], rho_rng[1])
    axs[0].set_ylabel(labels["rho"])
    axs[0].set_title(titles["rho"])

    axs[1].plot(x, vx)
    axs[1].set_ylim(v_rng[0], v_rng[1])
    axs[1].set_ylabel(labels["vx"])
    axs[1].set_title(titles["v-1D"])

    if log_plots["p"]: axs[2].semilogy(x, p)
    else: axs[2].plot(x, p)
    axs[2].set_ylim(p_rng[0], p_rng[1])
    axs[2].set_ylabel(labels["p"])
    axs[2].set_title(titles["p"])


    if log_plots["E"]: axs[3].semilogy(x, E)
    else: axs[3].plot(x, E)
    axs[3].set_ylim(E_rng[0], E_rng[1])
    axs[3].set_ylabel(labels["E"])
    axs[3].set_title(titles["E"])
    
    axs[4].plot(x, vy)
    axs[4].set_ylim(vy_rng[0], vy_rng[1])
    axs[4].set_ylabel(labels["vy"])
    axs[4].set_title(titles["vy"])
    
        
    axs[5].plot(x, By)
    axs[5].set_ylim(By_rng[0], By_rng[1])
    axs[5].set_ylabel(labels["By"])
    axs[5].set_title(titles["By"])

    axs[6].plot(x, vz)
    axs[6].set_ylim(vz_rng[0], vz_rng[1])
    axs[6].set_ylabel(labels["vz"])
    axs[6].set_xlabel(x_label)
    axs[6].set_title(titles["vz"])
    

    axs[7].plot(x, Bz)
    axs[7].set_ylim(Bz_rng[0], Bz_rng[1])
    axs[7].set_ylabel(labels["Bz"])
    axs[7].set_xlabel(x_label)
    axs[7].set_title(titles["Bz"])
    
    title = f"{plot_title} [t = {t:.4f}{time_unit}]"
    fig.suptitle(title, fontsize=30)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    
    png_name = imgFileName1D(n)
    fig.savefig(png_name, dpi=64)
    print(f"Wrote {png_name}")
    plt.close(fig)


rho_rng = range(key_rho)
v_rng = range(key_vx)
vy_rng = range(key_vy)
vz_rng = range(key_vz)
By_rng = range(key_By)
Bz_rng = range(key_Bz)
p_rng = range(key_p)
E_rng = range(key_E)

pad_range(rho_rng)
pad_range(v_rng)
pad_range(p_rng)
pad_range(E_rng)
pad_range(vy_rng)
pad_range(vz_rng)
pad_range(By_rng)
pad_range(Bz_rng)

n = 0
while fileExists(n):
    rho, vx,vy,vz,By,Bz,p, E, t, dx = readFile(n)
    plotFile(n,rho,vx,vy,By,p,E,t,dx, rho_rng, v_rng, p_rng, E_rng, vy_rng, By_rng, vz_rng, Bz_rng)
    n+=1

