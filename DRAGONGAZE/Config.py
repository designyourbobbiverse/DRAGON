'''
CONFIG.py
DRAGONGAZE
Created by Bobbie Markwick on 07/07/2026.
'''

hdf_dir = "/Users/bobbiemarkwick/DRAGON_OUT" #Where to find the HDF5 outputs from DRAGON
img_dir = "/Users/bobbiemarkwick/DRAGON_IMG" #Where to put the plots
h5_base_filename = "frame" #Specifieds HDF5 file names via <base_filename>_#####.h5
img_base_filename = "frame" #Specifieds image file names via <base_filename>_#####.png

log_plots = False #Set True to plot density, pressure, and energy density logarithmically

#Dragon is Inhernetly dimensionless, though CGS units are easy to add if desired
x_label = r"$x$"#\ (\mathrm{cm})$"
labels = {
    "rho": r"$\rho$",#\ (\mathrm{g \cdot cm^{-3}})$",
    "vx": r"$v_x$",#\ (\mathrm{cm \cdot s^{-1}})$",
    "vy": r"$v_y$",#\ (\mathrm{cm \cdot s^{-1}})$",
    "vz": r"$v_z$",#\ (\mathrm{cm \cdot s^{-1}})$",
    "Bx": r"$B_x$",
    "By": r"$B_y$",
    "Bz": r"$B_z$",
    "p": r"$p$",#\ (\mathrm{g \cdot cm^{-1} s^{-2}})$"
    "E": r"$E$"#\ (\mathrm{g \cdot cm^{2}\cdot s^{-2}})$",
}
time_unit = ""

titles = {
    "rho": "Mass Density",
    "v-1D": "Normal Velocity",
    "vx": "Velocity (x)",
    "p": "Pressure",
    "E": "Energy Density",
    "Bx": "Mangetic Field (x)",
    "By": "Mangetic Field (y)",
    "vy": "Velocity (y)",
    "Bz": "Mangetic Field (z)",
    "vz": "Velocity (z)",
}

x_mode = 1 #0 origin in middle, +1 = origin at left +x to right, -1 = origin at right +x to left

plot_title = "Ryu & Jones Shock Tube"

