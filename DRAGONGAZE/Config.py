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
rho_label = r"$\rho$"#\ (\mathrm{g \cdot cm^{-3}})$"
v_label = r"$v_x$"#\ (\mathrm{cm \cdot s^{-1}})$"
vy_label = r"$v_y$"#\ (\mathrm{cm \cdot s^{-1}})$"
vz_label = r"$v_z$"#\ (\mathrm{cm \cdot s^{-1}})$"
By_label = r"$B_y$"
Bz_label = r"$B_z$"
p_label = r"$p$"#\ (\mathrm{g \cdot cm^{-1} s^{-2}})$"
E_label = r"$E$"#\ (\mathrm{g \cdot cm^{2}\cdot s^{-2}})$"
time_unit = ""

rho_title = "Mass Density"
v_title = "Normal Velocity"
p_title = "Pressure"
E_title = "Energy Density"
By_title = "Transverse-y Mangetic Field"
vy_title = "Transverse-y Velocity"
Bz_title = "Transverse-z Mangetic Field"
vz_title = "Transverse-z Velocity"

x_mode = 1 #0 origin in middle, +1 = origin at left +x to right, -1 = origin at right +x to left

plot_title = "Ryu & Jones Shock Tube"


