'''
CONFIG.py
DRAGONGAZE
Created by Bobbie Markwick on 07/07/2026.
'''



'''
Files
'''
hdf_dir = "/Users/bobbiemarkwick/DRAGON_OUT" #Where to find the HDF5 outputs from DRAGON
img_dir = "/Users/bobbiemarkwick/DRAGON_IMG" #Where to put the plots
h5_base_filename = "frame" #Specifieds HDF5 file names via <base_filename>_#####.h5
img_base_filename = "frame" #Specifieds image file names via <base_filename>_#####.png
#In 2D, image file names will be <prefix>_<base_filename>_#####.png
#In 3D, image file names will be <prefix>_<axis_pair>_<base_filename>_#####.png where <axis_pair> is xy, xz, or yz
fileprefix = {
    "rho": "density",
    "vx": "vx",
    "vy": "vy",
    "vz": "vz",
    "Bx": "Bx",
    "By": "By",
    "Bz": "Bz",
    "p": "pressure",
    "E": "energy"
}

'''
Titles
'''

plot_title = ""#The top level title of your plot

#Subitles for individual plots
titles = {
    "rho": "Density",
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



'''
Axis Labels
'''

#Dragon is Inhernetly dimensionless, though CGS units are easy to add if desired
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
time_unit = ""#optional unit to display after the timestamp

x_label = r"$x$"#\ (\mathrm{cm})$"
y_label = r"$y$"#\ (\mathrm{cm})$"
z_label = r"$z$"#\ (\mathrm{cm})$"

#determines where on the axis the origin should be defined
#Doesn't alter the image, just changes how the axis is labelled
x_mode = 0 #0 origin in middle, +1 = origin at left +x to right, -1 = origin at right +x to left
y_mode = 0 #0 origin in middle, +1 = origin at left +y to right, -1 = origin at right +y to left
z_mode = 0 #0 origin in middle, +1 = origin at left +z to right, -1 = origin at right +z to left


'''
Plot content
'''

#Choose any of the options available in MatPlotLib (https://matplotlib.org/stable/gallery/color/colormap_reference.html)
cmaps = {
    "rho": "viridis",
    "vx": "PRGn",
    "vy": "PRGn",
    "vz": "PRGn",
    "p": "plasma",
    "E": "magma",
    "Bx": "seismic",
    "By": "seismic",
    "Bz": "seismic"
}

log_plots = { #If true, plots this quantity logarithmically. If false, plots it linearly
    "rho": False,
    "vx": False,
    "vy": False,
    "vz": False,
    "Bx": False,
    "By": False,
    "Bz": False,
    "p": True,
    "E": True
}
dpi = 100
