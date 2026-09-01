'''
Config.py
DRAGONGAZE
Created by Bobbie Markwick on 07/07/2026.

User-facing settings for the DRAGONGAZE plotting scripts: file locations, plot
titles/labels, colormaps, and axis conventions.

Each of the five dictionaries below (fileprefix, titles, labels, cmaps,
log_plots) is keyed by field name ("rho", "vx", "Bmag", etc). To add a new
plottable field, add a matching entry to all five here AND a corresponding
HDF5 dataset key in HDF5_keys.py's `keys` dict (see that file's header) --
GAZE2D.py/GAZE3D.py are field-generic and will pick it up automatically.
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
    "Bmag": "B",
    "p": "pressure",
    "E": "energy"
}

'''
Titles
'''

plot_title = ""#The top level title of your plot

#Subtitles for individual plots.
titles = {
    "rho": "Density",
    "v-1D": "Normal Velocity", #used for vx in GAZE1D.py/GAZE1DMHD.py
    "vx": "Velocity (x)",
    "p": "Pressure",
    "E": "Energy Density",
    "Bx": "Mangetic Field (x)",
    "By": "Mangetic Field (y)",
    "vy": "Velocity (y)",
    "Bz": "Mangetic Field (z)",
    "vz": "Velocity (z)",
    "Bmag": "Mangetic Field Strength", #|B| = sqrt(Bx^2+By^2+Bz^2)
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
    "Bmag": r"$|B|$",
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
    "Bz": "seismic",
    "Bmag": "plasma",
    "default": "viridis"
}

#If true, plots this quantity logarithmically (via LogNorm/semilogy). If false, plots it linearly.
#Note: In 1D, only rho/p/E respect this (velocity and transverse-B panels are always linear in 1D).
log_plots = {
    "rho": False,
    "vx": False,
    "vy": False,
    "vz": False,
    "Bx": False,
    "By": False,
    "Bz": False,
    "Bmag": False,
    "p": True,
    "E": True
}
dpi = 100
