#
#  GAZE1D.py
#  DRAGON
#
#  Created by Bobbie Markwick on 07/07/2026.
#


import h5py

hdf_dir = "/Users/bobbiemarkwick/DRAGON_OUT"
img_dir = "/Users/bobbiemarkwick/DRAGON_IMG"
base_filename = "frame"


with h5py.File(hdf_dir + "/" + base_filename + "_00000.h5", "r") as f:
    print(list(f.keys()))          # top-level groups/datasets
    print(dict(f.attrs))           # file-level attributes

    #rho = f["rho"][:]              # read dataset into NumPy array
    #vx  = f["vx"][:]
    #p   = f["p"][:]

    t = f.attrs["time"]
    cycle = f.attrs["cycle"]

