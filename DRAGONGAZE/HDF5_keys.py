'''
HDF5_keys.py
DRAGONGAZE
Created by Bobbie Markwick on 07/07/2026.

Dataset/attribute key strings for the HDF5 files DRAGONHOARD writes.
These are hand-copied from DRAGON/DRAGONHOARD/HDF5_Attrs.hpp and NOT kept in sync automatically.
If a key is renamed or added on the C++ side, this file must be updated to match or FileUtils.readField() will raise a KeyError.
'''

#Grid
key_fluid = "/fluid";
key_rho = key_fluid + "/rho";
key_p = key_fluid + "/p";
key_E = key_fluid + "/E";

key_v = key_fluid +"/v";
key_vx = key_v + "/x";
key_vy = key_v + "/y";
key_vz = key_v + "/z";

key_mom = key_fluid +"/momentum";
key_px = key_mom + "/x";
key_py = key_mom + "/y";
key_pz = key_mom + "/z";

key_B = key_fluid +"/B";
key_Bx = key_B + "/x";
key_By = key_B + "/y";
key_Bz = key_B + "/z";
key_A = key_fluid +"/A";
key_Ax = key_A + "/x";
key_Ay = key_A + "/y";
key_Az = key_A + "/z";

key_passive = "/passive";

#Metadata
key_fmt = "format_version";
key_wrt_opt = "write_option";
key_dim = "dim";
key_mhd = "MHD";

key_nx = "nx";
key_ny = "ny";
key_nz = "nz";
key_ng = "ng";
key_dx = "dx";
key_dy = "dy";
key_dz = "dz";
key_cyc = "cycle";
key_time = "time";

file_ext = ".h5";


#Maps a plottable field name (as used throughout DRAGONGAZE and Config.py) to its HDF5 dataset path.
#"Bmag" is deliberately absent -- it isn't stored, FileUtils.readField() computes it from Bx/By/Bz instead.
keys = {
    "rho": key_rho,
    "vx": key_vx,
    "vy": key_vy,
    "vz": key_vz,
    "Bx": key_Bx,
    "By": key_By,
    "Bz": key_Bz,
    "p": key_p,
    "E": key_E
}
