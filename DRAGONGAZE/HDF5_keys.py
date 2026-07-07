'''
HDF5_keys.py
DRAGONGAZE
Created by Bobbie Markwick on 07/07/2026.
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
