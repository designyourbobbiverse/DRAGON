//
//  HDF5_Attrs.hpp
//  DRAGON/IO
//
//  Created by Bobbie Markwick on 03/07/2026.
//
#ifndef HDF5_Attrs_hpp
#define HDF5_Attrs_hpp

#include <string>

//Put all of the attribute names here to ensure consistency
//Grid elements
inline const std::string key_fluid = "/fluid";
inline const std::string key_rho = key_fluid + "/rho";
inline const std::string key_p = key_fluid + "/p";
inline const std::string key_E = key_fluid + "/E";

inline const std::string key_v = key_fluid +"/v";
inline const std::string key_vx = key_v + "/x";
inline const std::string key_vy = key_v + "/y";
inline const std::string key_vz = key_v + "/z";

inline const std::string key_mom = key_fluid +"/momentum";
inline const std::string key_px = key_mom + "/x";
inline const std::string key_py = key_mom + "/y";
inline const std::string key_pz = key_mom + "/z";

inline const std::string key_B = key_fluid +"/B";
inline const std::string key_Bx = key_B + "/x";
inline const std::string key_By = key_B + "/y";
inline const std::string key_Bz = key_B + "/z";
inline const std::string key_A = key_fluid +"/A";
inline const std::string key_Ax = key_A + "/x";
inline const std::string key_Ay = key_A + "/y";
inline const std::string key_Az = key_A + "/z";

//Metadata
inline const std::string key_fmt = "format_version";
inline const std::string key_wrt_opt = "write_option";
inline const std::string key_dim = "dim";
inline const std::string key_mhd = "MHD";

inline const std::string key_nx = "nx";
inline const std::string key_ny = "ny";
inline const std::string key_nz = "nz";
inline const std::string key_ng = "ng";
inline const std::string key_dx = "dx";
inline const std::string key_dy = "dy";
inline const std::string key_dz = "dz";
inline const std::string key_cyc = "cycle";
inline const std::string key_time = "time";

//HDF5 File extension
inline const std::string file_ext = ".h5";


#endif
