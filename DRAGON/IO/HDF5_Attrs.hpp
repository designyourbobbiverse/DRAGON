//
//  HDF5_Attrs.hpp
//  DRAGON/IO
//
//  Created by Bobbie Markwick on 03/07/2026.
//

//Put all of the attribute names here to ensure consistency
//Grid elements
constexpr std::string key_fluid = "/fluid";
constexpr std::string key_rho = key_fluid + "/rho";
constexpr std::string key_p = key_fluid + "/p";
constexpr std::string key_E = key_fluid + "/E";

constexpr std::string key_v = key_fluid +"/v";
constexpr std::string key_vx = key_v + "/x";
constexpr std::string key_vy = key_v + "/y";
constexpr std::string key_vz = key_v + "/z";

constexpr std::string key_mom = key_fluid +"/momentum";
constexpr std::string key_px = key_mom + "/x";
constexpr std::string key_py = key_mom + "/y";
constexpr std::string key_pz = key_mom + "/z";

constexpr std::string key_B = key_fluid +"/B";
constexpr std::string key_Bx = key_B + "/x";
constexpr std::string key_By = key_B + "/y";
constexpr std::string key_Bz = key_B + "/z";
constexpr std::string key_A = key_fluid +"/A";
constexpr std::string key_Ax = key_A + "/x";
constexpr std::string key_Ay = key_A + "/y";
constexpr std::string key_Az = key_A + "/z";

//Metadata
constexpr std::string key_fmt = "format_version";
constexpr std::string key_wrt_opt = "write_option";
constexpr std::string key_dim = "dim";
constexpr std::string key_mhd = "MHD";

constexpr std::string key_nx = "nx";
constexpr std::string key_ny = "ny";
constexpr std::string key_nz = "nz";
constexpr std::string key_ng = "ng";
constexpr std::string key_dx = "dx";
constexpr std::string key_dy = "dy";
constexpr std::string key_dz = "dz";
constexpr std::string key_cyc = "cycle";
constexpr std::string key_time = "time";

//HDF5 File extension
constexpr std::string file_ext = ".h5";
