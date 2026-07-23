//
//
//  DRAGONHOARD_Config.h
//  DRAGONHOARD
//
//  Created by Bobbie Markwick on 23/07/2026.
//

#ifndef DRAGONHOARD_CONFIG_h
#define DRAGONHOARD_CONFIG_h

#include <string>


namespace DRAGONHOARD{
inline std::string output_base_name = "frame"; //File names will be <output_base_name>_#####.hdf
inline std::string output_dir = "/Users/bobbiemarkwick/DRAGON_OUT"; //Set this to your output directory

#define RESTART_FROM_FILE //Attempt to restart from the output of a previous run
    #define RESTART_FRAME 0 //Use a number <0 to automatically find the latest file and restart from that

//Output parameters
#define HDF5_COMPRESSION_LEVEL 4 //HDF5 compression level. Set to a nonpositive number for no compression

    #define HDF5_WRITE_PRIMITIVE 1 //Only write the primitive values to the output files
    #define HDF5_WRITE_CONSERVATIVE 6 //Only write the conservative values to the output files
    #define HDF5_WRITE_PRIMITIVE_AND_ENERGY 3 //Write the primitive values, plus the energy density
    #define HDF5_WRITE_PRIMITIVE_AND_CONSERVATIVE 7 //Write both Primitive & Conservative values (produces larger files)
#define HDF5_WRITE_OPTION HDF5_WRITE_PRIMITIVE_AND_ENERGY //Determines which fluid representation to use for output files

//Writes ghost cells to file. Useful for debugging, but typically you want to disable this
//#define WRITE_GHOSTS_TO_FILE

    #define HDF5_WRITE_OMIT 0 //Omits B values if they can be calculated from A
    #define HDF5_WRITE_FLOAT 1 //Writes B values as Floats if they can be calculated from A
    #define HDF5_WRITE_DOUBLE 2 //Writes all B values as Doubles
#define HDF5_REDUNDANT_VALS_OPTION HDF5_WRITE_FLOAT //Output precision of certain calculable values
    //Specifically, applies to B (2D+3D only), and E (HDF5_WRITE_PRIMITIVE_AND_ENERGY only)
}

#endif
