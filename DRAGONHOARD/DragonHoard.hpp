//
//  DRAGONHOARD.hpp
//  DRAGONHOARD
//
//  Created by Bobbie Markwick on 03/07/2026.
//

#ifndef HDF5Output_hpp
#define HDF5Output_hpp

#include "DRAGONHOARD_Config.h" //To give users of this module access to the output_base_name and output_dir settings

#include "Grid.hpp"
#include <string>

namespace DRAGONHOARD{

//Write the grid (and other metadata) to filename.h5
void writeToFile(DRAGON::Grid& grid, double t, int cycle, const std::string& filename);
void writeToFile(DRAGON::Grid1D& grid, double t, int cycle, const std::string& filename);
void writeToFile(DRAGON::Grid2D& grid, double t, int cycle, const std::string& filename);
void writeToFile(DRAGON::Grid3D& grid, double t, int cycle, const std::string& filename);

//Load the data in filename.h5 into grid
void loadFromFile(DRAGON::Grid& grid, double& t, int& cycle, const std::string& filename);
void loadFromFile(DRAGON::Grid1D& grid, double& t, int& cycle, const std::string& filename);
void loadFromFile(DRAGON::Grid2D& grid, double& t, int& cycle, const std::string& filename);
void loadFromFile(DRAGON::Grid3D& grid, double& t, int& cycle, const std::string& filename);

//Verify that the output directory exists, if not create it
void verifyOutputDirectory();

//Identify the filename for the restart frame
std::string restartFileName();

//Convert cycle integer to a 5-character string
std::string cycle_string(int n);

}


#endif /* HDF5Output_hpp */
