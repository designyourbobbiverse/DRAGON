//
//  DRAGONHOARD.hpp
//  DRAGONHOARD
//
//  Created by Bobbie Markwick on 03/07/2026.
//

#ifndef HDF5Output_hpp
#define HDF5Output_hpp

#include "Grid.hpp"
#include <string>

namespace DRAGONHOARD{
    

void writeToFile(Grid& grid, double t, int cycle, const std::string& filename);
void writeToFile(Grid1D& grid, double t, int cycle, const std::string& filename);
void writeToFile(Grid2D& grid, double t, int cycle, const std::string& filename);
void writeToFile(Grid3D& grid, double t, int cycle, const std::string& filename);

void loadFromFile(Grid& grid, double& t, int& cycle, const std::string& filename);
void loadFromFile(Grid1D& grid, double& t, int& cycle, const std::string& filename);
void loadFromFile(Grid2D& grid, double& t, int& cycle, const std::string& filename);
void loadFromFile(Grid3D& grid, double& t, int& cycle, const std::string& filename);

void verifyOutputDirectory();
std::string restartFileName();

std::string cycle_string(int n);

}


#endif /* HDF5Output_hpp */
