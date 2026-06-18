//
//  CFL.hpp
//  DRAGON
//
//  Created by Bobbie Markwick on 09/06/2026.
//
#include "Grid.hpp"

namespace CFL{
//Computes the minimum cfl time for the whole grid, with Coefficient applied
double cfl_time(const Grid1D& grid);
double cfl_time(const Grid2D& grid);
double cfl_time(const Grid3D& grid);

//Computes the CFL time for a single cell, Coefficient not applied
double cfl_speed(const PrimitiveState& W, double dx, double dy = 0.0, double dz = 0.0);
//Specific variants of computing the CFL time for a single cell
    //CFL time = highest (|v|+a)/dL of any one side
    double cfl_max_speed(const PrimitiveState& W, double dx, double dy=0.0, double dz = 0.0);
    //CFL time = sum of (|v|+a)/dL over all applicable sides
    double cfl_add_speed(const PrimitiveState& W, double dx, double dy=0.0, double dz = 0.0);
    //CFL time = [sum of ((|v|+a)/dL)^p]^(1/p)
    double cfl_pow_speed(const PrimitiveState& W, double p, double dx, double dy=0.0, double dz = 0.0);
}
