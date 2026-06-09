//
//  CFL.hpp
//  DRAGON
//
//  Created by Bobbie Markwick on 09/06/2026.
//
#include "Grid.hpp"

namespace CFL{
double cfl_time(const Grid1D& grid);
double cfl_time(const Grid2D& grid);
double cfl_time(const Grid3D& grid);

double cfl_speed(const PrimitiveState& W, double dx, double dy = 0.0, double dz = 0.0);

double cfl_max_speed(const PrimitiveState& W, double dx, double dy=0.0, double dz = 0.0);
double cfl_add_speed(const PrimitiveState& W, double dx, double dy=0.0, double dz = 0.0);
double cfl_pow_speed(const PrimitiveState& W, double p, double dx, double dy=0.0, double dz = 0.0);
}
