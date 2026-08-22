//
//  Grid.cpp
//  DRAGON/Hydro/ExtendedArray
//
//  Created by Bobbie Markwick on 5/07/2026.
//

#include "Hydro/Grid.hpp"
#include "Config.h"
using namespace DRAGON;

static int validGhosts(int g){
#if defined(MHD) && defined(CTU) //Boundary E fields use transverse fluxes in first ghost layer, so CTU+CT needs a 3rd ghost
    return std::max(g, 3);
#elif defined(MUSCL_Hancock) //MUSCL needs a 2nd ghost to calculate the boundary properly
    return std::max(g, 2);
#else //First order godunov only needs a single ghost
    return std::max(g, 1);
#endif
}

Grid1D::Grid1D(int s_, double dx_, int g_): w(s_, validGhosts(g_)), dx(dx_) { }
PrimitiveState& Grid1D::operator[](int k) { return w[k]; }
const PrimitiveState& Grid1D::operator[](int k) const { return w[k]; }
int Grid1D::getSize() const { return w.getSize(); }
int Grid1D::getGhosts() const { return w.getGhosts(); }

Grid2D::Grid2D(int nx_, int ny_, double dx_, double dy_, int g_):  w(nx_, ny_,validGhosts(g_)),
    #ifdef MHD //B lives on a staggered grid, needs one more face than the body
    B(nx_+1, ny_+1,w.getGhosts()),
    #endif
    dx(dx_), dy(dy_) { }
PrimitiveState& Grid2D::operator[](int i, int j) { return w[i,j]; }
const PrimitiveState& Grid2D::operator[](int i, int j) const { return w[i,j]; }
int Grid2D::getSizeX() const { return w.getSizeX(); }
int Grid2D::getSizeY() const { return w.getSizeY(); }
int Grid2D::getGhosts() const { return w.getGhosts(); }


Grid3D::Grid3D(int nx_, int ny_, int nz_, double dx_, double dy_, double dz_, int g_): w(nx_, ny_, nz_, validGhosts(g_)) ,
    #ifdef MHD //B lives on a staggered grid, needs one more face than the body
    B(nx_+1, ny_+1, nz_+1, w.getGhosts()),
    #endif
    dx(dx_), dy(dy_), dz(dz_) { }
PrimitiveState& Grid3D::operator[](int i, int j, int k) { return w[i,j,k]; }
const PrimitiveState& Grid3D::operator[](int i, int j, int k) const { return w[i,j,k]; }
int Grid3D::getSizeX() const { return w.getSizeX(); }
int Grid3D::getSizeY() const { return w.getSizeY(); }
int Grid3D::getSizeZ() const { return w.getSizeZ(); }
int Grid3D::getGhosts() const { return w.getGhosts(); }
