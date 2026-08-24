//  DRAGON/main/init.cpp
//  Created by Bobbie Markwick
//

#include "main/Problem.hpp" //The user's problem

#include "Config.h"
#include "DragonHoard.hpp"  //File restart
#include <iostream>         //Console messages
#include "MHD/CT.hpp"       //Magnetic Field initialization
using namespace DRAGON;

#ifndef TESTMODE
void Problem::load(Grid& problem, double& time, int& cycle){
    #ifdef RESTART_FROM_FILE
    std::string file = DRAGONHOARD::restartFileName();
    if (file.size() > 0) {
        try{
            DRAGONHOARD::loadFromFile(problem, time, cycle, file);
            return;
        } catch(std::exception& e) {
            std::cerr << e.what()<<"\nInitializing from scratch\n";
        }
    }
    #endif
    Grid1D* grid1D = dynamic_cast<Grid1D*>(&problem);
    Grid2D* grid2D = dynamic_cast<Grid2D*>(&problem);
    Grid3D* grid3D = dynamic_cast<Grid3D*>(&problem);
    if (grid1D) Problem::init(*grid1D);
    else if (grid2D) Problem::init(*grid2D);
    else if (grid3D) Problem::init(*grid3D);
    Problem::completeProblemInit(problem);

    DRAGONHOARD::writeToFile(problem, 0, 0, DRAGONHOARD::output_base_name + "_" + DRAGONHOARD::cycle_string(0));

}


//MARK: Problem initialization
void Problem::init(Grid1D& grid){
    const int nx = grid.getSize();
    const double dx = grid.dx;
    for (int i=0; i<nx; i++) {
        grid[i] = Problem::initialFluidState((i+0.5)*dx, 0, 0);
    }
    //1D MHD includes B on the cells, no need for special CT array
}
void Problem::init(Grid2D& grid){
    const int nx = grid.getSizeX(), ny = grid.getSizeY();
    const double dx = grid.dx, dy = grid.dy;
    for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++) {
            grid[i,j] = Problem::initialFluidState((i+0.5)*dx, (j+0.5)*dy, 0);
        }
    }
    #ifdef MHD
    ExtendedArray2D<vec3> A(nx+1, ny+1, 0);
    for (int i=0; i<=nx; i++) {
        for (int j=0; j<=ny; j++) {
            grid._B()[i,j] = {0,0,0};
            A[i,j].x = Problem::initialMagneticPotential((i+0.5)*dx, j*dy, 0).x;
            A[i,j].y = Problem::initialMagneticPotential(i*dx, (j+0.5)*dy, 0).y;
            A[i,j].z = Problem::initialMagneticPotential(i*dx, j*dy, 0).z;
        }
    }
    CT::Faraday(A, grid._B(), -1/dx, -1/dy, 0);
    grid.initialize_B_fields();
    #endif
}
void Problem::init(Grid3D& grid){
    const int nx = grid.getSizeX(), ny = grid.getSizeY(), nz = grid.getSizeZ();
    const double dx = grid.dx, dy = grid.dy, dz = grid.dz;
    for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++) {
            for (int k=0; k<nz; k++) {
                grid[i,j,k] = Problem::initialFluidState((i+0.5)*dx, (j+0.5)*dy, (k+0.5)*dz);
            }
        }
    }
    #ifdef MHD
    ExtendedArray3D<vec3> A(nx+1, ny+1, nz+1, 0);
    for (int i=0; i<=nx; i++) {
        for (int j=0; j<=ny; j++) {
            for (int k=0; k<=nz; k++) {
                grid._B()[i,j,k] = {0,0,0};
                A[i,j,k].x = Problem::initialMagneticPotential((i+0.5)*dx, j*dy, k*dz).x;
                A[i,j,k].y = Problem::initialMagneticPotential(i*dx, (j+0.5)*dy, k*dz).y;
                A[i,j,k].z = Problem::initialMagneticPotential(i*dx, j*dy, (k+0.5)*dz).z;
            }
        }
    }
    CT::Faraday(A, grid._B(), -1/dx, -1/dy, -1/dz, 0);
    grid.initialize_B_fields();
    #endif
}


#endif
