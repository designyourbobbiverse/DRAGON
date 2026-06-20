
//
//  Periodic_Tests.cpp
//  DRAGON/Testing/Core/Boundary
//
//  Created by Bobbie Markwick on 20/06/2026.
//

#include "Testing.hpp"
#include "Boundary.hpp"
#include "Grid.hpp"
#include <iostream>

using namespace DRAGON_Test;
using namespace Boundary;

//MARK: Helpers
static vec3 make_tagged_vec(double tag){ return vec3(tag + 100, tag + 200, tag + 300); }

static PrimitiveState G = make_tagged_state(-666);

void fill_1D(Grid1D& grid);
void fill_2D(Grid2D& grid);
void fill_3D(Grid3D& grid);

//MARK: Periodic
void DRAGON_Test::verify_boundary_periodic(bool output){
    if(output) std::cout<<"- 1D: ";
    verify_boundary_periodic_1D();
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- 2D: ";
    verify_boundary_periodic_2D();
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- 3D: ";
    verify_boundary_periodic_3D();
    if(output) std::cout<<"Passed\n";
}
void DRAGON_Test::verify_boundary_periodic_1D() {
    Grid1D grid(4, 1.0, 2);
    fill_1D(grid);

    Periodic(X).apply(grid);
    //Left
    expect_close(grid[-1], grid[3]);
    expect_close(grid[-2], grid[2]);
    //Right
    expect_close(grid[4], grid[0]);
    expect_close(grid[5], grid[1]);
}
//MARK: Periodic - 2D
void DRAGON_Test::verify_boundary_periodic_2D() {
    Grid2D grid(3, 4, 1.0, 1.0, 1);
    
    //X
    fill_2D(grid);
    Periodic(X).apply(grid);
    for (int j = 0; j < grid.getSizeY(); j++) {
        expect_close(grid[-1, j], grid[2, j]);
        expect_close(grid[3, j],  grid[0, j]);
#ifdef MHD
        expect_close(grid.getA()[-1, j], grid.getA()[2, j]);
        expect_close(grid.getA()[3, j],  grid.getA()[0, j]);
        expect_close(grid.getA()[5, j],  grid.getA()[2, j]);
#endif
    }
    //Y
    fill_2D(grid);
    Periodic(Y,false).apply(grid);
    for (int i = 0; i < grid.getSizeX(); i++) {
        expect_close(grid[i,-1], grid[i,3]);
        expect_close(grid[i,4],  grid[i,0]);
#ifdef MHD
        expect_close(grid.getA()[i, -1], grid.getA()[i, 3]);
        expect_close(grid.getA()[i, 4],  grid.getA()[i, 0]);
        expect_close(grid.getA()[i, 6],  grid.getA()[i, 2]);
#endif
    }
    //No corners = no corners
    expect_close(grid[-1,-1],G);
    //Corner
    fill_2D(grid);
    Periodic(X | Y).apply(grid);
    expect_close(grid[-1,-1], grid[2,3]);
}
//MARK: Periodic - 3D
void DRAGON_Test::verify_boundary_periodic_3D() {
    Grid3D grid(3, 4, 5, 1.0, 1.0, 1.0, 1);
    
    //X
    fill_3D(grid);
    Periodic("X",false).apply(grid);
    for (int j = 0; j < grid.getSizeY(); j++) {
        for (int k = 0; k < grid.getSizeZ(); k++) {
            expect_close(grid[-1, j,k], grid[2, j,k]);
            expect_close(grid[3, j,k],  grid[0, j,k]);
#ifdef MHD
            expect_close(grid.getA()[-1, j, k], grid.getA()[2, j, k]);
            expect_close(grid.getA()[3, j, k],  grid.getA()[0, j, k]);
            expect_close(grid.getA()[5, j, k],  grid.getA()[2, j, k]);
#endif
        }
    }
    //No corners = no corners
    expect_close(grid[-1,1,-1],G);
    //Y
    fill_3D(grid);
    Periodic(Y).apply(grid);
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int k = 0; k < grid.getSizeZ(); k++) {
            expect_close(grid[i,-1, k], grid[i,3, k]);
            expect_close(grid[i,4, k],  grid[i,0,k]);
#ifdef MHD
            expect_close(grid.getA()[i, -1, k], grid.getA()[i, 3, k]);
            expect_close(grid.getA()[i, 4, k],  grid.getA()[i, 0, k]);
            expect_close(grid.getA()[i, 6, k],  grid.getA()[i, 2, k]);
#endif
        }
    }
    //Z
    fill_3D(grid);
    Periodic("Z").apply(grid);
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int j = 0; j < grid.getSizeY(); j++) {
            expect_close(grid[i,j,-1], grid[i,j,4]);
            expect_close(grid[i,j,5],  grid[i,j,0]);
#ifdef MHD
            expect_close(grid.getA()[i, j, -1], grid.getA()[i, j, 4]);
            expect_close(grid.getA()[i, j, 5],  grid.getA()[i, j, 0]);
            expect_close(grid.getA()[i, j, 7],  grid.getA()[i, j, 2]);
#endif
        }
    }
    //Corner
    fill_3D(grid);
    Periodic(X | Y | Z).apply(grid);
    expect_close(grid[-1,-1, -1], grid[2,3,4]);
}
