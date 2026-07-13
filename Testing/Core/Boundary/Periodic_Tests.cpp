
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
#ifdef MHD
using namespace MagneticGrid;
#endif

//MARK: Helpers

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
    Grid2D grid(3, 4, 1.0, 1.0, 2);
    
    //X
    fill_2D(grid);
    Periodic(X).apply(grid);
    for (int j = 0; j < grid.getSizeY(); j++) {
        expect_close(grid[-1, j], grid[2, j]);
        expect_close(grid[3, j],  grid[0, j]);
        #ifdef MHD
        expect_close(magneticX(grid, -1, j), magneticX(grid, 2, j));
        expect_close(magneticY(grid, -1, j), magneticY(grid, 2, j));

        expect_close(magneticX(grid, 3, j), magneticX(grid, 0, j));
        expect_close(magneticY(grid, 3, j), magneticY(grid, 0, j));

        expect_close(magneticX(grid, 4, j), magneticX(grid, 1, j));
        expect_close(magneticY(grid, 4, j), magneticY(grid, 1, j));

        #endif
    }
    #ifdef MHD//A field corner check
    int ng = grid.getGhosts(), nx = grid.getSizeX(), ny = grid.getSizeY();
    vec3 dA = grid._A()[nx,ny+ng] - grid._A()[0,ny+ng];
    expect_close(grid._A()[nx+ng,ny+ng], grid._A()[ng,ny+ng] + dA);
    #endif
    //Y
    fill_2D(grid);
    Periodic(Y).apply(grid);
    for (int i = 0; i < grid.getSizeX(); i++) {
        expect_close(grid[i,-1], grid[i,3]);
        expect_close(grid[i,4],  grid[i,0]);
        #ifdef MHD
        expect_close(magneticX(grid, i,-1), magneticX(grid, i,3));
        expect_close(magneticY(grid, i,-1), magneticY(grid, i,3));

        expect_close(magneticX(grid, i,4), magneticX(grid, i,0));
        expect_close(magneticY(grid, i,4), magneticY(grid, i,0));

        expect_close(magneticX(grid, i,5), magneticX(grid, i,1));
        expect_close(magneticY(grid, i,5), magneticY(grid, i,1));
        #endif
    }
    //No corners = no corners
    fill_2D(grid);
    Periodic(Y,false).apply(grid);
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
    Periodic("X").apply(grid);
    for (int j = 0; j < grid.getSizeY(); j++) {
        for (int k = 0; k < grid.getSizeZ(); k++) {
            expect_close(grid[-1, j,k], grid[2, j,k]);
            expect_close(grid[3, j,k],  grid[0, j,k]);
            #ifdef MHD
            expect_close(magneticX(grid, -1, j,k), magneticX(grid, 2, j,k));
            expect_close(magneticY(grid, -1, j,k), magneticY(grid, 2, j,k));
            expect_close(magneticZ(grid, -1, j,k), magneticZ(grid, 2, j,k));

            expect_close(magneticX(grid, 3, j,k), magneticX(grid, 0, j,k));
            expect_close(magneticY(grid, 3, j,k), magneticY(grid, 0, j,k));
            expect_close(magneticZ(grid, 3, j,k), magneticZ(grid, 0, j,k));

            expect_close(magneticX(grid, 4, j,k), magneticX(grid, 1, j,k));
            expect_close(magneticY(grid, 4, j,k), magneticY(grid, 1, j,k));
            expect_close(magneticZ(grid, 4, j,k), magneticZ(grid, 1, j,k));
            #endif
        }
    }
    #ifdef MHD//A field corner check
    int ng = grid.getGhosts(), nx = grid.getSizeX(), ny = grid.getSizeY(), nz = grid.getSizeZ();
    vec3 dA = grid._A()[nx,ny+ng,nz+ng] - grid._A()[0,ny+ng,nz+ng];
    expect_close(grid._A()[nx+ng,ny+ng,nz+ng], grid._A()[ng,ny+ng,nz+ng] + dA);
    #endif
    //No corners = no corners
    fill_3D(grid);
    Periodic(X).apply(grid);
    expect_close(grid[-1,1,-1],G);
    //Y
    fill_3D(grid);
    Periodic(Y).apply(grid);
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int k = 0; k < grid.getSizeZ(); k++) {
            expect_close(grid[i,-1, k], grid[i,3, k]);
            expect_close(grid[i,4, k],  grid[i,0,k]);
            #ifdef MHD
            expect_close(magneticX(grid, i,-1,k), magneticX(grid, i,3,k));
            expect_close(magneticY(grid, i,-1,k), magneticY(grid, i,3,k));
            expect_close(magneticZ(grid, i,-1,k), magneticZ(grid, i,3,k));

            expect_close(magneticX(grid, i,4,k), magneticX(grid, i,0,k));
            expect_close(magneticY(grid, i,4,k), magneticY(grid, i,0,k));
            expect_close(magneticZ(grid, i,4,k), magneticZ(grid, i,0,k));

            expect_close(magneticX(grid, i,5,k), magneticX(grid, i,1,k));
            expect_close(magneticY(grid, i,5,k), magneticY(grid, i,1,k));
            expect_close(magneticZ(grid, i,5,k), magneticZ(grid, i,1,k));
            #endif
        }
    }
    #ifdef MHD//A field corner check
    dA = grid._A()[nx+ng,ny,nz+ng] - grid._A()[nx+ng,0,nz+ng];
    expect_close(grid._A()[nx+ng,ny+ng,nz+ng], grid._A()[nx+ng,ng,nz+ng] + dA);
    #endif
    //Z
    fill_3D(grid);
    Periodic("Z").apply(grid);
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int j = 0; j < grid.getSizeY(); j++) {
            expect_close(grid[i,j,-1], grid[i,j,4]);
            expect_close(grid[i,j,5],  grid[i,j,0]);
            #ifdef MHD
            expect_close(magneticX(grid, i,j,-1), magneticX(grid, i,j,4));
            expect_close(magneticY(grid, i,j,-1), magneticY(grid, i,j,4));
            expect_close(magneticZ(grid, i,j,-1), magneticZ(grid, i,j,4));

            expect_close(magneticX(grid, i,j,5), magneticX(grid, i,j,0));
            expect_close(magneticY(grid, i,j,5), magneticY(grid, i,j,0));
            expect_close(magneticZ(grid, i,j,5), magneticZ(grid, i,j,0));

            expect_close(magneticX(grid, i,j,6), magneticX(grid, i,j,1));
            expect_close(magneticY(grid, i,j,6), magneticY(grid, i,j,1));
            expect_close(magneticZ(grid, i,j,6), magneticZ(grid, i,j,1));
            #endif
        }
    }
    #ifdef MHD//A field corner check
    dA = grid._A()[nx+ng,ny+ng,nz] - grid._A()[nx+ng,ny+ng,0];
    expect_close(grid._A()[nx+ng,ny+ng,nz+ng], grid._A()[nx+ng,ny+ng,ng] + dA);
    #endif
    //Corner
    fill_3D(grid);
    Periodic(X | Y | Z).apply(grid);
    expect_close(grid[-1,-1, -1], grid[2,3,4]);
}
