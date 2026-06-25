
//
//  Outflow_Tests.cpp
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


void DRAGON_Test::verify_boundary_outflow(bool output){
    if(output) std::cout<<"- 1D: ";
    verify_boundary_outflow_1D();
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- 2D: ";
    verify_boundary_outflow_2D();
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- 3D: ";
    verify_boundary_outflow_3D();
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- Gated 1D: ";
    verify_boundary_outflow_1D_gated();
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- Gated 2D: ";
    verify_boundary_outflow_2D_gated();
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- Gated 3D: ";
    verify_boundary_outflow_3D_gated_X();
    verify_boundary_outflow_3D_gated_Y();
    verify_boundary_outflow_3D_gated_Z();
    if(output) std::cout<<"Passed\n";
}



//MARK: Helpers
static PrimitiveState G = make_tagged_state(-666);


void fill_1D(Grid1D& grid);
void fill_2D(Grid2D& grid);
void fill_3D(Grid3D& grid);

//MARK: Outflow - 1D
void DRAGON_Test::verify_boundary_outflow_1D() {
    Grid1D grid(4, 1.0, 2);
    fill_1D(grid);

    Outflow(X).apply(grid);
    //Left
    expect_close(grid[-1], grid[0]);
    expect_close(grid[-2], grid[0]);
    //Right
    expect_close(grid[4], grid[3]);
    expect_close(grid[5], grid[3]);
    //Single face filling: Left
    fill_1D(grid);
    Outflow("X-").apply(grid);
    expect_close(grid[-1], grid[0]);
    expect_close(grid[-2], grid[0]);
    expect_close(grid[4], G);
    expect_close(grid[5], G);
    //Single face filling: Right
    fill_1D(grid);
    Outflow("X+").apply(grid);
    expect_close(grid[4], grid[3]);
    expect_close(grid[5], grid[3]);
    expect_close(grid[-1], G);
    expect_close(grid[-2], G);
}
//MARK: Outflow - 2D
void DRAGON_Test::verify_boundary_outflow_2D() {
    Grid2D grid(3, 4, 1.0, 1.0, 1);
    
    //X
    fill_2D(grid);
    Outflow(X).apply(grid);
    for (int j = 0; j < grid.getSizeY(); j++) {
        expect_close(grid[-1, j], grid[0, j]);
        expect_close(grid[3, j],  grid[2, j]);
        #ifdef MHD //Check Transverse Magnetic Fields
        assert(magneticY(grid, -1,j) == magneticY(grid, 0,j));
        assert(magneticY(grid, 4,j) == magneticY(grid, 3,j));
        #endif
    }
    #ifdef MHD //A field Corners
    int ng = grid.getGhosts(), nx = grid.getSizeX(), ny = grid.getSizeY();
    vec3 expected = 2*grid.getA()[nx+ng-1,ny+ng] - grid.getA()[nx+ng-2,ny+ng];
    expected.x = grid.getA()[nx+ng-1,ny+ng].x;
    expect_close(grid.getA()[nx+ng,ny+ng], expected);
    #endif
    //Y
    fill_2D(grid);
    Outflow(Y,false).apply(grid);
    for (int i = 0; i < grid.getSizeX(); i++) {
        expect_close(grid[i,-1], grid[i,0]);
        expect_close(grid[i,4],  grid[i,3]);
        #ifdef MHD //Check Transverse Magnetic Fields
        assert(magneticX(grid, i,-1) == magneticX(grid, i,0));
        assert(magneticX(grid, i,5) == magneticX(grid, i,4));
        #endif
    }
    //No corners = no corners
    expect_close(grid[-1,-1],G);
    //Corner
    fill_2D(grid);
    Outflow("XY").apply(grid);
    expect_close(grid[-1,-1], grid[0,0]);
}
//MARK: Outflow - 3D
void DRAGON_Test::verify_boundary_outflow_3D() {
    Grid3D grid(3, 4, 5, 1.0, 1.0, 1.0, 1);
    
    //X
    fill_3D(grid);
    Outflow(X).apply(grid);
    for (int j = 0; j < grid.getSizeY(); j++) {
        for (int k = 0; k < grid.getSizeZ(); k++) {
            expect_close(grid[-1, j,k], grid[0, j,k]);
            expect_close(grid[3, j,k],  grid[2, j,k]);
            #ifdef MHD //Check Transverse Magnetic Fields
            assert(magneticY(grid, -1,j,k) == magneticY(grid, 0,j,k));
            assert(magneticZ(grid, -1,j,k) == magneticZ(grid, 0,j,k));
            assert(magneticY(grid, 3,j,k) == magneticY(grid, 2,j,k));
            assert(magneticZ(grid, 4,j,k) == magneticZ(grid, 3,j,k));
            #endif
        }
    }
    #ifdef MHD //A field corner check
    int ng = grid.getGhosts(), nx = grid.getSizeX(), ny = grid.getSizeY(), nz = grid.getSizeZ();
    vec3 expected = 2*grid.getA()[nx+ng-1,ny+ng,nz+ng] - grid.getA()[nx+ng-2,ny+ng,nz+ng];
    expected.x = grid.getA()[nx+ng-1,ny+ng,nz+ng].x;
    expect_close(grid.getA()[nx+ng,ny+ng,nz+ng], expected);
    #endif
    //Y
    fill_3D(grid);
    Outflow(Y_negative).apply(grid);
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int k = 0; k < grid.getSizeZ(); k++) {
            expect_close(grid[i,-1, k], grid[i,0, k]);
            expect_close(grid[i,4, k],  G);
            #ifdef MHD //Check Transverse Magnetic Fields
            assert(magneticX(grid, i,-1,k) == magneticX(grid, i,0,k));
            assert(magneticZ(grid, i,-1,k) == magneticZ(grid, i,0,k));
            #endif
        }
    }
    fill_3D(grid);
    Outflow(Y_positive).apply(grid);
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int k = 0; k < grid.getSizeZ(); k++) {
            expect_close(grid[i,-1, k], G);
            expect_close(grid[i,4, k],  grid[i,3,k]);
            #ifdef MHD //Check Transverse Magnetic Fields
            assert(magneticX(grid, i,4,k) == magneticX(grid, i,3,k));
            assert(magneticZ(grid, i,4,k) == magneticZ(grid, i,3,k));
            #endif
        }
    }
    #ifdef MHD//A field corner check
    expected = 2*grid.getA()[nx+ng,ny+ng-1,nz+ng] - grid.getA()[nx+ng,ny+ng-2,nz+ng];
    expected.y = grid.getA()[nx+ng,ny+ng-1,nz+ng].y;
    expect_close(grid.getA()[nx+ng,ny+ng,nz+ng], expected);
    #endif
    //Z
    fill_3D(grid);
    Outflow("Z").apply(grid);
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int j = 0; j < grid.getSizeY(); j++) {
            expect_close(grid[i,j,-1], grid[i,j,0]);
            expect_close(grid[i,j,5],  grid[i,j,4]);
            #ifdef MHD //Check Transverse Magnetic Fields
            assert(magneticX(grid, i,j,-1) == magneticX(grid, i,j,0));
            assert(magneticY(grid, i,j,-1) == magneticY(grid, i,j,0));
            assert(magneticX(grid, i,j,5) == magneticX(grid, i,j,4));
            assert(magneticY(grid, i,j,5) == magneticY(grid, i,j,4));

            #endif
        }
    }
    #ifdef MHD//A field corner check
    expected = 2*grid.getA()[nx+ng,ny+ng,nz+ng-1] - grid.getA()[nx+ng,ny+ng,nz+ng-2];
    expected.z = grid.getA()[nx+ng,ny+ng,nz+ng-1].z;
    expect_close(grid.getA()[nx+ng,ny+ng,nz+ng], expected);
    #endif
    //No corners = no corners
    Outflow("Z",false).apply(grid);
    expect_close(grid[-1,1,-1],G);
    //Corner
    fill_3D(grid);
    Outflow(X | Y | Z).apply(grid);
    expect_close(grid[-1,-1,-1], grid[0,0,0]);
}

