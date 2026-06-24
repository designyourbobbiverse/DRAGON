
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

#ifdef MHD
static void expect_outflow_A_z_2D(Grid2D& grid, int ghostI, int ghostJ, int sourceI, int sourceJ, int nextI, int nextJ) {
    double expected = 2 * grid.getA()[sourceI, sourceJ].z - grid.getA()[nextI, nextJ].z;
    assert(approx(grid.getA()[ghostI, ghostJ].z, expected));
}

static void expect_outflow_A_X_3D(Grid3D& grid, int ghostI, int j, int k, int sourceI, int nextI) {
    vec3 expected = 2 * grid.getA()[sourceI, j, k] - grid.getA()[nextI, j, k];
    expected.x = grid.getA()[sourceI, j, k].x;
    expect_close(grid.getA()[ghostI, j, k], expected);
}

static void expect_outflow_A_Y_3D(Grid3D& grid, int i, int ghostJ, int k, int sourceJ, int nextJ) {
    vec3 expected = 2 * grid.getA()[i, sourceJ, k] - grid.getA()[i, nextJ, k];
    expected.y = grid.getA()[i, sourceJ, k].y;
    expect_close(grid.getA()[i, ghostJ, k], expected);
}

static void expect_outflow_A_Z_3D(Grid3D& grid, int i, int j, int ghostK, int sourceK, int nextK) {
    vec3 expected = 2 * grid.getA()[i, j, sourceK] - grid.getA()[i, j, nextK];
    expected.z = grid.getA()[i, j, sourceK].z;
    expect_close(grid.getA()[i, j, ghostK], expected);
}
#endif

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
#ifdef MHD
        expect_outflow_A_z_2D(grid, -1, j, 0, j, 1, j);
        expect_outflow_A_z_2D(grid, 4, j, 3, j, 2, j);
#endif
    }
    //Y
    fill_2D(grid);
    Outflow(Y,false).apply(grid);
    for (int i = 0; i < grid.getSizeX(); i++) {
        expect_close(grid[i,-1], grid[i,0]);
        expect_close(grid[i,4],  grid[i,3]);
#ifdef MHD
        expect_outflow_A_z_2D(grid, i, -1, i, 0, i, 1);
        expect_outflow_A_z_2D(grid, i, 5, i, 4, i, 3);
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
#ifdef MHD
            expect_outflow_A_X_3D(grid, -1, j, k, 0, 1);
            expect_outflow_A_X_3D(grid, 4, j, k, 3, 2);
#endif
        }
    }
    //Y
    fill_3D(grid);
    Outflow(Y_negative).apply(grid);
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int k = 0; k < grid.getSizeZ(); k++) {
            expect_close(grid[i,-1, k], grid[i,0, k]);
            expect_close(grid[i,4, k],  G);
#ifdef MHD
            expect_outflow_A_Y_3D(grid, i, -1, k, 0, 1);
#endif
        }
    }
    fill_3D(grid);
    Outflow(Y_positive).apply(grid);
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int k = 0; k < grid.getSizeZ(); k++) {
            expect_close(grid[i,-1, k], G);
            expect_close(grid[i,4, k],  grid[i,3,k]);
#ifdef MHD
            expect_outflow_A_Y_3D(grid, i, 5, k, 4, 3);
#endif
        }
    }
    //Z
    fill_3D(grid);
    Outflow("Z",false).apply(grid);
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int j = 0; j < grid.getSizeY(); j++) {
            expect_close(grid[i,j,-1], grid[i,j,0]);
            expect_close(grid[i,j,5],  grid[i,j,4]);
#ifdef MHD
            expect_outflow_A_Z_3D(grid, i, j, -1, 0, 1);
            expect_outflow_A_Z_3D(grid, i, j, 6, 5, 4);
#endif
        }
    }
    //No corners = no corners
    expect_close(grid[-1,1,-1],G);
    //Corner
    fill_3D(grid);
    Outflow(X | Y | Z).apply(grid);
    expect_close(grid[-1,-1,-1], grid[0,0,0]);
}

