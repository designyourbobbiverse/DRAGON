
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



void DRAGON_Test::verify_boundary_fixed(bool output){
    if(output) std::cout<<"- 1D: ";
    verify_boundary_fixed_1D();
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- 2D: ";
    verify_boundary_fixed_2D();
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- 3D: ";
    verify_boundary_fixed_3D();
    if(output) std::cout<<"Passed\n";
}

//MARK: Helpers
static PrimitiveState G = make_tagged_state(-666);

#ifdef MHD
static double magneticX(Grid2D& grid, int i, int j) {
    return (grid.getA()[i, j + 1].z - grid.getA()[i, j].z) / grid.dy;
}

static double magneticY(Grid2D& grid, int i, int j) {
    return (grid.getA()[i, j].z - grid.getA()[i + 1, j].z) / grid.dx;
}

static double magneticX(Grid3D& grid, int i, int j, int k) {
    return (grid.getA()[i, j + 1, k].z - grid.getA()[i, j, k].z) / grid.dy
         - (grid.getA()[i, j, k + 1].y - grid.getA()[i, j, k].y) / grid.dz;
}

static double magneticY(Grid3D& grid, int i, int j, int k) {
    return (grid.getA()[i, j, k + 1].x - grid.getA()[i, j, k].x) / grid.dz
         - (grid.getA()[i + 1, j, k].z - grid.getA()[i, j, k].z) / grid.dx;
}

static double magneticZ(Grid3D& grid, int i, int j, int k) {
    return (grid.getA()[i + 1, j, k].y - grid.getA()[i, j, k].y) / grid.dx
         - (grid.getA()[i, j + 1, k].x - grid.getA()[i, j, k].x) / grid.dy;
}

#endif

void fill_1D(Grid1D& grid);
void fill_2D(Grid2D& grid);
void fill_3D(Grid3D& grid);

//MARK: Fixed - 1D
void DRAGON_Test::verify_boundary_fixed_1D() {
    Grid1D grid(4, 1.0, 2);
    PrimitiveState W = make_tagged_state(2112);
    
    fill_1D(grid);

    Boundary::Fixed(W,X).apply(grid);
    //Left
    expect_close(grid[-1], W);
    expect_close(grid[-2], W);
    //Right
    expect_close(grid[4], W);
    expect_close(grid[5], W);
    //Single face filling: Left
    fill_1D(grid);
    Boundary::Fixed(W,"X-").apply(grid);
    expect_close(grid[-1], W);
    expect_close(grid[-2], W);
    expect_close(grid[4], G);
    expect_close(grid[5], G);
    //Single face filling: Right
    fill_1D(grid);
    Boundary::Fixed(W,"X+").apply(grid);
    expect_close(grid[-1], G);
    expect_close(grid[-2], G);
    expect_close(grid[4], W);
    expect_close(grid[5], W);
}
//MARK: Fixed - 2D
void DRAGON_Test::verify_boundary_fixed_2D() {
    Grid2D grid(3, 4, 1.0, 1.0, 1);
    PrimitiveState W = make_tagged_state(2112);

    //X
    fill_2D(grid);
    Boundary::Fixed(W,X,false).apply(grid);
    for (int j = 0; j < grid.getSizeY(); j++) {
        expect_close(grid[-1, j], W);
        expect_close(grid[3, j],  W);
#ifdef MHD
        approx(magneticY(grid, -1, j), W.B.y);
        approx(magneticY(grid, 3, j), W.B.y);
#endif
    }
    //No Corners = No Corners
    expect_close(grid[-1,-1],G);
    //Y
    fill_2D(grid);
    Boundary::Fixed(W,"Y").apply(grid);
    for (int i = 0; i < grid.getSizeX(); i++) {
        expect_close(grid[i,-1], W);
        expect_close(grid[i,4],  W);
#ifdef MHD
        approx(magneticX(grid, i, -1), W.B.x);
        approx(magneticX(grid, i, 4), W.B.x);
#endif
    }
    //Corner
    fill_2D(grid);
    Boundary::Fixed(W,X | Y).apply(grid);
    expect_close(grid[-1,-1], W);
}
//MARK: Fixed - 3D
void DRAGON_Test::verify_boundary_fixed_3D() {
    Grid3D grid(3, 4, 5, 1.0, 1.0, 1.0, 1);
    PrimitiveState W = make_tagged_state(2112);

    //X
    fill_3D(grid);
    Boundary::Fixed(W,"X").apply(grid);
    for (int j = 0; j < grid.getSizeY(); j++) {
        for (int k = 0; k < grid.getSizeZ(); k++) {
            expect_close(grid[-1, j,k], W);
            expect_close(grid[3, j,k],  W);
#ifdef MHD
            approx(magneticY(grid, -1,j,k), W.B.y);
            approx(magneticZ(grid, -1,j,k), W.B.z);
            approx(magneticY(grid, 3,j,k), W.B.y);
            approx(magneticZ(grid, 3,j,k), W.B.z);
#endif
        }
    }
    //Y
    fill_3D(grid);
    Boundary::Fixed(W,"Y",false).apply(grid);
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int k = 0; k < grid.getSizeZ(); k++) {
            expect_close(grid[i,-1, k], W);
            expect_close(grid[i,4, k],  W);
#ifdef MHD
            if (k + 1 < grid.getSizeZ()) {
                approx(magneticX(grid, i,-1,k), W.B.x);
                approx(magneticX(grid, i,4,k), W.B.x);
            }
            if (i + 1 < grid.getSizeX()) {
                approx(magneticZ(grid, i,-1,k), W.B.z);
                approx(magneticZ(grid, i,4,k), W.B.z);
            }
#endif
        }
    }
    //No corners = no corners
    expect_close(grid[-1,-1,1],G);
    //Z
    fill_3D(grid);
    Boundary::Fixed(W,"Z").apply(grid);
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int j = 0; j < grid.getSizeY(); j++) {
            expect_close(grid[i,j,-1], W);
            expect_close(grid[i,j,5],  W);
#ifdef MHD
            approx(magneticX(grid, i, j, -1), W.B.x);
            approx(magneticY(grid, i, j, -1), W.B.y);
            approx(magneticX(grid, i, j, 5), W.B.x);
            approx(magneticY(grid, i, j, 5), W.B.y);
#endif
        }
    }
    //Corner
    fill_3D(grid);
    Boundary::Fixed(W,X | Y).apply(grid);
    expect_close(grid[-1,-1,-1], W);
}
