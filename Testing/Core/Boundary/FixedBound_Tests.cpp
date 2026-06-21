
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

static void expect_fixed_transverse_B_X(Grid2D& grid, const PrimitiveState& state, int i, int j) {
    assert(approx(magneticY(grid, i, j), state.B.y));
}

static void expect_fixed_transverse_B_Y(Grid2D& grid, const PrimitiveState& state, int i, int j) {
    assert(approx(magneticX(grid, i, j), state.B.x));
}

static void expect_fixed_transverse_B_X(Grid3D& grid, const PrimitiveState& state, int i, int j, int k) {
    assert(approx(magneticY(grid, i, j, k), state.B.y));
    assert(approx(magneticZ(grid, i, j, k), state.B.z));
}

static void expect_fixed_transverse_B_Y(Grid3D& grid, const PrimitiveState& state, int i, int j, int k) {
    assert(approx(magneticX(grid, i, j, k), state.B.x));
    assert(approx(magneticZ(grid, i, j, k), state.B.z));
}

static void expect_fixed_Bx(Grid3D& grid, const PrimitiveState& state, int i, int j, int k) {
    assert(approx(magneticX(grid, i, j, k), state.B.x));
}

static void expect_fixed_Bz(Grid3D& grid, const PrimitiveState& state, int i, int j, int k) {
    assert(approx(magneticZ(grid, i, j, k), state.B.z));
}

static void expect_fixed_transverse_B_Z(Grid3D& grid, const PrimitiveState& state, int i, int j, int k) {
    assert(approx(magneticX(grid, i, j, k), state.B.x));
    assert(approx(magneticY(grid, i, j, k), state.B.y));
}
#endif

void fill_1D(Grid1D& grid);
void fill_2D(Grid2D& grid);
void fill_3D(Grid3D& grid);

//MARK: Fixed
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
        expect_fixed_transverse_B_X(grid, W, -1, j);
        expect_fixed_transverse_B_X(grid, W, 3, j);
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
        expect_fixed_transverse_B_Y(grid, W, i, -1);
        expect_fixed_transverse_B_Y(grid, W, i, 4);
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
            expect_fixed_transverse_B_X(grid, W, -1, j, k);
            expect_fixed_transverse_B_X(grid, W, 3, j, k);
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
                expect_fixed_Bx(grid, W, i, -1, k);
                expect_fixed_Bx(grid, W, i, 4, k);
            }
            if (i + 1 < grid.getSizeX()) {
                expect_fixed_Bz(grid, W, i, -1, k);
                expect_fixed_Bz(grid, W, i, 4, k);
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
            expect_fixed_transverse_B_Z(grid, W, i, j, -1);
            expect_fixed_transverse_B_Z(grid, W, i, j, 5);
#endif
        }
    }
    //Corner
    fill_3D(grid);
    Boundary::Fixed(W,X | Y).apply(grid);
    expect_close(grid[-1,-1,-1], W);
}
