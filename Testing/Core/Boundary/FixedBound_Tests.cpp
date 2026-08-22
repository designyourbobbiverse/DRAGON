
//
//  Outflow_Tests.cpp
//  DRAGON/Testing/Core/Boundary
//
//  Created by Bobbie Markwick on 20/06/2026.
//

#include "Testing.hpp"
#include "Boundary/Boundary.hpp"

#include "Hydro/Grid.hpp"
#include <iostream> //For std::cout

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
        approx(grid._B()[-1,j].y, W.B.y);
        approx(grid._B()[3,j].y, W.B.y);
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
        approx(grid._B()[i,-1].x, W.B.x);
        approx(grid._B()[i,4].x, W.B.x);
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
            approx(grid._B()[-1,j,k].y, W.B.y);
            approx(grid._B()[-1,j,k].z, W.B.z);
            approx(grid._B()[3,j,k].y, W.B.y);
            approx(grid._B()[3,j,k].z, W.B.z);
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
                approx(grid._B()[i,-1,k].x, W.B.x);
                approx(grid._B()[i,4,k].x, W.B.x);
            }
            if (i + 1 < grid.getSizeX()) {
                approx(grid._B()[i,-1,k].z, W.B.z);
                approx(grid._B()[i,4,k].z, W.B.z);
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
            approx(grid._B()[i,j,-1].x, W.B.x);
            approx(grid._B()[i,j,-1].y, W.B.y);
            approx(grid._B()[i,j,5].x, W.B.x);
            approx(grid._B()[i,j,5].y, W.B.y);
            #endif
        }
    }
    //Corner
    fill_3D(grid);
    Boundary::Fixed(W,X | Y | Z).apply(grid);
    expect_close(grid[-1,-1,-1], W);
}
