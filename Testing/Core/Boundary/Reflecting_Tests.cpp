
//
//  Reflecting_Tests.cpp
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

void DRAGON_Test::verify_boundary_reflective(bool output){
    if(output) std::cout<<"- 1D: ";
    verify_boundary_reflective_1D();
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- 2D: ";
    verify_boundary_reflective_2D();
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- 3D: ";
    verify_boundary_reflective_3D();
    if(output) std::cout<<"Passed\n";
}


//MARK: Helpers
static vec3 make_tagged_vec(double tag){ return vec3(tag + 100, tag + 200, tag + 300); }

static PrimitiveState G = make_tagged_state(-666);

void fill_1D(Grid1D& grid);
void fill_2D(Grid2D& grid);
void fill_3D(Grid3D& grid);

//MARK: Reflecting - 1D

void reflectX1D(PrimitiveState& W){
    W.v.x *= -1;
#ifdef MHD
    W.B *= -1;
    W.B.x *= -1;
#endif
}
void DRAGON_Test::verify_boundary_reflective_1D() {

    Grid1D grid(4, 1.0, 2);
    fill_1D(grid);

    Reflective(X).apply(grid);
    //Left
    PrimitiveState w = grid[0]; reflectX1D(w);
    expect_close(grid[-1], w);
    w = grid[1]; reflectX1D(w);
    expect_close(grid[-2], w);
    //Right
    w = grid[3]; reflectX1D(w);
    expect_close(grid[4], w);
    w = grid[2]; reflectX1D(w);
    expect_close(grid[5], w);
    //Single face filling: Left
    fill_1D(grid);
    Reflective("X-").apply(grid);
    w = grid[0]; reflectX1D(w);
    expect_close(grid[-1], w);
    w = grid[1]; reflectX1D(w);
    expect_close(grid[-2], w);
    //Single face filling: Right
    fill_1D(grid);
    Reflective("X+").apply(grid);
    w = grid[3]; reflectX1D(w);
    expect_close(grid[4], w);
    w = grid[2]; reflectX1D(w);
    expect_close(grid[5], w);
}
//MARK: Reflecting - 2D
void reflectX2D(PrimitiveState& W){
    W.v.x *= -1;
#ifdef MHD
    W.B.z *= -1;
#endif
}

void reflectY2D(PrimitiveState& W){
    W.v.y *= -1;
#ifdef MHD
    W.B.z *= -1;
#endif
}

void DRAGON_Test::verify_boundary_reflective_2D(){
    Grid2D grid(3, 4, 1.0, 1.0, 1);
    //X
    fill_2D(grid);
    Reflective(X,false).apply(grid);
    for (int j = 0; j < grid.getSizeY(); j++) {
        PrimitiveState w = grid[0,j]; reflectX2D(w);
        expect_close(grid[-1, j], w);
#ifdef MHD
        expect_close(grid.getA()[-1, j], grid.getA()[1, j]);
#endif
       
        w = grid[2,j]; reflectX2D(w);
        expect_close(grid[3, j], w);
#ifdef MHD
        expect_close(grid.getA()[4, j], grid.getA()[2, j]);
#endif
    }
    //No corners = no corners
    expect_close(grid[-1,-1],G);
    //Y
    fill_2D(grid);
    Reflective(Y).apply(grid);
    for (int i = 0; i < grid.getSizeX(); i++) {
        PrimitiveState w = grid[i,0]; reflectY2D(w);
        expect_close(grid[i,-1], w);
#ifdef MHD
        expect_close(grid.getA()[i, -1], grid.getA()[i, 1]);
#endif
       
        w = grid[i,3]; reflectY2D(w);
        expect_close(grid[i,4], w);
#ifdef MHD
        expect_close(grid.getA()[i, 5], grid.getA()[i, 3]);
#endif
    }
    //Corner
    fill_2D(grid);
    Reflective(X | Y).apply(grid);
    auto w = grid[0,0]; reflectX2D(w); reflectY2D(w);
    expect_close(grid[-1,-1], w);
}
//MARK: Reflecting - 3D

void reflectX3D(PrimitiveState& W){
    W.v.x *= -1;
}

void reflectY3D(PrimitiveState& W){
    W.v.y *= -1;
}

void reflectZ3D(PrimitiveState& W){
    W.v.z *= -1;
}

void DRAGON_Test::verify_boundary_reflective_3D(){
    Grid3D grid(3, 4,5,1.0, 1.0, 1.0, 1);
    //X
    fill_3D(grid);
    Reflective(X).apply(grid);
    for (int j = 0; j < grid.getSizeY(); j++) {
        for (int k = 0; k < grid.getSizeZ(); k++) {
            PrimitiveState w = grid[0,j,k]; reflectX3D(w);
            expect_close(grid[-1, j,k], w);
#ifdef MHD
            
#endif
            
            w = grid[2,j,k]; reflectX3D(w);
            expect_close(grid[3, j,k], w);
#ifdef MHD
            
#endif
        }
    }
    //Y
    fill_3D(grid);
    Reflective(Y).apply(grid);
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int k = 0; k < grid.getSizeZ(); k++) {
            PrimitiveState w = grid[i,0,k]; reflectY3D(w);
            expect_close(grid[i,-1,k], w);
#ifdef MHD
           
#endif
            
            w = grid[i,3,k]; reflectY3D(w);
            expect_close(grid[i,4,k], w);
#ifdef MHD
           
#endif
        }
    }
    //Z
    fill_3D(grid);
    Reflective("Z",false).apply(grid);
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int j = 0; j < grid.getSizeY(); j++) {
            PrimitiveState w = grid[i,j,0]; reflectZ3D(w);
            expect_close(grid[i,j,-1], w);
#ifdef MHD
            
#endif
            
            w = grid[i,j,4]; reflectZ3D(w);
            expect_close(grid[i,j,5], w);
#ifdef MHD
           
#endif
        }
    }
    //No corners = no corners
    expect_close(grid[1,-1,-1],G);
    //Corner
    fill_3D(grid);
    Reflective().apply(grid);
    auto w = grid[0,0,0]; reflectX3D(w); reflectY3D(w); reflectZ3D(w);
    expect_close(grid[-1,-1,-1], w);
}


