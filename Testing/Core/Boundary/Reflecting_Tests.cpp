
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
#ifdef MHD
using namespace MagneticGrid;
#endif

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
#ifdef MHD
    if(output) std::cout<<"- 1D Conductive: ";
    verify_boundary_reflective_conductive_1D();
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- 2D Conductive: ";
    verify_boundary_reflective_conductive_2D();
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- 3D Conductive: ";
    verify_boundary_reflective_conductive_3D();
    if(output) std::cout<<"Passed\n";
#endif
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

void reflectX1DNonconductive(PrimitiveState& W){
    W.v.x *= -1;
}

void DRAGON_Test::verify_boundary_reflective_1D() {

    Grid1D grid(4, 1.0, 2);
    fill_1D(grid);

    #ifdef MHD
    Reflective(X, true, false).apply(grid);
    #else
    Reflective(X).apply(grid);
    #endif
    //Left
    PrimitiveState w = grid[0]; reflectX1DNonconductive(w);
    expect_close(grid[-1], w);
    w = grid[1]; reflectX1DNonconductive(w);
    expect_close(grid[-2], w);
    //Right
    w = grid[3]; reflectX1DNonconductive(w);
    expect_close(grid[4], w);
    w = grid[2]; reflectX1DNonconductive(w);
    expect_close(grid[5], w);
    //Single face filling: Left
    fill_1D(grid);
    #ifdef MHD
    Reflective("X-", true, false).apply(grid);
    #else
    Reflective("X-").apply(grid);
    #endif
    w = grid[0]; reflectX1DNonconductive(w);
    expect_close(grid[-1], w);
    w = grid[1]; reflectX1DNonconductive(w);
    expect_close(grid[-2], w);
    //Single face filling: Right
    fill_1D(grid);
    #ifdef MHD
    Reflective("X+", true, false).apply(grid);
    #else
    Reflective("X+").apply(grid);
    #endif
    w = grid[3]; reflectX1DNonconductive(w);
    expect_close(grid[4], w);
    w = grid[2]; reflectX1DNonconductive(w);
    expect_close(grid[5], w);
}
//MARK: Reflecting - 2D
void reflectX2D(PrimitiveState& W, bool conductive=true){
    W.v.x *= -1;
#ifdef MHD
    if(conductive) W.B.z *= -1;
#endif
}

void reflectY2D(PrimitiveState& W, bool conductive=true){
    W.v.y *= -1;
#ifdef MHD
    if(conductive) W.B.z *= -1;
#endif
}

void DRAGON_Test::verify_boundary_reflective_2D(){
    Grid2D grid(3, 4, 1.0, 1.0, 1);
    //X
    fill_2D(grid);
    #ifdef MHD
    Reflective(X, false, false).apply(grid);
    #else
    Reflective(X, false).apply(grid);
    #endif
    for (int j = 0; j < grid.getSizeY(); j++) {
        PrimitiveState w = grid[0,j]; reflectX2D(w,false);
        expect_close(grid[-1, j], w);
        #ifdef MHD //Transverse Magnetic Fields
        assert(magneticY(grid, -1, j) == magneticY(grid, 0, j));
        #endif
       
        w = grid[2,j]; reflectX2D(w,false);
        expect_close(grid[3, j], w);
        #ifdef MHD //Transverse Magnetic Fields
        assert(magneticY(grid, 3, j) == magneticY(grid, 2, j));
        #endif
    }
    //No corners = no corners
    expect_close(grid[-1,-1],G);
    //Y
    fill_2D(grid);
    #ifdef MHD
    Reflective(Y, true, false).apply(grid);
    #else
    Reflective(Y).apply(grid);
    #endif
    for (int i = 0; i < grid.getSizeX(); i++) {
        PrimitiveState w = grid[i,0]; reflectY2D(w,false);
        expect_close(grid[i,-1], w);
        #ifdef MHD //Transverse Magnetic Fields
        assert(magneticX(grid, i,-1) == magneticX(grid, i,0));
        #endif
       
        w = grid[i,3]; reflectY2D(w,false);
        expect_close(grid[i,4], w);
        #ifdef MHD //Transverse Magnetic Fields
        assert(magneticX(grid, i,4) == magneticX(grid, i,3));
        #endif
    }
    //Corner
    fill_2D(grid);
    #ifdef MHD
    Reflective(X | Y, true, false).apply(grid);
    #else
    Reflective(X | Y).apply(grid);
    #endif
    auto w = grid[0,0]; reflectX2D(w,false); reflectY2D(w,false);
    expect_close(grid[-1,-1], w);
    #ifdef MHD //A field corner checks
    int ng = grid.getGhosts(), nx = grid.getSizeX(), ny = grid.getSizeY();
    fill_2D(grid);
    Reflective(X, true, false).apply(grid);
    assert(approx(grid.getA()[nx+ng,ny+ng].z, 2*grid.getA()[nx+ng-1,ny+ng].z - grid.getA()[nx+ng-2,ny+ng].z));

    fill_2D(grid);
    Reflective(Y, true, false).apply(grid);
    assert(approx(grid.getA()[nx+ng,ny+ng].z, 2*grid.getA()[nx+ng,ny+ng-1].z - grid.getA()[nx+ng,ny+ng-2].z));
    #endif
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
    #ifdef MHD
    Reflective(X, true, false).apply(grid);
    #else
    Reflective(X).apply(grid);
    #endif
    for (int j = 0; j < grid.getSizeY(); j++) {
        for (int k = 0; k < grid.getSizeZ(); k++) {
            PrimitiveState w = grid[0,j,k]; reflectX3D(w);
            expect_close(grid[-1, j,k], w);
            #ifdef MHD //Transverse Magnetic Fields
            assert(magneticY(grid, -1, j,k) == magneticY(grid, 0, j,k));
            assert(magneticZ(grid, -1, j,k) == magneticZ(grid, 0, j,k));
            #endif
            
            w = grid[2,j,k]; reflectX3D(w);
            expect_close(grid[3, j,k], w);
            #ifdef MHD //Transverse Magnetiic Fields
            assert(magneticY(grid, 3, j,k) == magneticY(grid, 2, j,k));
            assert(magneticZ(grid, 3, j,k) == magneticZ(grid, 2, j,k));
            #endif
        }
    }
    //Y
    fill_3D(grid);
    #ifdef MHD
    Reflective(Y, true, false).apply(grid);
    #else
    Reflective(Y).apply(grid);
    #endif
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int k = 0; k < grid.getSizeZ(); k++) {
            PrimitiveState w = grid[i,0,k]; reflectY3D(w);
            expect_close(grid[i,-1,k], w);
            #ifdef MHD //Transverse Magnetic Fields
            assert(magneticX(grid, i,-1,k) == magneticX(grid, i,0,k));
            assert(magneticZ(grid, i,-1,k) == magneticZ(grid, i,0,k));
            #endif
            
            w = grid[i,3,k]; reflectY3D(w);
            expect_close(grid[i,4,k], w);
            #ifdef MHD //Transverse Magnetic Fields
            assert(magneticX(grid, i,4,k) == magneticX(grid, i,3,k));
            assert(magneticZ(grid, i,4,k) == magneticZ(grid, i,3,k));
            #endif
        }
    }
    //Z
    fill_3D(grid);
    #ifdef MHD
    Reflective("Z", true, false).apply(grid);
    #else
    Reflective("Z", false).apply(grid);
    #endif
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int j = 0; j < grid.getSizeY(); j++) {
            PrimitiveState w = grid[i,j,0]; reflectZ3D(w);
            expect_close(grid[i,j,-1], w);
            #ifdef MHD //Transverse Magnetic Fields
            assert(magneticX(grid, i,j,-1) == magneticX(grid, i,j,0));
            assert(magneticY(grid, i,j,-1) == magneticY(grid, i,j,0));
            #endif
            
            w = grid[i,j,4]; reflectZ3D(w);
            expect_close(grid[i,j,5], w);
            #ifdef MHD //Transverse Magnetic Fields
            assert(magneticX(grid, i,j,5) == magneticX(grid, i,j,4));
            assert(magneticY(grid, i,j,5) == magneticY(grid, i,j,4));
            #endif
        }
    }
    //No corners = no corners
    fill_3D(grid);
    Reflective("Z", false).apply(grid);
    expect_close(grid[1,-1,-1],G);
    //Corner
    fill_3D(grid);
    #ifdef MHD
    Reflective(X | Y | Z, true, false).apply(grid);
    #else
    Reflective().apply(grid);
    #endif
    auto w = grid[0,0,0]; reflectX3D(w); reflectY3D(w); reflectZ3D(w);
    expect_close(grid[-1,-1,-1], w);
    #ifdef MHD//A field corner checks
    int ng = grid.getGhosts(), nx = grid.getSizeX(), ny = grid.getSizeY(), nz = grid.getSizeZ();
    fill_3D(grid);
    Reflective(X, true, false).apply(grid);
    vec3 expected = 2*grid.getA()[nx+ng-1,ny+ng,nz+ng] - grid.getA()[nx+ng-2,ny+ng,nz+ng];
    expected.x = grid.getA()[nx+ng-1,ny+ng,nz+ng].x;
    expect_close(grid.getA()[nx+ng,ny+ng,nz+ng], expected);

    fill_3D(grid);
    Reflective(Y, true, false).apply(grid);
    expected = 2*grid.getA()[nx+ng,ny+ng-1,nz+ng] - grid.getA()[nx+ng,ny+ng-2,nz+ng];
    expected.y = grid.getA()[nx+ng,ny+ng-1,nz+ng].y;
    expect_close(grid.getA()[nx+ng,ny+ng,nz+ng], expected);

    fill_3D(grid);
    Reflective(Z, true, false).apply(grid);
    expected = 2*grid.getA()[nx+ng,ny+ng,nz+ng-1] - grid.getA()[nx+ng,ny+ng,nz+ng-2];
    expected.z = grid.getA()[nx+ng,ny+ng,nz+ng-1].z;
    expect_close(grid.getA()[nx+ng,ny+ng,nz+ng], expected);
    #endif
}
//MARK: Conductive - 1D
#ifdef MHD
void DRAGON_Test::verify_boundary_reflective_conductive_1D() {
    Grid1D grid(4, 1.0, 2);
    fill_1D(grid);

    Reflective(X).apply(grid);
    PrimitiveState w = grid[0]; reflectX1D(w);
    expect_close(grid[-1], w);
    w = grid[1]; reflectX1D(w);
    expect_close(grid[-2], w);
    w = grid[3]; reflectX1D(w);
    expect_close(grid[4], w);
    w = grid[2]; reflectX1D(w);
    expect_close(grid[5], w);
}
//MARK: Conductive - 2D

void DRAGON_Test::verify_boundary_reflective_conductive_2D(){
    Grid2D grid(3, 4, 1.0, 1.0, 1);

    fill_2D(grid);
    Reflective(X).apply(grid);
    for (int j = 0; j < grid.getSizeY(); j++) {
        PrimitiveState w = grid[0,j]; reflectX2D(w);
        expect_close(grid[-1, j], w);
        expect_close(grid.getA()[-1, j], grid.getA()[1, j]);

        w = grid[2,j]; reflectX2D(w);
        expect_close(grid[3, j], w);
        expect_close(grid.getA()[4, j], grid.getA()[2, j]);
    }

    fill_2D(grid);
    Reflective(Y).apply(grid);
    for (int i = 0; i < grid.getSizeX(); i++) {
        PrimitiveState w = grid[i,0]; reflectY2D(w);
        expect_close(grid[i,-1], w);
        expect_close(grid.getA()[i, -1], grid.getA()[i, 1]);

        w = grid[i,3]; reflectY2D(w);
        expect_close(grid[i,4], w);
        expect_close(grid.getA()[i, 5], grid.getA()[i, 3]);
    }
}


//MARK: Conductive - 3D
void DRAGON_Test::verify_boundary_reflective_conductive_3D(){
    Grid3D grid(3, 4, 5, 1.0, 1.0, 1.0, 1);

    fill_3D(grid);
    Reflective(X).apply(grid);
    for (int j = 0; j < grid.getSizeY(); j++) {
        for (int k = 0; k < grid.getSizeZ(); k++) {
            PrimitiveState w = grid[0,j,k]; reflectX3D(w);
            expect_close(grid[-1, j,k], w);
            vec3 A = grid.getA()[1, j, k];
            A.x = grid.getA()[0, j, k].x;
            expect_close(grid.getA()[-1, j, k], A);

            w = grid[2,j,k]; reflectX3D(w);
            expect_close(grid[3, j,k], w);
            A = grid.getA()[2, j, k];
            A.x = -grid.getA()[1, j, k].x;
            expect_close(grid.getA()[4, j, k], A);
            assert(approx(grid.getA()[3, j, k].x, -grid.getA()[2, j, k].x));
        }
    }

    fill_3D(grid);
    Reflective(Y).apply(grid);
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int k = 0; k < grid.getSizeZ(); k++) {
            PrimitiveState w = grid[i,0,k]; reflectY3D(w);
            expect_close(grid[i,-1,k], w);
            vec3 A = grid.getA()[i, 1, k];
            A.y = -grid.getA()[i, 0, k].y;
            expect_close(grid.getA()[i, -1, k], A);

            w = grid[i,3,k]; reflectY3D(w);
            expect_close(grid[i,4,k], w);
            A = grid.getA()[i, 3, k];
            A.y = -grid.getA()[i, 2, k].y;
            expect_close(grid.getA()[i, 5, k], A);
            assert(approx(grid.getA()[i, 4, k].y, -grid.getA()[i, 3, k].y));
        }
    }

    fill_3D(grid);
    Reflective("Z", false).apply(grid);
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int j = 0; j < grid.getSizeY(); j++) {
            PrimitiveState w = grid[i,j,0]; reflectZ3D(w);
            expect_close(grid[i,j,-1], w);
            vec3 A = grid.getA()[i, j, 1];
            A.z = -grid.getA()[i, j, 0].z;
            expect_close(grid.getA()[i, j, -1], A);

            w = grid[i,j,4]; reflectZ3D(w);
            expect_close(grid[i,j,5], w);
            A = grid.getA()[i, j, 4];
            A.z = -grid.getA()[i, j, 3].z;
            expect_close(grid.getA()[i, j, 6], A);
            assert(approx(grid.getA()[i, j, 5].z, -grid.getA()[i, j, 4].z));
        }
    }
}
#endif
