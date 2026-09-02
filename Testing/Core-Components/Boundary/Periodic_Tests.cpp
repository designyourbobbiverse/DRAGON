
//
//  Periodic_Tests.cpp
//  DRAGON/Testing/Core-Components/Boundary
//
//  Created by Bobbie Markwick on 20/06/2026.
//

#include "Testing.hpp"
#include "Boundary/Boundary.hpp"

#include "Hydro/Grid.hpp"
#include <iostream> //For std::cout

using namespace DRAGON_Test;
using namespace Boundary;

//MARK: Helpers

static PrimitiveState G = make_tagged_state(-666);

void fill_1D(Grid1D& grid);
void fill_2D(Grid2D& grid);
void fill_3D(Grid3D& grid);

//MARK: Periodic
void DRAGON_Test::verify_boundary_periodic(bool output){
    if (output) std::cout<<"- 1D: ";
    verify_boundary_periodic_1D();
    if (output) std::cout<<"Passed\n";
    if (output) std::cout<<"- 2D: ";
    verify_boundary_periodic_2D();
    if (output) std::cout<<"Passed\n";
    if (output) std::cout<<"- 3D: ";
    verify_boundary_periodic_3D();
    if (output) std::cout<<"Passed\n";
}
void DRAGON_Test::verify_boundary_periodic_1D(){
    Grid1D grid(4, 1.0, 2);
    fill_1D(grid);

    Periodic(X).apply(grid);
    //Left
    expect_close(grid[-1], grid[3]);
    expect_close(grid[-2], grid[2]);
    auto p1 = grid.passives()[-1,"Test"], p2 = grid.passives()[3,"Test"];
    assert(p1 == p2);
    //Right
    expect_close(grid[4], grid[0]);
    expect_close(grid[5], grid[1]);
    auto p3 = grid.passives()[4,"Test"], p4 = grid.passives()[0,"Test"];
    assert(p3 == p4);
}
//MARK: Periodic - 2D
void DRAGON_Test::verify_boundary_periodic_2D(){
    Grid2D grid(3, 4, 1.0, 1.0, 2);
    
    //X
    fill_2D(grid);
    Periodic(X).apply(grid);
    for (int j = 0; j < grid.getSizeY(); j++) {
        expect_close(grid[-1, j], grid[2, j]);
        expect_close(grid[3, j],  grid[0, j]);
        #ifdef MHD
        expect_close(grid._B()[-1,j], grid._B()[2,j]);
        expect_close(grid._B()[3,j], grid._B()[0,j]);
        expect_close(grid._B()[4,j], grid._B()[1,j]);
        #endif
        auto p1 = grid.passives()[-1,j,"Test"], p2 = grid.passives()[2,j,"Test"];
        assert(p1 == p2);
        auto p3 = grid.passives()[0,j,"Test"], p4 = grid.passives()[3,j,"Test"];
        assert(p3 == p4);
    }
    //Y
    fill_2D(grid);
    Periodic(Y).apply(grid);
    for (int i = 0; i < grid.getSizeX(); i++) {
        expect_close(grid[i,-1], grid[i,3]);
        expect_close(grid[i,4],  grid[i,0]);
        #ifdef MHD
        expect_close(grid._B()[i,-1], grid._B()[i,3]);
        expect_close(grid._B()[i,4], grid._B()[i,0]);
        expect_close(grid._B()[i,5], grid._B()[i,1]);
        #endif
        auto p1 = grid.passives()[i,-1,"Test"], p2 = grid.passives()[i,3,"Test"];
        assert(p1 == p2);
        auto p3 = grid.passives()[i,0,"Test"], p4 = grid.passives()[i,4,"Test"];
        assert(p3 == p4);
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
void DRAGON_Test::verify_boundary_periodic_3D(){
    Grid3D grid(3, 4, 5, 1.0, 1.0, 1.0, 1);
    
    //X
    fill_3D(grid);
    Periodic("X").apply(grid);
    for (int j = 0; j < grid.getSizeY(); j++) {
        for (int k = 0; k < grid.getSizeZ(); k++) {
            expect_close(grid[-1, j,k], grid[2, j,k]);
            expect_close(grid[3, j,k],  grid[0, j,k]);
            #ifdef MHD
            expect_close(grid._B()[-1,j,k], grid._B()[2,j,k]);

            expect_close(grid._B()[3,j,k].x, grid._B()[0,j,k].x);
            expect_close(grid._B()[3,j,k].y, grid._B()[0,j,k].y);
            expect_close(grid._B()[3,j,k].z, grid._B()[0,j,k].z);

            expect_close(grid._B()[4,j,k].x, grid._B()[1,j,k].x);
            expect_close(grid._B()[4,j,k].y, grid._B()[1,j,k].y);
            expect_close(grid._B()[4,j,k].z, grid._B()[1,j,k].z);
            #endif
            auto p1 = grid.passives()[-1,j,k,"Test"], p2 = grid.passives()[2,j,k,"Test"];
            assert(p1 == p2);
            auto p3 = grid.passives()[0,j,k,"Test"], p4 = grid.passives()[3,j,k,"Test"];
            assert(p3 == p4);
        }
    }
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
            expect_close(grid._B()[i,-1,k].x, grid._B()[i,3,k].x);
            expect_close(grid._B()[i,-1,k].y, grid._B()[i,3,k].y);
            expect_close(grid._B()[i,-1,k].z, grid._B()[i,3,k].z);

            expect_close(grid._B()[i,4,k].x, grid._B()[i,0,k].x);
            expect_close(grid._B()[i,4,k].y, grid._B()[i,0,k].y);
            expect_close(grid._B()[i,4,k].z, grid._B()[i,0,k].z);

            expect_close(grid._B()[i,5,k].x, grid._B()[i,1,k].x);
            expect_close(grid._B()[i,5,k].y, grid._B()[i,1,k].y);
            expect_close(grid._B()[i,5,k].z, grid._B()[i,1,k].z);
            #endif
            auto p1 = grid.passives()[i,-1,k,"Test"], p2 = grid.passives()[i,3,k,"Test"];
            assert(p1 == p2);
            auto p3 = grid.passives()[i,0,k,"Test"], p4 = grid.passives()[i,4,k,"Test"];
            assert(p3 == p4);
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
            expect_close(grid._B()[i,j,-1].x, grid._B()[i,j,4].x);
            expect_close(grid._B()[i,j,-1].y, grid._B()[i,j,4].y);
            expect_close(grid._B()[i,j,-1].z, grid._B()[i,j,4].z);

            expect_close(grid._B()[i,j,5].x, grid._B()[i,j,0].x);
            expect_close(grid._B()[i,j,5].y, grid._B()[i,j,0].y);
            expect_close(grid._B()[i,j,5].z, grid._B()[i,j,0].z);

            expect_close(grid._B()[i,j,6].x, grid._B()[i,j,1].x);
            expect_close(grid._B()[i,j,6].y, grid._B()[i,j,1].y);
            expect_close(grid._B()[i,j,6].z, grid._B()[i,j,1].z);
            #endif
            auto p1 = grid.passives()[i,j,-1,"Test"], p2 = grid.passives()[i,j,4,"Test"];
            assert(p1 == p2);
            auto p3 = grid.passives()[i,j,0,"Test"], p4 = grid.passives()[i,j,5,"Test"];
            assert(p3 == p4);
        }
    }
    //Corner
    fill_3D(grid);
    Periodic(X | Y | Z).apply(grid);
    expect_close(grid[-1,-1, -1], grid[2,3,4]);
}
