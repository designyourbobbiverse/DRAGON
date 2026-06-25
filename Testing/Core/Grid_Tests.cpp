//
//  Grid_Tests.cpp
//  DRAGON/Testing/Core
//
//  Created by Bobbie Markwick on 09/06/2026.
//

#include "Testing.hpp"
#include "Boundary.hpp"
#include "Grid.hpp"
#include "Config.h"
#include <iostream>

using namespace DRAGON_Test;
using namespace Boundary;

void DRAGON_Test::verify_grid(bool output){
    if(output) std::cout << "Grid Structure:\n";
    if(output) std::cout << "- 1D: ";
    verify_grid1D();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "- 2D: ";
    verify_grid2D();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "- 3D: ";
    verify_grid3D();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "All Grid Tests Passed.\n\n";
}

PrimitiveState DRAGON_Test::make_tagged_state(double tag){
    PrimitiveState W = make_state(tag, tag + 10, tag + 20, tag + 30, tag + 40);
#ifdef MHD
    W.B = {tag + 50, tag + 60, tag + 70};
#endif
    return W;
}

void DRAGON_Test::verify_grid1D(){
    Grid1D g(5, 0.1, 2);
    assert(g.getSize() == 5);
    assert(g.getGhosts() == 2);
    for(int i = -2; i < 5 + 2; i++) {
        g[i] = make_tagged_state(i*0.1);
    }

    for(int i = -2; i < 5 + 2; i++) {
        expect_close(g[i],  make_tagged_state(i*0.1));
    }
}

void DRAGON_Test::verify_grid2D(){
    Grid2D g(3,4, 0.1,0.1, 2);
    assert(g.getSizeX() == 3);
    assert(g.getSizeY() == 4);
    assert(g.getGhosts() == 2);
    for(int i = -2; i < 3 + 2; i++) {
        for(int j = -2; j < 4 + 2; j++) {
            g[i,j] = make_tagged_state(i*0.1 + j*0.01);
        }
    }

    for(int i = -2; i < 3 + 2; i++) {
        for(int j = -2; j < 4 + 2; j++) {
            expect_close(g[i,j],  make_tagged_state(i*0.1 + j*0.01));
        }
    }
}
void DRAGON_Test::verify_grid3D(){
    Grid3D g(3,4,5, 0.1,0.1,0.1, 2);
    assert(g.getSizeX() == 3);
    assert(g.getSizeY() == 4);
    assert(g.getSizeZ() == 5);
    assert(g.getGhosts() == 2);
    for(int i = -2; i < 3 + 2; i++) {
        for(int j = -2; j < 4 + 2; j++) {
            for(int k = -2; k < 5 + 2; k++) {
                g[i,j,k] = make_tagged_state(i*0.1 + j*0.01 + k*0.001);
            }
        }
    }

    for(int i = -2; i < 3 + 2; i++) {
        for(int j = -2; j < 4 + 2; j++) {
            for(int k = -2; k < 5 + 2; k++) {
                expect_close(g[i,j,k],  make_tagged_state(i*0.1 + j*0.01 + k*0.001));
            }
        }
    }
}
//MARK: MHD Helpers
#ifdef MHD
double MagneticGrid::magneticX(Grid2D& grid, int i, int j) {
    return (grid.getA()[i, j + 1].z - grid.getA()[i, j].z) / grid.dy;
}

double MagneticGrid::magneticY(Grid2D& grid, int i, int j) {
    return (grid.getA()[i, j].z - grid.getA()[i + 1, j].z) / grid.dx;
}

double MagneticGrid::magneticX(Grid3D& grid, int i, int j, int k) {
    return (grid.getA()[i, j + 1, k].z - grid.getA()[i, j, k].z) / grid.dy
    - (grid.getA()[i, j, k + 1].y - grid.getA()[i, j, k].y) / grid.dz;
}

double MagneticGrid::magneticY(Grid3D& grid, int i, int j, int k) {
    return (grid.getA()[i, j, k + 1].x - grid.getA()[i, j, k].x) / grid.dz
    - (grid.getA()[i + 1, j, k].z - grid.getA()[i, j, k].z) / grid.dx;
}

double MagneticGrid::magneticZ(Grid3D& grid, int i, int j, int k) {
    return (grid.getA()[i + 1, j, k].y - grid.getA()[i, j, k].y) / grid.dx
    - (grid.getA()[i, j + 1, k].x - grid.getA()[i, j, k].x) / grid.dy;
}
#endif
