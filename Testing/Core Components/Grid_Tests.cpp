//
//  Grid_Tests.cpp
//  DRAGON
//
//  Created by Bobbie Markwick on 09/06/2026.
//

#include "Testing.hpp"
#include "Boundary.hpp"
#include "Grid.hpp"
#include <iostream>

using namespace DRAGON_Test;
using namespace Boundary;

void DRAGON_Test::verify_grid(){
    verify_grid1D();
    verify_grid2D();
    verify_grid3D();
}

PrimitiveState DRAGON_Test::make_tagged_state(double tag){
    return make_state(tag, tag + 10, tag + 20, tag + 30, tag + 40);
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
