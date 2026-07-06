//
//  Grid_Tests.cpp
//  DRAGON/Testing/Core
//
//  Created by Bobbie Markwick on 09/06/2026.
//

#include "Testing.hpp"
#include "Boundary.hpp"
#include "HDF5Output.hpp"
#include "Config.h"
#include <iostream>

using namespace DRAGON_Test;
using namespace Boundary;

void DRAGON_Test::verify_IO(bool output){
    if(output) std::cout << "File I/O:\n";
    if(output) std::cout << "- 1D: ";
    verify_IO1D();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "- 2D: ";
    verify_IO2D();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "- 3D: ";
    verify_IO3D();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "All Grid Tests Passed.\n\n";
}

void DRAGON_Test::verify_IO1D(){
    Grid1D g(5, 0.1, 2);
    assert(g.getSize() == 5);
    assert(g.getGhosts() == 2);
    for(int i = -2; i < 5 + 2; i++) {
        g[i] = make_tagged_state(i*0.1);
    }

    IO::writeToFile(g, 0.666, 666, "TEST");
    //Write doesn't change grid
    for(int i = 0; i < 5; i++) {
        expect_close(g[i],  make_tagged_state(i*0.1));
    }
    
    //Read
    Grid1D g2(5, 0.1, 2);

    double t; int n;
    IO::loadFromFile(g2, t, n, "TEST.h5");
    

    for(int i = 0; i < 5; i++) {
        expect_close(g[i],  g2[i]);
    }
    
    assert(t==0.666);
    assert(n==666);
    
}

void DRAGON_Test::verify_IO2D(){
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
void DRAGON_Test::verify_IO3D(){
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
