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
#include <fstream>
#include <iostream>

using namespace DRAGON_Test;
using namespace Boundary;

#ifdef MHD
namespace {
vec3 make_tagged_A(double tag){
    return {tag + 100, tag + 200, tag + 300};
}
}
#endif

constexpr std::string filename = "__TESTING__";

void DRAGON_Test::verify_IO(bool output){
    if(output) std::cout << "File I/O:\n";
    if(output) std::cout << "- 1D Restart: ";
    verify_IO1D();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "- 2D Restart: ";
    verify_IO2D();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "- 3D Restart: ";
    verify_IO3D();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "- Dimension Checks: ";
    verify_IO_dim_assert();
    verify_IO_size_assert();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "All File I/O Tests Passed.\n\n";

    std::remove((CONFIG::output_dir + "/" + filename + ".h5").c_str());//Clean up
}

void DRAGON_Test::verify_IO1D(){
    Grid1D g(5, 0.1, 2);
    for(int i = -2; i < 5 + 2; i++) {
        g[i] = make_tagged_state(i*0.1);
    }

    IO::writeToFile(g, 0.666, 666, filename);
    //Write doesn't change grid
    for(int i = -2; i < 5 + 2; i++) {
        expect_close(g[i],  make_tagged_state(i*0.1));
    }
    
    //Read
    Grid1D g2(5, 0.1, 2);

    double t; int n;
    IO::loadFromFile(g2, t, n, filename);
    
    assert(t==0.666);
    assert(n==666);
    assert(g2.dx==g.dx);

    
    #ifdef WRITE_GHOSTS_TO_FILE
    const int ng = g.getGhosts();
    #else
    const int ng = 0;
    #endif
    
    for(int i = -ng; i < 5+ng; i++) {
        expect_close(g[i],  g2[i]);
    }
    
}

void DRAGON_Test::verify_IO2D(){
    Grid2D g(3,4, 0.1,0.2, 2);
    
    for(int i = -2; i < 3 + 2; i++) {
        for(int j = -2; j < 4 + 2; j++) {
            g[i,j] = make_tagged_state(i*0.1 + j*0.01);
        }
    }
    #ifdef MHD
    for(int i = -2; i <= 3 + 2; i++) {
        for(int j = -2; j <= 4 + 2; j++) {
            g._A()[i,j] = make_tagged_A(i*0.1 + j*0.01);
        }
    }
    #endif

    IO::writeToFile(g, 0.666, 666, filename);
    //Write doesn't change grid
    for(int i = -2; i < 3 + 2; i++) {
        for(int j = -2; j < 4 + 2; j++) {
            expect_close(g[i,j],  make_tagged_state(i*0.1 + j*0.01));
            
        }
    }
    #ifdef MHD
    for(int i = -2; i <= 3 + 2; i++) {
        for(int j = -2; j <= 4 + 2; j++) {
            expect_close(g._A()[i,j], make_tagged_A(i*0.1 + j*0.01));
        }
    }
    #endif
    
    //Read
    Grid2D g2(3,4,0.0, 0.0, 2);

    double t; int n;
    IO::loadFromFile(g2, t, n, filename + ".h5");
    

    assert(t==0.666);
    assert(n==666);
    assert(g2.dx==g.dx);
    assert(g2.dy==g.dy);
    
    #ifdef WRITE_GHOSTS_TO_FILE
    const int ng = g.getGhosts();
    #else
    const int ng = 0;
    #endif
    
    for(int i = -ng; i < 3 + ng; i++) {
        for(int j = -ng; j < 4 + ng; j++) {
            expect_close(g[i,j],  g2[i,j]);
        }
    }
    #ifdef MHD
    for(int i = -ng; i <= 3 + ng; i++) {
        for(int j = -ng; j <= 4 + ng; j++) {
            expect_close(g._A()[i,j].z,  g2._A()[i,j].z);
        }
    }
    #endif
}
void DRAGON_Test::verify_IO3D(){
    Grid3D g(3,4,5, 0.1,0.2,0.3, 2);
    for(int i = -2; i < 3 + 2; i++) {
        for(int j = -2; j < 4 + 2; j++) {
            for(int k = -2; k < 5 + 2; k++) {
                g[i,j,k] = make_tagged_state(i*0.1 + j*0.01 + k*0.001);
            }
        }
    }
    #ifdef MHD
    for(int i = -2; i <= 3 + 2; i++) {
        for(int j = -2; j <= 4 + 2; j++) {
            for(int k = -2; k <= 5 + 2; k++) {
                g._A()[i,j,k] = make_tagged_A(i*0.1 + j*0.01 + k*0.001);
            }
        }
    }
    #endif

    IO::writeToFile(g, 0.666, 666, filename);
    //Write doesn't change grid
    for(int i = -2; i < 3 + 2; i++) {
        for(int j = -2; j < 4 + 2; j++) {
            for(int k = -2; k < 5 + 2; k++) {
                expect_close(g[i,j,k],  make_tagged_state(i*0.1 + j*0.01 + k*0.001));
            }
        }
    }
    #ifdef MHD
    for(int i = -2; i <= 3 + 2; i++) {
        for(int j = -2; j <= 4 + 2; j++) {
            for(int k = -2; k <= 5 + 2; k++) {
                expect_close(g._A()[i,j,k], make_tagged_A(i*0.1 + j*0.01 + k*0.001));
            }
        }
    }
    #endif
    
    //Read
    Grid3D g2(3,4,5,0.0, 0.0,0.0, 2);

    double t; int n;
    IO::loadFromFile(g2, t, n, filename + ".h5");
    

    assert(t==0.666);
    assert(n==666);
    assert(g2.dx==g.dx);
    assert(g2.dy==g.dy);
    assert(g2.dz==g.dz);

    #ifdef WRITE_GHOSTS_TO_FILE
    const int ng = g.getGhosts();
    #else
    const int ng = 0;
    #endif
    
    for(int i = -ng; i < 3 + ng; i++) {
        for(int j = -ng; j < 4 + ng; j++) {
            for(int k = -ng; k < 5 + ng; k++) {
                expect_close(g[i,j,k],  g2[i,j,k]);
            }
        }
    }
    #ifdef MHD
    for(int i = -ng; i <= 3 + ng; i++) {
        for(int j = -ng; j <= 4 + ng; j++) {
            for(int k = -ng; k <= 5 + ng; k++) {
                expect_close(g._A()[i,j,k],  g2._A()[i,j,k]);
            }
        }
    }
    #endif
}
void DRAGON_Test::verify_IO_dim_assert(){
    Grid2D g(3,4, 0.1,0.2, 2);
    
    for(int i = -2; i < 3 + 2; i++) {
        for(int j = -2; j < 4 + 2; j++) {
            g[i,j] = make_tagged_state(i*0.1 + j*0.01);
        }
    }

    IO::writeToFile(g, 0.666, 666, filename);
   
    
    try {
        Grid3D g2(3,4,5,0.0, 0.0,0.0, 2);

        double t; int n;
        IO::loadFromFile(g2, t, n, filename + ".h5");
    } catch (...) { return; } //Test passes on caught exception
    
    assert(false); //If exception isn't raised, test fails
}
void DRAGON_Test::verify_IO_size_assert(){
    Grid3D g(3,4,5,0.1, 0.1,0.2, 2);

    for(int i = -2; i < 3 + 2; i++) {
        for(int j = -2; j < 4 + 2; j++) {
            for(int k = -2; k < 5 + 2; k++) {
                g[i,j,k] = make_tagged_state(i*0.1 + j*0.01 + k*0.001);
            }
        }
    }

    IO::writeToFile(g, 0.666, 666, filename);


    try {
        Grid3D g2(3,5,4,0.0, 0.0,0.0, 2);

        double t; int n;
        IO::loadFromFile(g2, t, n, filename + ".h5");
    } catch (...) { return; } //Test passes on caught exception

    assert(false); //If exception isn't raised, test fails
}
