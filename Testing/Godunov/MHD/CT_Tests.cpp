//
//  CT_Tests.cpp
//  DRAGON/Testing/Godunov/MHD
//
//  Created by Bobbie Markwick on 25/06/2026.
//

#include "Testing.hpp"
#include "MHD/CT.hpp"

#include "Refinement/DistGrid.hpp"
#include "Constants.h"      //For _pi, _1_8pi
#include <iostream>

#ifdef MHD
using namespace DRAGON_Test;



//MARK: Stationary Field
void DRAGON_Test::verify_ct_stationary_2D(){
    double dx = _pi/5;
    Grid2D grid(10,10,dx,dx, 2), expected(10,10,dx,dx, 2);
    double p0 = 5.0;
    PrimitiveState W = make_state(1.0, 0.0, 0.0, 0.0, p0);
    W.B = {0.2, -0.3, 0.4};

    for (int i = 0; i <= grid.getSizeX(); i++) {
        for (int j = 0; j <= grid.getSizeY(); j++) {
            grid[i,j] = W;
            grid._B()[i,j] = W.B;
        }
    }
    grid.boundary = Boundary::Periodic();
    grid.boundary.apply(grid);
    grid.initialize_B_fields();

    for (int i = 0; i <= grid.getSizeX(); i++) {
        for (int j = 0; j <= grid.getSizeY(); j++) {
            grid[i,j].p = p0 - grid[i,j].B * grid[i,j].B * _1_8pi;
            expected[i,j] = grid[i,j];
        }
    }
    assert_divergenceless(grid._B(),dx,dx);

    
    grid.advance(1.0);
    
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int j = 0; j < grid.getSizeY(); j++) {
            expect_close(grid[i,j].B, expected[i,j].B);
        }
    }
    assert_divergenceless(grid._B(),dx,dx);
}

void DRAGON_Test::verify_ct_stationary_3D(){
    double dx = _pi/5;
    Grid3D grid(10,10,10,dx,dx,dx, 2), expected(10,10,10,dx,dx,dx, 2);
    double p0 = 5.0;
    PrimitiveState W = make_state(1.0, 0.0, 0.0, 0.0, p0);
    W.B = {0.2, -0.3, 0.4};

    for (int i = 0; i <= grid.getSizeX(); i++) {
        for (int j = 0; j <= grid.getSizeY(); j++) {
            for (int k = 0; k <= grid.getSizeZ(); k++) {
                grid[i,j,k] = W;
                grid._B()[i,j,k] = W.B;
            }
        }
    }
    grid.boundary = Boundary::Periodic();
    grid.boundary.apply(grid);
    grid.initialize_B_fields();

    for (int i = 0; i <= grid.getSizeX(); i++) {
        for (int j = 0; j <= grid.getSizeY(); j++) {
            for (int k = 0; k <= grid.getSizeZ(); k++) {
                expected[i,j,k] = grid[i,j,k];
            }
        }
    }
    assert_divergenceless(grid._B(),dx,dx,dx);

    
    grid.advance(0.1);
    
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int j = 0; j < grid.getSizeY(); j++) {
            for (int k = 0; k < grid.getSizeZ(); k++) {
                expect_close(grid[i,j,k].B, expected[i,j,k].B);
            }
        }
    }
    assert_divergenceless(grid._B(),dx,dx,dx);
}


//MARK: d/dz = 0 -> Bz stays zero
void DRAGON_Test::verify_ct_uniform_Bz0_2D(){
    constexpr int n = 32;
    constexpr double dx = 2.0/n;
    
    auto grid = Grid2D(n,n, dx, dx);
    grid.boundary = Boundary::Periodic();
    
    for (int i=0; i<=grid.getSizeX();i++) {
        for (int j=0; j<=grid.getSizeY(); j++) {
            grid[i,j] = make_state(1.0, 2.0, 3.0, 4.0, 5.0);
            grid._B()[i,j] = {0.2,0.3,0.0};
        }
    }    
    grid.initialize_B_fields();
    
    grid.advance(1.0);
    for (int i=0; i<grid.getSizeX();i++) {
        for (int j=0; j<grid.getSizeY(); j++) {
            assert(approx(grid[i,j].B.z,0));
        }
    }
    for (int i=0; i<=grid.getSizeX();i++) {
        for (int j=0; j<=grid.getSizeY(); j++) {
            grid._B()[i,j] = {0.0,0.2,0.3};
        }
    }
    grid.initialize_B_fields();
    
    grid.advance(0.1);
    for (int i=0; i<grid.getSizeX();i++) {
        for (int j=0; j<grid.getSizeY(); j++) {
            assert(approx(grid[i,j].B.x,0));
        }
    }
    for (int i=0; i<=grid.getSizeX();i++) {
        for (int j=0; j<=grid.getSizeY(); j++) {
            grid._B()[i,j] = {0.2,0.0,0.3};
        }
    }
    grid.initialize_B_fields();
    
    grid.advance(0.1);
    for (int i=0; i<grid.getSizeX();i++) {
        for (int j=0; j<grid.getSizeY(); j++) {
            assert(approx(grid[i,j].B.y,0));
        }
    }
}


void DRAGON_Test::verify_ct_uniform_Bz0_3D(){
    constexpr int n = 16;
    constexpr double dx = 2.0/n;
    
    auto grid = Grid3D(n,n,n, dx, dx, dx);
    grid.boundary = Boundary::Periodic();
    
    for (int i=0; i<=grid.getSizeX();i++) {
        for (int j=0; j<=grid.getSizeY(); j++) {
            for (int k=0; k<=grid.getSizeZ(); k++) {
                grid[i,j,k] = make_state(1.0, 2.0, 3.0, 4.0, 5.0);
                grid._B()[i,j,k] = {0.2,0.3,0.0};
            }
        }
    }
    grid.initialize_B_fields();
    
    grid.advance(1.0);
    for (int i=0; i<grid.getSizeX();i++) {
        for (int j=0; j<grid.getSizeY(); j++) {
            for (int k=0; k<grid.getSizeZ(); k++) {
                assert(approx(grid[i,j,k].B.z,0));
            }
        }
    }
}
#endif
