//
//  HydroGodunov_Tests.cpp
//  DRAGON/Testing/Core/Godunov
//
//  Created by Bobbie Markwick on 26/07/2026.
//

#include "Testing.hpp"
#include "Hydro/Grid.hpp"

#include <iostream> //For std::cout

using namespace DRAGON_Test;
using namespace Boundary;

//MARK: Uniform Flow Tests
void DRAGON_Test::verify_god_uniform_stationary_1D(){
    Grid1D grid(10, 1.0, 2);
    PrimitiveState W = make_state(1.0, 0.0, 0.0, 0.0, 5.0);
    for (int i = 0; i < grid.getSize(); i++) grid[i] = W;
    grid.boundary = Reflective();
    grid.advance(1.0);
    for (int i = 0; i < grid.getSize(); i++) expect_close(grid[i], W);
}
void DRAGON_Test::verify_god_uniform_moving_1D(){
    Grid1D grid(10, 1.0, 2);
    PrimitiveState W = make_state(1.0, 1.0, 2.0, 3.0, 5.0);
    for (int i = 0; i < grid.getSize(); i++) grid[i] = W;
    grid.boundary = Outflow();
    grid.advance(1.0);
    for (int i = 0; i < grid.getSize(); i++) expect_close(grid[i], W);
}

void DRAGON_Test::verify_god_uniform_stationary_2D(bool split){
    Grid2D grid(10,10,1.0, 1.0, 2);
    PrimitiveState W = make_state(1.0, 0.0, 0.0, 0.0, 5.0);
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int j = 0; j < grid.getSizeY(); j++) {
            grid[i,j] = W;
        }
    }
    grid.boundary = Reflective();
    
    if (split) grid.advance_split(1.0);
    else grid.advance_unsplit(1.0);
    
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int j = 0; j < grid.getSizeY(); j++) {
            expect_close(grid[i,j], W);
        }
    }
}
void DRAGON_Test::verify_god_uniform_moving_2D(bool split){
    Grid2D grid(10,10,1.0, 1.0, 2);
    PrimitiveState W = make_state(1.0, 1.0, 2.0, 3.0, 5.0);
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int j = 0; j < grid.getSizeY(); j++) {
            grid[i,j] = W;
        }
    }
    grid.boundary = Outflow();
    
    if (split) grid.advance_split(1.0);
    else grid.advance_unsplit(1.0);
    
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int j = 0; j < grid.getSizeY(); j++) {
            expect_close(grid[i,j], W);
        }
    }
}

void DRAGON_Test::verify_god_uniform_stationary_3D(bool split){
    Grid3D grid(10,10,10, 1.0,1.0,1.0,2);
    PrimitiveState W = make_state(1.0, 0.0, 0.0, 0.0, 5.0);
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int j = 0; j < grid.getSizeY(); j++) {
            for (int k=0; k<grid.getSizeZ(); k++) {
                grid[i,j,k] = W;
            }
        }
    }
    grid.boundary = Reflective();
    
    if (split) grid.advance_split(1.0);
    else grid.advance_unsplit(1.0);
    
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int j = 0; j < grid.getSizeY(); j++) {
            for (int k=0; k<grid.getSizeZ(); k++) {
                expect_close(grid[i,j,k], W);
            }
        }
    }
}
void DRAGON_Test::verify_god_uniform_moving_3D(bool split){
    Grid3D grid(10,10,10, 1.0,1.0, 1.0, 2);
    PrimitiveState W = make_state(1.0, 1.0, 2.0, 3.0, 5.0);
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int j = 0; j < grid.getSizeY(); j++) {
            for (int k = 0; k < grid.getSizeZ(); k++) {
                grid[i,j,k] = W;
            }
        }
    }
    grid.boundary = Outflow();
    
    if (split) grid.advance_split(1.0);
    else grid.advance_unsplit(1.0);
    
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int j = 0; j < grid.getSizeY(); j++) {
            for (int k = 0; k < grid.getSizeZ(); k++) {
                expect_close(grid[i,j,k], W);
            }
        }
    }
}

//MARK: Periodic Conservation Tests
void DRAGON_Test::verify_god_periodic_conservation_1D(){
    Grid1D grid(10, 1.0, 2);
    for (int i = 0; i < grid.getSize(); i++)
        grid[i] = make_state(1.0+0.1*i, 1.0+0.1*i, -0.1*i, 0.1*i, 10.0-0.1*i);
    grid.boundary = Periodic();
    
    ConservativeState expected = ConservativeState();
    for (int i = 0; i < grid.getSize(); i++) expected += ConservativeState(grid[i]);
    
    grid.advance(10.0);
    
    ConservativeState got = ConservativeState();
    for (int i = 0; i < grid.getSize(); i++) got += ConservativeState(grid[i]);
    
    expect_close(expected, got, 1e-14);

}

void DRAGON_Test::verify_god_periodic_conservation_2D(bool split){
    Grid2D grid(10,10,1.0, 1.0, 2);
    grid.boundary = Periodic();
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int j = 0; j < grid.getSizeY(); j++) {
            grid[i,j] = make_state(1.0+0.1*i+0.1*j, 1.0+0.1*i, 1.0-0.1*j, 0.1*i*j, 10.0-0.1*i+0.1*j);
        }
    }
    
    ConservativeState expected = ConservativeState();
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int j = 0; j < grid.getSizeY(); j++) {
            expected += ConservativeState(grid[i,j]);
        }
    }
    
    if (split) grid.advance_split(10.0);
    else grid.advance_unsplit(10.0);
    
    ConservativeState got = ConservativeState();
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int j = 0; j < grid.getSizeY(); j++) {
            got += ConservativeState(grid[i,j]);
            assert((grid[i,j].rho > 0));
        }
    }
    expect_close(expected, got, 1e-14);
}
void DRAGON_Test::verify_god_periodic_conservation_3D(bool split){
    Grid3D grid(10,10,10,1.0, 1.0,1.0, 2);
    grid.boundary = Periodic();
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int j = 0; j < grid.getSizeY(); j++) {
            for (int k=0; k<grid.getSizeZ(); k++) {
                grid[i,j,k] = make_state(1.0+0.1*i+0.1*j, 1.0+0.1*i, 1.0-0.1*j, 0.1*k, 10.0-0.1*i+0.1*j-0.1*k);
            }
        }
    }
    
    ConservativeState expected = ConservativeState();
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int j = 0; j < grid.getSizeY(); j++) {
            for (int k=0; k<grid.getSizeZ(); k++) {
                expected += ConservativeState(grid[i,j,k]);
            }
        }
    }
    
    if (split) grid.advance_split(2.0);
    else grid.advance_unsplit(2.0);
    
    ConservativeState got = ConservativeState();
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int j = 0; j < grid.getSizeY(); j++) {
            for (int k=0; k<grid.getSizeZ(); k++) {
                got += ConservativeState(grid[i,j,k]);
                assert((grid[i,j,k].rho > 0));
            }
        }
    }
    expect_close(expected, got, 1e-13);
}
//MARK: dt=0 Tests
void DRAGON_Test::verify_god_dt0_1D(){
    Grid1D grid(10, 1.0, 2), expected(10, 1.0, 2);
    grid.boundary = Periodic();

    for (int i = 0; i < grid.getSize(); i++) {
        grid[i] = make_state(1.0+0.1*i, 1.0+0.1*i, -0.1*i, 0.1*i, 10.0-0.1*i);
        expected[i] = grid[i];
    }
        
    grid.advance(0.0);
        
    for (int i = 0; i < grid.getSize(); i++) {
        expect_close(expected[i],grid[i]);
    }

}

void DRAGON_Test::verify_god_dt0_2D(bool split){
    Grid2D grid(10,10,1.0, 1.0, 2), expected(10,10,1.0, 1.0, 2);
    grid.boundary = Periodic();
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int j = 0; j < grid.getSizeY(); j++) {
            grid[i,j] = make_state(1.0+0.1*i+0.1*j, 1.0+0.1*i, 1.0-0.1*j, 0.1*i*j, 10.0-0.1*i+0.1*j);
            expected[i,j] = grid[i,j];
        }
    }

    if (split) grid.advance_split(0.0);
    else grid.advance_unsplit(0.0);
    
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int j = 0; j < grid.getSizeY(); j++) {
            expect_close(expected[i,j],grid[i,j]);
        }
    }
}
void DRAGON_Test::verify_god_dt0_3D(bool split){
    Grid3D grid(10,10,10,1.0, 1.0,1.0, 2), expected(10,10,10,1.0, 1.0,1.0, 2);
    grid.boundary = Periodic();
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int j = 0; j < grid.getSizeY(); j++) {
            for (int k=0; k<grid.getSizeZ(); k++) {
                grid[i,j,k] = make_state(1.0+0.1*i+0.1*j, 1.0+0.1*i, 1.0-0.1*j, 0.1*k, 10.0-0.1*i+0.1*j-0.1*k);
                expected[i,j,k] = grid[i,j,k];
            }
        }
    }
    if (split) grid.advance_split(0.0);
    else grid.advance_unsplit(0.0);
    
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int j = 0; j < grid.getSizeY(); j++) {
            for (int k=0; k<grid.getSizeZ(); k++) {
                expect_close(expected[i,j,k],grid[i,j,k]);
            }
        }
    }
}
