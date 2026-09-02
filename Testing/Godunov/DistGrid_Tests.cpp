//
//  Godunov_Tests.cpp
//  DRAGON/Testing/Godunov
//
//  Created by Bobbie Markwick on 09/06/2026.
//

#include "Testing.hpp"
#include "Refinement/DistGrid.hpp"

#include <iostream>       //For std::cout

using namespace DRAGON_Test;
using namespace Boundary;

void DRAGON_Test::verify_god_dist_grid_1D(){
    Grid1D grid(100, 1.0);
    DistGrid1D dgrid(100, 1.0);
    
    grid.boundary = Periodic();
    dgrid.boundary = Periodic();
    grid.passives().add("Test");
    dgrid.passives().add("Test");
    
    for (int i = 0; i < grid.getSize(); i++) {
        grid[i] = make_state(1.0+0.01*i, 1.0+0.01*i, -0.01*i, 0.01*i, 10.0-0.01*i);
        #ifdef MHD
        grid[i].B = vec3{0, -0.2*i, -0.3*i};
        #endif
        dgrid[i] = grid[i];
        
        grid.passives()[i,"Test"] = 10.0-0.1*i;
        dgrid.passives()[i,"Test"] = grid.passives()[i,"Test"];
    }
        
    
    grid.advance(1.0);
    dgrid.advance(1.0);
    
    for (int i = 0; i < grid.getSize(); i++) {
        expect_close(grid[i], dgrid[i], 1e-26);
        assert(approx(grid.passives()[i,"Test"],  dgrid.passives()[i,"Test"], 1e-26));
    }
}


void DRAGON_Test::verify_god_dist_grid_2D(){
    Grid2D grid(32,32,1.0, 1.0, 4);
    DistGrid2D dgrid(32,32,1.0, 1.0, 4);
    grid.boundary = Periodic();
    dgrid.boundary = Periodic();
    grid.passives().add("Test");
    dgrid.passives().add("Test");
    
    for (int i = 0; i <= grid.getSizeX(); i++) {
        for (int j = 0; j <= grid.getSizeY(); j++) {
            grid[i,j] = make_state(1.0+0.1*i+0.1*j, 1.0+0.1*i, 1.0-0.1*j, 0.1*i*j, 10.0-0.1*i+0.1*j);
            dgrid[i,j] = grid[i,j];
            
            grid.passives()[i,j,"Test"] = 10.0-0.1*i+0.1*j;
            dgrid.passives()[i,j,"Test"] = grid.passives()[i,j,"Test"];
        }
    }
    
    grid.advance(1.0);
    dgrid.advance(1.0);
    
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int j = 0; j < grid.getSizeY(); j++) {
            expect_close(grid[i,j], dgrid[i,j], 1e-26);
            assert(approx(grid.passives()[i,j,"Test"],  dgrid.passives()[i,j,"Test"], 1e-26));
        }
    }
}
#ifdef MHD
void DRAGON_Test::verify_god_dist_grid_2D_MHD(){
    Grid2D grid(32,32,1.0, 1.0);
    DistGrid2D dgrid(32,32,1.0, 1.0);
    grid.boundary = Periodic();
    dgrid.boundary = Periodic();
    grid.passives().add("Test");
    dgrid.passives().add("Test");
    
    for (int i = 0; i <= grid.getSizeX(); i++) {
        for (int j = 0; j <= grid.getSizeY(); j++) {
            grid[i,j] = make_state(1.0+0.1*i+0.1*j, 1.0+0.1*i, 1.0-0.1*j, 0.1*i*j, 10.0-0.1*i+0.1*j);
            grid[i,j].B = vec3{0.1, -0.2, -0.3};
            grid._B()[i,j] = vec3{0.1, -0.2, -0.3};
            
            dgrid[i,j] = grid[i,j];
            dgrid._B()[i,j] = grid._B()[i,j];
            
            grid.passives()[i,j,"Test"] = 10.0-0.1*i+0.1*j;
            dgrid.passives()[i,j,"Test"] = grid.passives()[i,j,"Test"];
        }
    }
    grid.advance(1.0);
    dgrid.advance(1.0);

    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int j = 0; j < grid.getSizeY(); j++) {
            expect_close(grid[i,j], dgrid[i,j], 1e-26);
            assert(approx(grid.passives()[i,j,"Test"],  dgrid.passives()[i,j,"Test"], 1e-26));
        }
    }
    for (int i = 0; i <= grid.getSizeX(); i++) {
        for (int j = 0; j <= grid.getSizeY(); j++) {
            expect_close(grid._B()[i,j], dgrid._B()[i,j], 1e-26);
        }
    }
}
#endif

void DRAGON_Test::verify_god_dist_grid_3D(){
    Grid3D grid(16,16,16,1.0, 1.0,1.0,4);
    DistGrid3D dgrid(16,16,16, 1.0,1.0,1.0,4);
    grid.boundary = Periodic();
    dgrid.boundary = Periodic();
    grid.passives().add("Test");
    dgrid.passives().add("Test");
    
    for (int i = 0; i <= grid.getSizeX(); i++) {
        for (int j = 0; j <= grid.getSizeY(); j++) {
            for (int k = 0; k <= grid.getSizeZ(); k++) {
                grid[i,j,k] = make_state(1.0+0.1*i+0.1*j, 1.0+0.1*i, 1.0-0.1*j, 0.1*k, 10.0-0.1*i+0.1*j-0.1*k);
                dgrid[i,j,k] = grid[i,j,k];
                
                grid.passives()[i,j,k,"Test"] = 10.0-0.1*i+0.01*j-0.1*k;
                dgrid.passives()[i,j,k,"Test"] = grid.passives()[i,j,k,"Test"];
            }
        }
    }
    
    grid.advance(1.0);
    dgrid.advance(1.0);
    
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int j = 0; j < grid.getSizeY(); j++) {
            for (int k = 0; k < grid.getSizeZ(); k++) {
                expect_close(grid[i,j,k], dgrid[i,j,k], 1e-26);
                assert(approx(grid.passives()[i,j,k,"Test"],  dgrid.passives()[i,j,k,"Test"], 1e-26));
            }
        }
    }
}

#ifdef MHD
void DRAGON_Test::verify_god_dist_grid_3D_MHD(){
    Grid3D grid(16,16,16,1.0, 1.0,1.0,3);
    DistGrid3D dgrid(16,16,16,1.0, 1.0,1.0,3);
    grid.boundary = Periodic();
    dgrid.boundary = Periodic();
    grid.passives().add("Test");
    dgrid.passives().add("Test");
    
    for (int i = 0; i <= grid.getSizeX(); i++) {
        for (int j = 0; j <= grid.getSizeY(); j++) {
            for (int k=0; k <= grid.getSizeZ(); k++) {
                grid[i,j,k] = make_state(1.0+0.1*i+0.1*j, 1.0+0.1*i, 1.0-0.1*j, 0.1*k, 10.0-0.1*i+0.1*j-0.1*k);
                grid[i,j,k].B = vec3{0.1, -0.2, -0.3};
                grid._B()[i,j,k] = vec3{0.1, -0.2, -0.3};
                dgrid[i,j,k] = grid[i,j,k];
                dgrid._B()[i,j,k] = grid._B()[i,j,k];
                
                grid.passives()[i,j,k,"Test"] = 10.0-0.1*i+0.1*j-0.1*k;
                dgrid.passives()[i,j,k,"Test"] = grid.passives()[i,j,k,"Test"];
            }
        }
    }
        
    grid.advance(1.0);
    dgrid.advance(1.0);
    
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int j = 0; j < grid.getSizeY(); j++) {
            for (int k = 0; k < grid.getSizeZ(); k++) {
                expect_close(grid[i,j,k], dgrid[i,j,k], 1e-26);
                assert(approx(grid.passives()[i,j,k,"Test"],  dgrid.passives()[i,j,k,"Test"], 1e-26));
            }
        }
    }
    for (int i = 0; i <= grid.getSizeX(); i++) {
        for (int j = 0; j <= grid.getSizeY(); j++) {
            for (int k = 0; k <= grid.getSizeZ(); k++) {
                expect_close(grid._B()[i,j,k], dgrid._B()[i,j,k], 1e-26);
            }
        }
    }
}
#endif
