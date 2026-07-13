//
//  Godunov_Tests.cpp
//  DRAGON/Testing/Core/Godunov
//
//  Created by Bobbie Markwick on 09/06/2026.
//

#include "Testing.hpp"
#include "Boundary.hpp"
#include "DistGrid.hpp"
#include "DragonWing.hpp"

using namespace DRAGON_Test;
using namespace Boundary;

void DRAGON_Test::verify_god_dist_grid_1D(){
    Grid1D grid(100, 1.0);
    DistGrid1D dgrid(100, 1.0);
    for (int i = 0; i < grid.getSize(); i++) {
        grid[i] = make_state(1.0+0.01*i, 1.0+0.01*i, -0.01*i, 0.01*i, 10.0-0.01*i);
        #ifdef MHD
        grid[i].B = vec3{0, -0.2*i, -0.3*i};
        #endif
        dgrid[i] = grid[i];
    }
        
    grid.boundary = Periodic();
    dgrid.boundary = Periodic();
    
    grid.advance(1.0);
    dgrid.advance(1.0);
    
    for (int i = 0; i < grid.getSize(); i++){
        expect_close(grid[i], dgrid[i]);
    }
    
    DRAGONWING::initialize(0);
}


void DRAGON_Test::verify_god_dist_grid_2D(){
    Grid2D grid(32,32,1.0, 1.0, 4);
    DistGrid2D dgrid(32,32,1.0, 1.0, 4);
    grid.boundary = Periodic();
    dgrid.boundary = Periodic();
    
    for (int i = 0; i <= grid.getSizeX(); i++){
        for (int j = 0; j <= grid.getSizeY(); j++){
            grid[i,j] = make_state(1.0+0.1*i+0.1*j, 1.0+0.1*i, 1.0-0.1*j, 0.1*i*j, 10.0-0.1*i+0.1*j);
            dgrid[i,j] = grid[i,j];
        }
    }
    
    grid.advance(1.0);
    dgrid.advance(1.0);
    
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            expect_close(grid[i,j], dgrid[i,j]);
        }
    }
    
    DRAGONWING::initialize(0);
}
#ifdef MHD
void DRAGON_Test::verify_god_dist_grid_2D_MHD(){
    Grid2D grid(32,32,1.0, 1.0,3);
    DistGrid2D dgrid(32,32,1.0, 1.0,3);
    grid.boundary = Periodic();
    dgrid.boundary = Periodic();
    
    for (int i = 0; i <= grid.getSizeX(); i++){
        for (int j = 0; j <= grid.getSizeY(); j++){
            grid[i,j] = make_state(1.0+0.1*i+0.1*j, 1.0+0.1*i, 1.0-0.1*j, 0.1*i*j, 10.0-0.1*i+0.1*j);
            grid[i,j].B = {0,0.03,1};
            grid._A()[i,j] = vec3{0,0, -0.03*i};
            dgrid[i,j] = grid[i,j];
            dgrid._A()[i,j] = grid._A()[i,j];
        }
    }
    grid.advance(1.0);
    dgrid.advance(1.0);
    
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            expect_close(grid[i,j], dgrid[i,j]);
            expect_close(grid._A()[i,j].z, dgrid._A()[i,j].z);
        }
    }
    
    DRAGONWING::initialize(0);
}
#endif

void DRAGON_Test::verify_god_dist_grid_3D(){
    Grid3D grid(16,16,16,1.0, 1.0,1.0,4);
    DistGrid3D dgrid(16,16,16, 1.0,1.0,1.0,4);
    grid.boundary = Periodic();
    dgrid.boundary = Periodic();
    
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            for(int k = 0; k < grid.getSizeZ(); k++){
                grid[i,j,k] = make_state(1.0+0.1*i+0.1*j, 1.0+0.1*i, 1.0-0.1*j, 0.1*k, 10.0-0.1*i+0.1*j-0.1*k);
                dgrid[i,j,k] = grid[i,j,k];
            }
        }
    }
    
    grid.advance(1.0);
    dgrid.advance(1.0);
    
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            for(int k = 0; k < grid.getSizeZ(); k++){
                expect_close(grid[i,j,k], dgrid[i,j,k]);
            }
        }
    }
    
    DRAGONWING::initialize(0);
}

#ifdef MHD
void DRAGON_Test::verify_god_dist_grid_3D_MHD(){
    Grid3D grid(16,16,16,1.0, 1.0,1.0,3);
    DistGrid3D dgrid(16,16,16,1.0, 1.0,1.0,3);
    grid.boundary = Periodic();
    dgrid.boundary = Periodic();
    
    for (int i = 0; i <= grid.getSizeX(); i++){
        for (int j = 0; j <= grid.getSizeY(); j++){
            for(int k=0; k <= grid.getSizeZ(); k++){
                grid[i,j,k] = make_state(1.0+0.1*i+0.1*j, 1.0+0.1*i, 1.0-0.1*j, 0.1*k, 10.0-0.1*i+0.1*j-0.1*k);
                grid[i,j,k].B = vec3{0.1*i, -0.2*j, -0.3*k};
                grid._A()[i,j,k] = vec3{0.1*j, -0.2*k, -0.3*i};
                dgrid[i,j,k] = grid[i,j,k];
                dgrid._A()[i,j,k] = grid._A()[i,j,k];
            }
        }
    }
    
    grid.advance(0.00001);
    dgrid.advance(0.00001);
    
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            for(int k = 0; k < grid.getSizeZ(); k++){
                expect_close(grid[i,j,k], dgrid[i,j,k]);
                expect_close(grid._A()[i,j,k], dgrid._A()[i,j,k]);
            }
        }
    }
    
    DRAGONWING::initialize(0);
}
#endif
