//
//  MHDGodunov_Tests.cpp
//  DRAGON/Testing/Core/Godunov
//
//  Created by Bobbie Markwick on 09/06/2026.
//

#include "Testing.hpp"
#include "Boundary.hpp"
#include "Grid.hpp"
#include "CFL.hpp"
#include <iostream>

using namespace DRAGON_Test;
using namespace Boundary;

#ifdef MHD

void DRAGON_Test::verify_godunov_1D_MHD(bool output){
    int prev = CONFIG::riemann_choice;
    CONFIG::riemann_choice = RIEMANN_HLLD;
    
    //if(output) std::cout<<"1D MHD Godunov Scheme: \n";
    if(output) std::cout<<"- MHD Uniform Flows: ";
    verify_god_uniform_stationary_1D_MHD();
    verify_god_uniform_moving_1D_MHD();
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- MHD Periodic Conservation: ";
    verify_god_periodic_conservation_1D_MHD();
    if(output) std::cout<<"Passed\n";
    
    CONFIG::riemann_choice = prev;
}
void DRAGON_Test::verify_godunov_2D_MHD(bool output){
    int prev = CONFIG::riemann_choice;
    CONFIG::riemann_choice = RIEMANN_HLLD;
    
    if(output) std::cout<<"- MHD Uniform Flows: ";
    verify_god_uniform_stationary_2D_MHD();
    verify_god_uniform_moving_2D_MHD();
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- MHD Periodic Conservation: ";
    verify_god_periodic_conservation_2D_MHD();
    if(output) std::cout<<"Passed\n";
    
    CONFIG::riemann_choice = prev;
}
void DRAGON_Test::verify_godunov_3D_MHD(bool output){
    int prev = CONFIG::riemann_choice;
    CONFIG::riemann_choice = RIEMANN_HLLD;
    
    if(output) std::cout<<"- MHD Uniform Flows: ";
    verify_god_uniform_stationary_3D_MHD();
    verify_god_uniform_moving_3D_MHD();
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- MHD Periodic Conservation: ";
    verify_god_periodic_conservation_3D_MHD();
    if(output) std::cout<<"Passed\n";
    
    CONFIG::riemann_choice = prev;
}

//MARK: Uniform Flow Tests
void DRAGON_Test::verify_god_uniform_stationary_1D_MHD(){
    Grid1D grid(10, 1.0, 2);
    PrimitiveState W = make_state(1.0, 0.0, 0.0, 0.0, 5.0);
    W.B = {0,2,3};
    for (int i = 0; i < grid.getSize(); i++) grid[i] = W;
    grid.boundary = Reflective();
    grid.advance(1.0);
    for (int i = 0; i < grid.getSize(); i++) expect_close(grid[i], W);
}
void DRAGON_Test::verify_god_uniform_moving_1D_MHD(){
    Grid1D grid(10, 1.0, 2);
    PrimitiveState W = make_state(1.0, 1.0, 2.0, 3.0, 5.0);
    W.B = {1,2,3};
    for (int i = 0; i < grid.getSize(); i++) grid[i] = W;
    grid.boundary = Outflow();
    grid.advance(1.0);
    for (int i = 0; i < grid.getSize(); i++) expect_close(grid[i], W);
}

void DRAGON_Test::verify_god_uniform_stationary_2D_MHD(){
    Grid2D grid(10,10,1.0, 1.0, 2);
    PrimitiveState W = make_state(1.0, 0.0, 0.0, 0.0, 5.0);
    //W.B.z = 1.0;
    for (int i = 0; i <= grid.getSizeX(); i++){
        for (int j = 0; j <= grid.getSizeY(); j++){
            grid[i,j] = W;
            grid.getA()[i,j] = vec3{0,0,0.1*i};
        }
    }
    grid.boundary = Outflow();
    grid.boundary.apply(grid);
    grid.computeBodyAveragedFields();
    W = grid[1,1];
    
    
    grid.advance(1.0);
    
    grid.computeBodyAveragedFields();

    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            expect_close(grid[i,j], W);
        }
    }
}
void DRAGON_Test::verify_god_uniform_moving_2D_MHD(){
    Grid2D grid(10,10,1.0, 1.0, 2);
    PrimitiveState W = make_state(1.0, 1.0, 2.0, 3.0, 5.0);
    W.B.z = 1.0;
    for (int i = 0; i <= grid.getSizeX(); i++){
        for (int j = 0; j <= grid.getSizeY(); j++){
            grid[i,j] = W;
            grid.getA()[i,j] = vec3{0,0,0.2*i};
        }
    }
    grid.boundary = Outflow();
    grid.boundary.apply(grid);
    grid.computeBodyAveragedFields();
    W = grid[1,1];
    
    grid.advance(1.0);
    
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            expect_close(grid[i,j], W);
        }
    }
}

void DRAGON_Test::verify_god_uniform_stationary_3D_MHD(){
    Grid3D grid(10,10,10, 1.0,1.0,1.0,2);
    PrimitiveState W = make_state(1.0, 0.0, 0.0, 0.0, 5.0);
    for (int i = 0; i <= grid.getSizeX(); i++){
        for (int j = 0; j <= grid.getSizeY(); j++){
            for(int k = 0; k <= grid.getSizeZ(); k++){
                grid[i,j,k] = W;
                grid.getA()[i,j,k] = vec3{0,0,0.1*i};
            }
        }
    }
    grid.boundary = Reflective();
    grid.boundary.apply(grid);
    grid.computeBodyAveragedFields();
    W = grid[1,1,1];
    
    grid.advance(1.0);
    
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            for(int k=0; k<grid.getSizeZ(); k++){
                expect_close(grid[i,j,k], W);
            }
        }
    }
}
void DRAGON_Test::verify_god_uniform_moving_3D_MHD(){
    Grid3D grid(10,10,10, 1.0,1.0, 1.0, 2);
    PrimitiveState W = make_state(1.0, 1.0, 2.0, 3.0, 5.0);
    for (int i = 0; i <= grid.getSizeX(); i++){
        for (int j = 0; j <= grid.getSizeY(); j++){
            for (int k = 0; k <= grid.getSizeZ(); k++){
                grid[i,j,k] = W;
                grid.getA()[i,j,k] = vec3{0,0,0.125*i};
            }
        }
    }
    grid.computeBodyAveragedFields();
    W = grid[2,2,2];

    grid.boundary = Fixed(W);
    
    grid.advance(1.0);
    
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            for (int k = 0; k < grid.getSizeZ(); k++){
                expect_close(grid[i,j,k], W);
            }
        }
    }
}

//MARK: Periodic Conservation Tests
void DRAGON_Test::verify_god_periodic_conservation_1D_MHD(){
    Grid1D grid(10, 1.0, 2);
    for (int i = 0; i < grid.getSize(); i++){
        grid[i] = make_state(1.0+0.1*i, 1.0+0.1*i, -0.1*i, 0.1*i, 10.0-0.1*i);
        grid[i].B = vec3{0, -0.2*i, -0.3*i};
    }
    grid.boundary = Periodic();
    grid.boundary.apply(grid);
    
    ConservativeState expected = ConservativeState();
    for (int i = 0; i < grid.getSize(); i++) expected += ConservativeState(grid[i]);
    
    grid.advance(10.0);
    
    ConservativeState got = ConservativeState();
    for (int i = 0; i < grid.getSize(); i++) got += ConservativeState(grid[i]);
    
    expect_close(expected, got);

}

void DRAGON_Test::verify_god_periodic_conservation_2D_MHD(){
    Grid2D grid(10,10,1.0, 1.0, 2);
    grid.boundary = Periodic();
    for (int i = 0; i <= grid.getSizeX(); i++){
        for (int j = 0; j <= grid.getSizeY(); j++){
            grid[i,j] = make_state(1.0+0.1*i+0.1*j, 1.0+0.1*i, 1.0-0.1*j, 0.1*i*j, 10.0-0.1*i+0.1*j);
            grid[i,j].B = vec3{0.1*i, -0.2*j, -0.3*i};
            grid.getA()[i,j] = vec3{0.1*i, -0.2*j, -0.3*i};
            
        }
    }
    grid.boundary.apply(grid);
    grid.computeBodyAveragedFields();
    
    ConservativeState expected = ConservativeState();
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            expected += ConservativeState(grid[i,j]);
        }
    }
    expected.B = {0,0,0};//B isn't a conserved quanitty
    
    grid.advance(10.0);
    
    ConservativeState got = ConservativeState();
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            got += ConservativeState(grid[i,j]);
            assert((grid[i,j].rho > 0));
        }
    }
    got.B = {0,0,0};//B isn't a conserved quanitty
    
    expect_close(expected, got);
}
void DRAGON_Test::verify_god_periodic_conservation_3D_MHD(){
    Grid3D grid(10,10,10,1.0, 1.0,1.0, 2);
    grid.boundary = Periodic();
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            for(int k = 0; k < grid.getSizeZ(); k++){
                grid[i,j,k] = make_state(1.0+0.1*i+0.1*j, 1.0+0.1*i, 1.0-0.1*j, 0.1*k, 10.0-0.1*i+0.1*j-0.1*k);
                grid.getA()[i,j,k] = vec3{0.1*i, -0.2*j, -0.3*k};
            }
        }
    }
    grid.boundary.apply(grid);
    grid.computeBodyAveragedFields();

    
    ConservativeState expected = ConservativeState();
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            for(int k=0; k<grid.getSizeZ(); k++){
                expected += ConservativeState(grid[i,j,k]);
            }
        }
    }
    expected.B = {0,0,0};//B isn't a conserved quanitty

    
    grid.advance(2.0);
    
    ConservativeState got = ConservativeState();
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            for(int k=0; k<grid.getSizeZ(); k++){
                got += ConservativeState(grid[i,j,k]);
                assert((grid[i,j,k].rho > 0));
            }
        }
    }
    got.B = {0,0,0};//B isn't a conserved quanitty

    
    expect_close(expected, got);
}
#endif
