//
//  Godunov_Tests.cpp
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

void DRAGON_Test::verify_godunov_1D(bool output){
    if(output) std::cout<<"1D Godunov Scheme: \n";
    if(output) std::cout<<"- Uniform Flows: ";
    verify_god_uniform_stationary_1D();
    verify_god_uniform_moving_1D();
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- Periodic Conservation: ";
    verify_god_periodic_conservation_1D();
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- Zero Time: ";
    verify_god_dt0_1D();
    if(output) std::cout<<"Passed\n";
}
void DRAGON_Test::verify_godunov_2D_Split(bool output){
    if(output) std::cout<<"2D Split Scheme: \n";
    if(output) std::cout<<"- Uniform Flows: ";
    verify_god_uniform_stationary_2D(true);
    verify_god_uniform_moving_2D(true);
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- Periodic Conservation: ";
    verify_god_periodic_conservation_2D(true);
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- Zero Time: ";
    verify_god_dt0_2D(true);
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- 1D Match (X): ";
    verify_2D_X_match_1D(true);
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- 1D Match (Y): ";
    verify_2D_Y_match_1D(true);
    if(output) std::cout<<"Passed\n";
}
void DRAGON_Test::verify_godunov_2D_Unsplit(bool output){
    CONFIG::riemann_choice = RIEMANN_EXACT;
    if(output) std::cout<<"2D Unsplit Scheme: \n";
    if(output) std::cout<<"- Uniform Flows: ";
    verify_god_uniform_stationary_2D(false);
    verify_god_uniform_moving_2D(false);
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- Periodic Conservation: ";
    verify_god_periodic_conservation_2D(false);
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- Zero Time: ";
    verify_god_dt0_2D(false);
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- 1D Match (X): ";
    verify_2D_X_match_1D(false);
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- 1D Match (Y): ";
    verify_2D_Y_match_1D(false);
    if(output) std::cout<<"Passed\n";
#ifdef CTU
    if(output) std::cout<<"- Diagonal Contact: ";
    verify_ctu_diagonal_contact_2D();
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- Blast: ";
    verify_ctu_blast_2D();
    if(output) std::cout<<"Passed\n";
#endif
}
void DRAGON_Test::verify_godunov_3D_Split(bool output){
    if(output) std::cout<<"3D Split Scheme: \n";
    if(output) std::cout<<"- Uniform Flows: ";
    verify_god_uniform_stationary_3D(true);
    verify_god_uniform_moving_3D(true);
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- Periodic Conservation: ";
    verify_god_periodic_conservation_3D(true);
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- Zero Time: ";
    verify_god_dt0_3D(true);
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- 1D Match (X): ";
    verify_3D_X_match_1D(true);
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- 1D Match (Y): ";
    verify_3D_Y_match_1D(true);
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- 1D Match (Z): ";
    verify_3D_Z_match_1D(true);
    if(output) std::cout<<"Passed\n";
}
void DRAGON_Test::verify_godunov_3D_Unsplit(bool output){
    CONFIG::riemann_choice = RIEMANN_HLLC;
    if(output) std::cout<<"3D Unsplit Scheme: \n";
    if(output) std::cout<<"- Uniform Flows: ";
    verify_god_uniform_stationary_3D(false);
    verify_god_uniform_moving_3D(false);
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- Periodic Conservation: ";
    verify_god_periodic_conservation_3D(false);
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- Zero Time: ";
    verify_god_dt0_3D(false);
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- 1D Match (X): ";
    verify_3D_X_match_1D(false);
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- 1D Match (Y): ";
    verify_3D_Y_match_1D(false);
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- 1D Match (Z): ";
    verify_3D_Z_match_1D(false);
    if(output) std::cout<<"Passed\n";
#ifdef CTU
    if(output) std::cout<<"- Diagonal Contact: ";
    verify_ctu_diagonal_contact_3D();
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- Blast: ";
    verify_ctu_blast_3D();
    if(output) std::cout<<"Passed\n";
#endif
}

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
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            grid[i,j] = W;
        }
    }
    grid.boundary = Reflective();
    
    if(split) grid.advance_split(1.0);
    else grid.advance_unsplit(1.0);
    
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            expect_close(grid[i,j], W);
        }
    }
}
void DRAGON_Test::verify_god_uniform_moving_2D(bool split){
    Grid2D grid(10,10,1.0, 1.0, 2);
    PrimitiveState W = make_state(1.0, 1.0, 2.0, 3.0, 5.0);
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            grid[i,j] = W;
        }
    }
    grid.boundary = Outflow();
    
    if(split) grid.advance_split(1.0);
    else grid.advance_unsplit(1.0);
    
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            expect_close(grid[i,j], W);
        }
    }
}

void DRAGON_Test::verify_god_uniform_stationary_3D(bool split){
    Grid3D grid(10,10,10, 1.0,1.0,1.0,2);
    PrimitiveState W = make_state(1.0, 0.0, 0.0, 0.0, 5.0);
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            for(int k=0; k<grid.getSizeZ(); k++){
                grid[i,j,k] = W;
            }
        }
    }
    grid.boundary = Reflective();
    
    if(split) grid.advance_split(1.0);
    else grid.advance_unsplit(1.0);
    
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            for(int k=0; k<grid.getSizeZ(); k++){
                expect_close(grid[i,j,k], W);
            }
        }
    }
}
void DRAGON_Test::verify_god_uniform_moving_3D(bool split){
    Grid3D grid(10,10,10, 1.0,1.0, 1.0, 2);
    PrimitiveState W = make_state(1.0, 1.0, 2.0, 3.0, 5.0);
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            for (int k = 0; k < grid.getSizeZ(); k++){
                grid[i,j,k] = W;
            }
        }
    }
    grid.boundary = Outflow();
    
    if(split) grid.advance_split(1.0);
    else grid.advance_unsplit(1.0);
    
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            for (int k = 0; k < grid.getSizeZ(); k++){
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
    
    expect_close(expected, got);

}

void DRAGON_Test::verify_god_periodic_conservation_2D(bool split){
    Grid2D grid(10,10,1.0, 1.0, 2);
    grid.boundary = Periodic();
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            grid[i,j] = make_state(1.0+0.1*i+0.1*j, 1.0+0.1*i, 1.0-0.1*j, 0.1*i*j, 10.0-0.1*i+0.1*j);
        }
    }
    
    ConservativeState expected = ConservativeState();
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            expected += ConservativeState(grid[i,j]);
        }
    }
    
    if(split) grid.advance_split(10.0);
    else grid.advance_unsplit(10.0);
    
    ConservativeState got = ConservativeState();
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            got += ConservativeState(grid[i,j]);
            assert((grid[i,j].rho > 0));
        }
    }
    expect_close(expected, got);
}
void DRAGON_Test::verify_god_periodic_conservation_3D(bool split){
    Grid3D grid(10,10,10,1.0, 1.0,1.0, 2);
    grid.boundary = Periodic();
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            for(int k=0; k<grid.getSizeZ(); k++){
                grid[i,j,k] = make_state(1.0+0.1*i+0.1*j, 1.0+0.1*i, 1.0-0.1*j, 0.1*k, 10.0-0.1*i+0.1*j-0.1*k);
            }
        }
    }
    
    ConservativeState expected = ConservativeState();
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            for(int k=0; k<grid.getSizeZ(); k++){
                expected += ConservativeState(grid[i,j,k]);
            }
        }
    }
    
    if(split) grid.advance_split(2.0);
    else grid.advance_unsplit(2.0);
    
    ConservativeState got = ConservativeState();
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            for(int k=0; k<grid.getSizeZ(); k++){
                got += ConservativeState(grid[i,j,k]);
                assert((grid[i,j,k].rho > 0));
            }
        }
    }
    expect_close(expected, got);
}
//MARK: dt=0 Tests
void DRAGON_Test::verify_god_dt0_1D(){
    Grid1D grid(10, 1.0, 2), expected(10, 1.0, 2);
    grid.boundary = Periodic();

    for (int i = 0; i < grid.getSize(); i++){
        grid[i] = make_state(1.0+0.1*i, 1.0+0.1*i, -0.1*i, 0.1*i, 10.0-0.1*i);
        expected[i] = grid[i];
    }
        
    grid.advance(0.0);
        
    for (int i = 0; i < grid.getSize(); i++){
        expect_close(expected[i],grid[i]);
    }

}

void DRAGON_Test::verify_god_dt0_2D(bool split){
    Grid2D grid(10,10,1.0, 1.0, 2), expected(10,10,1.0, 1.0, 2);
    grid.boundary = Periodic();
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            grid[i,j] = make_state(1.0+0.1*i+0.1*j, 1.0+0.1*i, 1.0-0.1*j, 0.1*i*j, 10.0-0.1*i+0.1*j);
            expected[i,j] = grid[i,j];
        }
    }

    if(split) grid.advance_split(0.0);
    else grid.advance_unsplit(0.0);
    
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            expect_close(expected[i,j],grid[i,j]);
        }
    }
}
void DRAGON_Test::verify_god_dt0_3D(bool split){
    Grid3D grid(10,10,10,1.0, 1.0,1.0, 2), expected(10,10,10,1.0, 1.0,1.0, 2);
    grid.boundary = Periodic();
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            for(int k=0; k<grid.getSizeZ(); k++){
                grid[i,j,k] = make_state(1.0+0.1*i+0.1*j, 1.0+0.1*i, 1.0-0.1*j, 0.1*k, 10.0-0.1*i+0.1*j-0.1*k);
                expected[i,j,k] = grid[i,j,k];
            }
        }
    }
    if(split) grid.advance_split(0.0);
    else grid.advance_unsplit(0.0);
    
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            for(int k=0; k<grid.getSizeZ(); k++){
                expect_close(expected[i,j,k],grid[i,j,k]);
            }
        }
    }
}
