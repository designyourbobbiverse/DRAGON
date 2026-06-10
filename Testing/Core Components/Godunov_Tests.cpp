//
//  Godunov_Tests.cpp
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

void DRAGON_Test::verify_godunov_1D(bool output){
    verify_god_uniform_stationary_1D();
    verify_god_uniform_moving_1D();
    verify_god_periodic_conservation_1D();
    verify_god_dt0_1D();
    if(output) std::cout<<"All 1D Godunov Scheme Tests Passed\n";
}
void DRAGON_Test::verify_godunov_2D_Split(bool output){
    verify_god_uniform_stationary_2D_Split();
    verify_god_uniform_moving_2D_Split();
    verify_god_periodic_conservation_2D_Split();
    verify_god_dt0_2D_Split();
    verify_split_2D_X_match_1D();
    verify_split_2D_Y_match_1D();
    if(output) std::cout<<"All 2D Split Scheme Tests Passed\n";
}
void DRAGON_Test::verify_godunov_3D_Split(bool output){
    verify_god_uniform_stationary_3D_Split();
    verify_god_uniform_moving_3D_Split();
    verify_god_periodic_conservation_3D_Split();
    verify_god_dt0_3D_Split();
    verify_split_3D_X_match_1D();
    verify_split_3D_Y_match_1D();
    verify_split_3D_Z_match_1D();
    if(output) std::cout<<"All 3D Split Scheme Tests Passed\n";
    
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

void DRAGON_Test::verify_god_uniform_stationary_2D_Split(){
    Grid2D grid(10,10,1.0, 1.0, 2);
    PrimitiveState W = make_state(1.0, 0.0, 0.0, 0.0, 5.0);
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            grid[i,j] = W;
        }
    }
    grid.boundary = Reflective();
    grid.advance(1.0);
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            expect_close(grid[i,j], W);
        }
    }
}
void DRAGON_Test::verify_god_uniform_moving_2D_Split(){
    Grid2D grid(10,10,1.0, 1.0, 2);
    PrimitiveState W = make_state(1.0, 1.0, 2.0, 3.0, 5.0);
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            grid[i,j] = W;
        }
    }
    grid.boundary = Outflow();
    
    grid.advance(1.0);
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            expect_close(grid[i,j], W);
        }
    }
}

void DRAGON_Test::verify_god_uniform_stationary_3D_Split(){
    Grid3D grid(10,10,10,1.0, 1.0,1.0, 2);
    PrimitiveState W = make_state(1.0, 0.0, 0.0, 0.0, 5.0);
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            for(int k=0; k<grid.getSizeZ(); k++){
                grid[i,j,k] = W;
            }
        }
    }
    grid.boundary = Reflective();
    grid.advance(1.0);
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            for(int k=0; k<grid.getSizeZ(); k++){
                expect_close(grid[i,j,k], W);
            }
        }
    }
}
void DRAGON_Test::verify_god_uniform_moving_3D_Split(){
    Grid2D grid(10,10,1.0, 1.0, 2);
    PrimitiveState W = make_state(1.0, 1.0, 2.0, 3.0, 5.0);
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            grid[i,j] = W;
        }
    }
    grid.boundary = Outflow();
    
    grid.advance(1.0);
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            expect_close(grid[i,j], W);
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

void DRAGON_Test::verify_god_periodic_conservation_2D_Split(){
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
    grid.advance(10.0);
    
    ConservativeState got = ConservativeState();
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            got += ConservativeState(grid[i,j]);
            assert((grid[i,j].rho > 0));
        }
    }
    expect_close(expected, got);
}
void DRAGON_Test::verify_god_periodic_conservation_3D_Split(){
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
    grid.advance(1.0);
    
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

void DRAGON_Test::verify_god_dt0_2D_Split(){
    Grid2D grid(10,10,1.0, 1.0, 2), expected(10,10,1.0, 1.0, 2);
    grid.boundary = Periodic();
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            grid[i,j] = make_state(1.0+0.1*i+0.1*j, 1.0+0.1*i, 1.0-0.1*j, 0.1*i*j, 10.0-0.1*i+0.1*j);
            expected[i,j] = grid[i,j];
        }
    }

    grid.advance(0.0);
    
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            expect_close(expected[i,j],grid[i,j]);
        }
    }
}
void DRAGON_Test::verify_god_dt0_3D_Split(){
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
    grid.advance(0.0);
    
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            for(int k=0; k<grid.getSizeZ(); k++){
                expect_close(expected[i,j,k],grid[i,j,k]);
            }
        }
    }
}

//MARK: Split-scheme Match Tests
void DRAGON_Test::verify_split_2D_X_match_1D(){
    Grid2D grid(10,10,1.0,1.0,2);
    Grid1D expected(10, 1.0, 2);
    grid.boundary = Outflow() + Periodic(Y);
    expected.boundary = Outflow();;

    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            grid[i,j] = make_state(1.0+0.1*i, 1.0+0.1*i, -0.1*i, 0.1*i, 10.0-0.1*i);
        }
        expected[i] = grid[i,0];
    }
        
    double dt = 0.01;
    grid.advance(dt);
    expected.advance(dt/2);expected.advance(dt/2);
        
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            expect_close(grid[i,j], expected[i]);
        }
    }

}
void DRAGON_Test::verify_split_2D_Y_match_1D(){
    Grid2D grid(10,10,1.0,1.0,2);
    Grid1D expected(10, 1.0, 2);
    grid.boundary = Outflow() + Periodic(X);
    expected.boundary = Outflow();;

    for (int j = 0; j < grid.getSizeY(); j++){
        for (int i = 0; i < grid.getSizeX(); i++){
            grid[i,j] = make_state(1.0+0.1*j, 1.0+0.1*j, -0.1*j, 0.1*j, 10.0-0.1*j);
        }
        expected[j] = grid[0,j];
        expected[j].swapXY();
    }
        
    double dt = 0.01;
    grid.advance(dt);
    expected.advance(dt);
        
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            expected[j].swapXY();
            expect_close(grid[i,j], expected[j]);
        }
    }

}
void DRAGON_Test::verify_split_3D_X_match_1D(){
    Grid3D grid(6,3,3,1.0,1.0,1.0,2);
    Grid1D expected(6, 1.0, 2);
    grid.boundary = Outflow() + Periodic(Y|Z);
    expected.boundary = Outflow();;

    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            for (int k = 0; k < grid.getSizeZ(); k++){
                grid[i,j,k] = make_state(1.0+0.1*i, 1.0+0.1*i, -0.1*i, 0.1*i, 10.0-0.1*i);
            }
        }
        expected[i] = grid[i,0,0];
    }
        
    double dt = 0.01;
    grid.advance(dt);
    expected.advance(dt/2);expected.advance(dt/2);
        
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            for (int k = 0; k < grid.getSizeZ(); k++){
                expect_close(grid[i,j,k], expected[i]);
            }
        }
    }

}
void DRAGON_Test::verify_split_3D_Y_match_1D(){
    Grid3D grid(3,6,3,1.0,1.0,1.0,2);
    Grid1D expected(6, 1.0, 2);
    grid.boundary = Outflow() + Periodic(X|Z);
    expected.boundary = Outflow();;

    for (int j = 0; j < grid.getSizeY(); j++){
        for (int i = 0; i < grid.getSizeX(); i++){
            for (int k = 0; k < grid.getSizeZ(); k++){
                grid[i,j,k] = make_state(1.0+0.1*j, 1.0+0.1*j, -0.1*j, 0.1*j, 10.0-0.1*j);
            }
        }
        expected[j] = grid[0,j,0];
        expected[j].swapXY();
    }
        
    double dt = 0.01;
    grid.advance(dt);
    expected.advance(dt/2);expected.advance(dt/2);
        
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            for (int k = 0; k < grid.getSizeZ(); k++){
                expected[j].swapXY();
                expect_close(grid[i,j,k], expected[j]);
            }
        }
    }
}
void DRAGON_Test::verify_split_3D_Z_match_1D(){
    Grid3D grid(3,3,6,1.0,1.0,1.0,2);
    Grid1D expected(6, 1.0, 2);
    grid.boundary = Outflow() + Periodic(X|Y);
    expected.boundary = Outflow();;

    for (int k = 0; k < grid.getSizeZ(); k++){
        for (int j = 0; j < grid.getSizeY(); j++){
            for (int i = 0; i < grid.getSizeX(); i++){
                grid[i,j,k] = make_state(1.0+0.1*k, 1.0+0.1*k, -0.1*k, 0.1*k, 10.0-0.1*k);
            }
        }
        expected[k] = grid[0,0,k];
        expected[k].swapXZ();
    }
        
    double dt = 0.001;
    grid.advance(dt);
    expected.advance(dt);
        
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            for (int k = 0; k < grid.getSizeZ(); k++){
                expected[k].swapXZ();
                expect_close(grid[i,j,k], expected[k]);
            }
        }
    }

}
