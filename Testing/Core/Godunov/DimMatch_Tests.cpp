//
//  Core/Godunov/DimMatch_Tests.cpp
//  DRAGON/Testing
//
//  Created by Bobbie Markwick on 09/06/2026.
//

#include "Testing.hpp"
#include "Grid.hpp"
#include "CFL.hpp"
#include <iostream>

using namespace DRAGON_Test;
using namespace Boundary;

//MARK: 2D Matches 1D
void DRAGON_Test::verify_2D_X_match_1D(bool split){
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
    if(split){
        grid.advance_split(dt,false);
        expected.advance(dt/2,false);expected.advance(dt/2,false);
    } else{
        grid.advance_unsplit(dt,false);
        expected.advance(dt,false);
    }
        
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            expect_close(grid[i,j], expected[i]);
        }
    }
    if(split){
        grid.advance_split(dt,false);
        expected.advance(dt,false);
    } else{
        grid.advance_unsplit(dt,false);
        expected.advance(dt,false);
    }
        
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            expect_close(grid[i,j], expected[i], 1e-8);
        }
    }

}
void DRAGON_Test::verify_2D_Y_match_1D(bool split){
    Grid2D grid(10,10,1.0,1.0,2);
    Grid1D expected(10, 1.0, 2);
    grid.boundary = Outflow() + Periodic(X);
    expected.boundary = Outflow();;

    for (int j = 0; j < grid.getSizeY(); j++){
        for (int i = 0; i < grid.getSizeX(); i++){
            grid[i,j] = make_state(1.0+0.1*j, 1.0+0.1*j, -0.1*j, 0.1*j, 10.0-0.1*j);
        }
        expected[j] = grid[0,j].swappedXY();
    }
        
    double dt = 0.01;
    if(split){
        grid.advance_split(dt,false);
        expected.advance(dt,false);
    } else{
        grid.advance_unsplit(dt,false);
        expected.advance(dt,false);
    }
        
    for (int j = 0; j < grid.getSizeY(); j++){
        for (int i = 0; i < grid.getSizeX(); i++){
            expect_close(grid[i,j], expected[j].swappedXY());
        }
    }
    
    if(split){
        grid.advance_split(dt,false);
        expected.advance(dt/2,false);expected.advance(dt/2,false);
    } else{
        grid.advance_unsplit(dt,false);
        expected.advance(dt,false);
    }
        
    for (int j = 0; j < grid.getSizeY(); j++){
        for (int i = 0; i < grid.getSizeX(); i++){
            expect_close(grid[i,j], expected[j].swappedXY(), 1e-8);
        }
    }

}

//MARK: 3D Matches 1D
void DRAGON_Test::verify_3D_X_match_1D(bool split){
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
    if(split){
        grid.advance_split(dt,false);
        expected.advance(dt/2,false);expected.advance(dt/2,false);
    } else {
        grid.advance_unsplit(dt,false);
        expected.advance(dt,false);
    }
        
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            for (int k = 0; k < grid.getSizeZ(); k++){
                expect_close(grid[i,j,k], expected[i]);
            }
        }
    }
    if(split){
        grid.advance_split(dt,false);
        expected.advance(dt/2,false);expected.advance(dt/2,false);
        grid.advance_split(dt,false);
        expected.advance(dt,false);
        grid.advance_split(dt,false);
        expected.advance(dt,false);
        grid.advance_split(dt,false);
        expected.advance(dt/2,false);expected.advance(dt/2,false);
        grid.advance_split(dt,false);
        expected.advance(dt/2,false);expected.advance(dt/2,false);
    } else {
        grid.advance_unsplit(dt,false);
        expected.advance(dt,false);
        grid.advance_unsplit(dt,false);
        expected.advance(dt,false);
    }
    
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            for (int k = 0; k < grid.getSizeZ(); k++){
                expect_close(grid[i,j,k], expected[i], 1e-8);
            }
        }
    }
    

}
void DRAGON_Test::verify_3D_Y_match_1D(bool split){
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
        expected[j] = grid[0,j,0].swappedXY();
    }
        
    double dt = 0.01;
    if(split){
        grid.advance_split(dt,false);
        expected.advance(dt/2,false);expected.advance(dt/2,false);
    } else {
        grid.advance_unsplit(dt,false);
        expected.advance(dt,false);
    }
        
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            for (int k = 0; k < grid.getSizeZ(); k++){
                expect_close(grid[i,j,k], expected[j].swappedXY());
            }
        }
    }
    if(split){
        grid.advance_split(dt,false);
        expected.advance(dt,false);
        grid.advance_split(dt,false);
        expected.advance(dt/2,false);expected.advance(dt/2,false);
        grid.advance_split(dt,false);
        expected.advance(dt/2,false);expected.advance(dt/2,false);
        grid.advance_split(dt,false);
        expected.advance(dt,false);
        grid.advance_split(dt,false);
        expected.advance(dt/2,false);expected.advance(dt/2,false);
    } else {
        grid.advance_unsplit(dt,false);
        expected.advance(dt,false);
    }
    
    
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            for (int k = 0; k < grid.getSizeZ(); k++){
                expect_close(grid[i,j,k], expected[j].swappedXY(), 1e-8);
            }
        }
    }
}
void DRAGON_Test::verify_3D_Z_match_1D(bool split){
    Grid3D grid(3,3,6,1.0,1.0,1.0,2);
    Grid1D expected(6, 1.0, 2);
    grid.boundary = Outflow() + Periodic(X|Y);
    expected.boundary = Outflow();

    for (int k = 0; k < grid.getSizeZ(); k++){
        for (int j = 0; j < grid.getSizeY(); j++){
            for (int i = 0; i < grid.getSizeX(); i++){
                grid[i,j,k] = make_state(1.0+0.1*k, 1.0+0.1*k, -0.1*k, 0.1*k, 10.0-0.1*k);
            }
        }
        expected[k] = grid[0,0,k].swappedXZ();
    }
        
    double dt = 0.01;
    if(split){
        grid.advance_split(dt,false);
        expected.advance(dt,false);
    } else {
        grid.advance_unsplit(dt,false);
        expected.advance(dt,false);
    }
        
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            for (int k = 0; k < grid.getSizeZ(); k++){
                expect_close(grid[i,j,k], expected[k].swappedXZ());
            }
        }
    }
    
    if(split){
        grid.advance_split(dt,false);
        expected.advance(dt/2,false);expected.advance(dt/2,false);
        grid.advance_split(dt,false);
        expected.advance(dt/2,false);expected.advance(dt/2,false);
        grid.advance_split(dt,false);
        expected.advance(dt/2,false);expected.advance(dt/2,false);
        grid.advance_split(dt,false);
        expected.advance(dt/2,false);expected.advance(dt/2,false);
        grid.advance_split(dt,false);
        expected.advance(dt,false);
    } else {
        grid.advance_unsplit(dt,false);
        expected.advance(dt,false);
        grid.advance_unsplit(dt,false);
        expected.advance(dt,false);
    }
        
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            for (int k = 0; k < grid.getSizeZ(); k++){
                expect_close(grid[i,j,k], expected[k].swappedXZ(),1e-8);
            }
        }
    }

}
