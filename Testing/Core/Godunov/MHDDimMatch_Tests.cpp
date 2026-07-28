//
//  DimMatch_Tests.cpp
//  DRAGON/Testing/Core/Godunov
//
//  Created by Bobbie Markwick on 09/06/2026.
//

#include "Testing.hpp"
#include "Grid.hpp"
#include "CFL.hpp"
#include <iostream>

using namespace DRAGON_Test;
using namespace Boundary;

//Dai-Woodward Shock Tube Parameters
constexpr double rho_L = 1.08;
constexpr double rho_R = 1.0;
constexpr double vx = 1.2;
constexpr double vy = 0.01;
constexpr double vz = 0.5;
constexpr double Bx = 2.0;
constexpr double ByL = 3.6;
constexpr double ByR = 4.0;
constexpr double Bz = 2.0;
constexpr double p_L = 0.95;
constexpr double p_R = 1.0;
constexpr int nx = 32;
constexpr int ny = 8;
constexpr double dx = 1.0/nx;

//MARK: 2D - X
void DRAGON_Test::verify_2D_X_match_1D_MHD(){
    /*
    auto expected = Grid2D(nx, ny, dx, dx);
    expected.boundary =  Boundary::Outflow("X") + Boundary::Periodic("YZ");

    for(int i = 0; i<nx; i++){
        for(int j = 0; j<ny; j++){
            expected[i,j].rho = i < nx/2 ? rho_L : rho_R;
            expected[i,j].p = i < nx/2 ? p_L : p_R;
            expected[i,j].v = i < nx/2 ? vec3{vx,vy,vz} : vec3{0,0,0};
            expected[i,j].B = i < nx/2 ? vec3{Bx,ByL,Bz} : vec3{Bx,ByR,Bz};
            expected._A()[i,j] = vec3{0, Bz * i, Bx * j} * dx;
            expected._A()[i,j].z -= (i < nx/2 ? ByL : ByR) * (i - nx/2) * dx;
        }
    }
    
    auto grid = Grid3D(nx, ny, ny, dx, dx, dx);
    grid.boundary = Boundary::Outflow("X") + Boundary::Periodic("YZ");
    for(int i = 0; i<=nx; i++){
        for(int j = 0; j<=ny; j++){
            for(int k = 0; k<=ny; k++){
                grid[i,j,k] = expected[i,j];
                grid._A()[i,j,k] = expected._A()[i,j];
            }
        }
    }
    
    const double dt = 0.1*fmin(CFL::cfl_time(grid), CFL::cfl_time(expected));

    
    grid.advance(dt);
    expected.advance(dt);
    
    
    for (int i = 0; i < nx; i++){
        for (int j = 0; j < ny; j++){
            for(int k = 0; k<ny; k++){
                expect_close(grid[i,j,k], expected[i,j]);
            }
        }
    }
     */
    
}
//MARK: 2D - Y
void DRAGON_Test::verify_2D_Y_match_1D_MHD(){
    
}
//MARK: 3D - X
void DRAGON_Test::verify_3D_X_match_1D_MHD(){
    
}
//MARK: 3D - Y
void DRAGON_Test::verify_3D_Y_match_1D_MHD(){
    
}
//MARK: 3D - Z
void DRAGON_Test::verify_3D_Z_match_1D_MHD(){
    
}
