//
//  DimMatch_Tests.cpp
//  DRAGON/Testing/Core/Godunov
//
//  Created by Bobbie Markwick on 09/06/2026.
//

#include "Testing.hpp"
#include "Grid.hpp"

#include <iostream> //For std::cout
#include "CFL.hpp"  //These tests only work on a single timestep

#ifdef MHD

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

void DRAGON_Test::verify_godunov_2D_MHD_Match(bool output){
    int prev = CONFIG::riemann_choice;
    CONFIG::riemann_choice = RIEMANN_HLLD;

    if(output) std::cout<<"- 1D Match (MHD-X): ";
    verify_2D_X_match_1D_MHD();
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- 1D Match (MHD-Y): ";
    verify_2D_Y_match_1D_MHD();
    if(output) std::cout<<"Passed\n";
    
    CONFIG::riemann_choice = prev;
}


void DRAGON_Test::verify_godunov_3D_MHD_Match(bool output){
    int prev = CONFIG::riemann_choice;
    CONFIG::riemann_choice = RIEMANN_HLLD;

    if(output) std::cout<<"- 1D Match (MHD-X): ";
    verify_3D_X_match_1D_MHD();
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- 1D Match (MHD-Y): ";
    verify_3D_Y_match_1D_MHD();
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- 1D Match (MHD-Z): ";
    verify_3D_Z_match_1D_MHD();
    if(output) std::cout<<"Passed\n";
    
    CONFIG::riemann_choice = prev;
}


//MARK: 2D - X
void DRAGON_Test::verify_2D_X_match_1D_MHD(){
    auto expected = Grid1D(nx, dx);
    expected.boundary =  Boundary::Outflow("X") + Boundary::Periodic("YZ");

    for(int i = 0; i<nx; i++){
        expected[i].rho = i < nx/2 ? rho_L : rho_R;
        expected[i].p = i < nx/2 ? p_L : p_R;
        expected[i].v = i < nx/2 ? vec3{vx,vy,vz} : vec3{0,0,0};
        expected[i].B = i < nx/2 ? vec3{Bx,ByL,Bz} : vec3{Bx,ByR,Bz};
    }
    
    auto grid = Grid2D(nx, ny, dx, dx);
    grid.boundary = Boundary::Outflow("X") + Boundary::Periodic("YZ");
    for(int i = 0; i<=nx; i++){
        for(int j = 0; j<=ny; j++){
            grid[i,j] = expected[i];
            grid._B()[i,j] = i < nx/2 ? vec3{Bx,ByL,Bz} : vec3{Bx,ByR,Bz};
        }
    }
    
    const double dt = std::min(CFL::cfl_time(grid), CFL::cfl_time(expected));
    grid.advance(dt);
    expected.advance(dt);
    
    
    for (int i = 0; i < nx; i++){
        for (int j = 0; j < ny; j++){
            assert(approx(grid[i,j].rho, expected[i].rho));
            assert(approx(grid[i,j].energy(), expected[i].energy()));
            assert(approx(grid[i,j].B.x, expected[i].B.x));
            expect_close(grid[i,j].v, expected[i].v);
        }
    }
}
//MARK: 2D - Y
void DRAGON_Test::verify_2D_Y_match_1D_MHD(){
    auto expected = Grid1D(nx, dx);
    expected.boundary =  Boundary::Periodic("X") + Boundary::Periodic("XZ");

    for(int i = 0; i<nx; i++){
        expected[i].rho = i < nx/2 ? rho_L : rho_R;
        expected[i].p = i < nx/2 ? p_L : p_R;
        expected[i].v = i < nx/2 ? vec3{vx,vy,vz} : vec3{0,0,0};
        expected[i].B = i < nx/2 ? vec3{Bx,ByL,Bz} : vec3{Bx,ByR,Bz};
    }
    
    auto grid = Grid2D(ny, nx, dx, dx);
    grid.boundary = Boundary::Periodic("Y") + Boundary::Periodic("XZ");
    for(int i = 0; i<=ny; i++){
        for(int j = 0; j<=nx; j++){
            grid[i,j] = expected[j].swappedXY();
            grid._B()[i,j] = j < nx/2 ? vec3{ByL,Bx,Bz} : vec3{ByR,Bx,Bz};
        }
    }
    
    const double dt = std::min(CFL::cfl_time(grid), CFL::cfl_time(expected));
    grid.advance(dt);
    expected.advance(dt);
    
    
    for (int i = 0; i < ny; i++){
        for (int j = 0; j < nx; j++){
            assert(approx(grid[i,j].rho, expected[j].rho));
            assert(approx(grid[i,j].energy(), expected[j].energy()));
            assert(approx(grid[i,j].B.y, expected[j].B.x));
            expect_close(grid[i,j].v.swappedXY(), expected[j].v);
        }
    }
}
//MARK: 3D - X
void DRAGON_Test::verify_3D_X_match_1D_MHD(){
    auto expected = Grid1D(nx, dx);
    expected.boundary =  Boundary::Outflow("X") + Boundary::Periodic("YZ");

    for(int i = 0; i<nx; i++){
        expected[i].rho = i < nx/2 ? rho_L : rho_R;
        expected[i].p = i < nx/2 ? p_L : p_R;
        expected[i].v = i < nx/2 ? vec3{vx,vy,vz} : vec3{0,0,0};
        expected[i].B = i < nx/2 ? vec3{Bx,ByL,Bz} : vec3{Bx,ByR,Bz};
    }
    
    auto grid = Grid3D(nx, ny, ny, dx, dx, dx);
    grid.boundary = Boundary::Outflow("X") + Boundary::Periodic("YZ");
    for(int i = 0; i<=nx; i++){
        for(int j = 0; j<=ny; j++){
            for(int k = 0; k<=ny; k++){
                grid[i,j,k] = expected[i];
                grid._B()[i,j,k] = i < nx/2 ? vec3{Bx,ByL,Bz} : vec3{Bx,ByR,Bz};
            }
        }
    }
    
    const double dt = std::min(CFL::cfl_time(grid), CFL::cfl_time(expected));
    grid.advance(dt);
    expected.advance(dt);
    
    
    for (int i = 0; i < nx; i++){
        for (int j = 0; j < ny; j++){
            for(int k = 0; k<ny; k++){
                assert(approx(grid[i,j,k].rho, expected[i].rho));
                assert(approx(grid[i,j,k].energy(), expected[i].energy()));
                assert(approx(grid[i,j,k].B.x, expected[i].B.x));
                expect_close(grid[i,j,k].v, expected[i].v);
            }
        }
    }
}
//MARK: 3D - Y
void DRAGON_Test::verify_3D_Y_match_1D_MHD(){
    auto expected = Grid1D(nx, dx);
    expected.boundary =  Boundary::Outflow("X") + Boundary::Periodic("YZ");

    for(int i = 0; i<nx; i++){
        expected[i].rho = i < nx/2 ? rho_L : rho_R;
        expected[i].p = i < nx/2 ? p_L : p_R;
        expected[i].v = i < nx/2 ? vec3{vx,vy,vz} : vec3{0,0,0};
        expected[i].B = i < nx/2 ? vec3{Bx,ByL,Bz} : vec3{Bx,ByR,Bz};
    }
    
    auto grid = Grid3D(ny, nx, ny, dx, dx, dx);
    grid.boundary = Boundary::Outflow("Y") + Boundary::Periodic("XZ");
    for(int i = 0; i<=ny; i++){
        for(int j = 0; j<=nx; j++){
            for(int k = 0; k<=ny; k++){
                grid[i,j,k] = expected[j].swappedXY();
                grid._B()[i,j,k] = j < nx/2 ? vec3{ByL,Bx,Bz} : vec3{ByR,Bx,Bz};
            }
        }
    }
    
    const double dt = std::min(CFL::cfl_time(grid), CFL::cfl_time(expected));
    grid.advance(dt);
    expected.advance(dt);
    
    
    for (int i = 0; i < ny; i++){
        for (int j = 0; j < nx; j++){
            for(int k = 0; k< ny; k++){
                assert(approx(grid[i,j,k].rho, expected[j].rho));
                assert(approx(grid[i,j,k].energy(), expected[j].energy()));
                expect_close(grid[i,j,k].v.swappedXY(), expected[j].v);
                assert(approx(grid[i,j,k].B.y, expected[j].B.x));
            }
        }
    }
}
//MARK: 3D - Z
void DRAGON_Test::verify_3D_Z_match_1D_MHD(){
    auto expected = Grid1D(nx, dx);
    expected.boundary =  Boundary::Outflow("X") + Boundary::Periodic("YZ");

    for(int i = 0; i<nx; i++){
        expected[i].rho = i < nx/2 ? rho_L : rho_R;
        expected[i].p = i < nx/2 ? p_L : p_R;
        expected[i].v = i < nx/2 ? vec3{vx,vy,vz} : vec3{0,0,0};
        expected[i].B = i < nx/2 ? vec3{Bx,ByL,Bz} : vec3{Bx,ByR,Bz};
    }
    
    auto grid = Grid3D(ny, ny, nx, dx, dx, dx);
    grid.boundary = Boundary::Outflow("Z") + Boundary::Periodic("XY");
    for(int i = 0; i<=ny; i++){
        for(int j = 0; j<=ny; j++){
            for(int k = 0; k<=nx; k++){
                grid[i,j,k] = expected[k].swappedXZ();
                grid._B()[i,j,k] = k < nx/2 ? vec3{Bz,ByL,Bx} : vec3{Bz,ByR,Bx};
            }
        }
    }
    
    const double dt = std::min(CFL::cfl_time(grid), CFL::cfl_time(expected));
    grid.advance(dt);
    expected.advance(dt);
    
    
    for (int i = 0; i < ny; i++){
        for (int j = 0; j < ny; j++){
            for(int k = 0; k < nx; k++){
                assert(approx(grid[i,j,k].rho, expected[k].rho));
                assert(approx(grid[i,j,k].energy(), expected[k].energy()));
                expect_close(grid[i,j,k].v.swappedXZ(), expected[k].v);
                assert(approx(grid[i,j,k].B.z, expected[k].B.x));
            }
        }
    }
}
#endif
