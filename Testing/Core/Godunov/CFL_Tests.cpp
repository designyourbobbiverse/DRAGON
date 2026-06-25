//
//  Core/Godunov/CFL_Tests.cpp
//  DRAGON/Testing
//
//  Created by Bobbie Markwick on 09/06/2026.
//

#include "Testing.hpp"
#include "CFL.hpp"
#include "Config.h"
#include "Constants.h"
#include <math.h>
#include <iostream>

using namespace DRAGON_Test;

void DRAGON_Test::verify_cfl(bool output){
    std::cout << "CFL: \n";
    //Algorithm Choices
    std::cout << "- Individual Cell: ";
    verify_cfl_add_speed_1D();
    verify_cfl_add_speed_2D();
    verify_cfl_add_speed_3D();
    //std::cout << "Passed\n";
    //std::cout << "- Max Speed Variant: ";
    verify_cfl_max_speed_1D();
    verify_cfl_max_speed_2D();
    verify_cfl_max_speed_3D();
    //std::cout << "Passed\n";
    //std::cout << "- Power-Normed Speed: ";
    verify_cfl_pow_speed_1D();
    verify_cfl_pow_speed_2D();
    verify_cfl_pow_speed_3D();
    //std::cout << "Passed\n";
    //std::cout << "- Dispatch: ";
    verify_cfl_dispatch_add_3D();
    verify_cfl_dispatch_max_3D();
    verify_cfl_dispatch_pow_3D();
    verify_cfl_mhd_speed_3D();
    std::cout << "Passed\n";
    //Grid
    std::cout << "- 1D Grid: ";
    verify_cfl_time_1D_uniform();
    verify_cfl_time_1D_uses_fastest_cell();
    verify_cfl_time_1D_ignores_ghost_cells();
    std::cout << "Passed\n";
    std::cout << "- 2D Grid: ";
    verify_cfl_time_2D_uniform();
    verify_cfl_time_2D_visits_last_cell();
    verify_cfl_time_2D_ignores_ghost_cells();
    std::cout << "Passed\n";
    std::cout << "- 3D Grid: ";
    verify_cfl_time_3D_uniform();
    verify_cfl_time_3D_uses_fastest_cell();
    verify_cfl_time_3D_ignores_ghost_cells();
    verify_cfl_mhd_time_3D();
    std::cout << "Passed\n";
    std::cout << "All CFL Tests Passed.\n\n";

}


//MARK: CFL Speed = Added Speeds
void DRAGON_Test::verify_cfl_add_speed_1D() {
    PrimitiveState W = make_state(1.0, -3.0, 0.0, 0.0, 1.0);
    double dx = 5.0;
    double a = sqrt(_gamma * W.p / W.rho);
    
    double expected = (std::fabs(W.v.x) + a) / dx;
    double got = CFL::cfl_add_speed(W, dx);
    assert(approx(got, expected));
}

void DRAGON_Test::verify_cfl_add_speed_2D() {
    PrimitiveState W = make_state(2.0, 3.0, -4.0, 0.0, 5.0);
    double dx = 0.25, dy = 0.10;
    double a = sqrt(_gamma * W.p / W.rho);

    double sx = (std::fabs(W.v.x) + a) / dx;
    double sy = (std::fabs(W.v.y) + a) / dy;
    double got = CFL::cfl_add_speed(W, dx, dy);
    assert(approx(got, sx + sy));
}

void DRAGON_Test::verify_cfl_add_speed_3D() {
    PrimitiveState W = make_state(2.0, 3.0, -4.0, 5.0, 5.0);
    double dx = 0.25, dy = 0.10, dz  = 0.50;
    double a = sqrt(_gamma * W.p / W.rho);

    double sx = (std::fabs(W.v.x) + a) / dx;
    double sy = (std::fabs(W.v.y) + a) / dy;
    double sz = (std::fabs(W.v.z) + a) / dz;
    double expected = sx + sy + sz;
    double got = CFL::cfl_add_speed(W, dx, dy, dz);

    assert(approx(got, expected));
}
//MARK: CFL Speed = Max Speed
void DRAGON_Test::verify_cfl_max_speed_1D() {
    PrimitiveState W = make_state(1.0, 3.0, 0.0, 0.0, 1.0);
    double dx = 0.25;
    double a = sqrt(_gamma * W.p / W.rho);
    
    double expected = (std::fabs(W.v.x) + a) / dx;
    double got = CFL::cfl_max_speed(W, dx);
    assert(approx(got, expected));
}

void DRAGON_Test::verify_cfl_max_speed_2D() {
    PrimitiveState W = make_state(2.0, 3.0, -4.0, 0.0, 5.0);
    double dx = 0.25, dy = 0.10;
    double a = sqrt(_gamma * W.p / W.rho);

    double sx = (std::fabs(W.v.x) + a) / dx;
    double sy = (std::fabs(W.v.y) + a) / dy;
    double got = CFL::cfl_max_speed(W, dx, dy);
    assert(approx(got, std::max(sx, sy)));
}

void DRAGON_Test::verify_cfl_max_speed_3D() {
    PrimitiveState W = make_state(2.0, 3.0, -4.0, 5.0, 5.0);
    double dx = 0.25, dy = 0.10, dz = 0.50;
    double a = sqrt(_gamma * W.p / W.rho);

    double sx = (std::fabs(W.v.x) + a) / dx;
    double sy = (std::fabs(W.v.y) + a) / dy;
    double sz = (std::fabs(W.v.z) + a) / dz;
    double expected = std::max(sx, std::max(sy, sz));
    double got = CFL::cfl_max_speed(W, dx, dy, dz);

    assert(approx(got, expected));
}

//MARK: CFL Speed = L2 Norm
void DRAGON_Test::verify_cfl_pow_speed_1D() {
    PrimitiveState W = make_state(1.0, -3.0, 0.0, 0.0, 1.0);
    double dx = 5.0;
    double a = sqrt(_gamma * W.p / W.rho);
    
    double expected = (std::fabs(W.v.x) + a) / dx;
    double got = CFL::cfl_pow_speed(W,2, dx);
    assert(approx(got, expected));
}

void DRAGON_Test::verify_cfl_pow_speed_2D() {
    PrimitiveState W = make_state(2.0, 3.0, -4.0, 0.0, 5.0);
    double dx = 0.25, dy = 0.10;
    double a = sqrt(_gamma * W.p / W.rho);

    double sx = (std::fabs(W.v.x) + a) / dx;
    double sy = (std::fabs(W.v.y) + a) / dy;
    double got = CFL::cfl_pow_speed(W,2, dx, dy);
    assert(approx(got, sqrt(sx*sx + sy*sy)));
}

void DRAGON_Test::verify_cfl_pow_speed_3D() {
    PrimitiveState W = make_state(2.0, 3.0, -4.0, 5.0, 5.0);
    double dx = 0.25, dy = 0.10, dz  = 0.50;
    double a = sqrt(_gamma * W.p / W.rho);

    double sx = (std::fabs(W.v.x) + a) / dx;
    double sy = (std::fabs(W.v.y) + a) / dy;
    double sz = (std::fabs(W.v.z) + a) / dz;
    double expected = sqrt(sx*sx + sy*sy + sz*sz);
    double got = CFL::cfl_pow_speed(W,2, dx, dy, dz);

    assert(approx(got, expected));
}
//MARK: CFL Dispatch

void DRAGON_Test::verify_cfl_dispatch_add_3D() {
    PrimitiveState W = make_state(2.0, 3.0, -4.0, 5.0, 5.0);
    double dx = 0.25, dy = 0.10, dz = 0.50;

    int previous = CONFIG::cfl_choice;
    CONFIG::cfl_choice = CFL_ADD;

    double expected = CFL::cfl_add_speed(W, dx, dy, dz);
    double got = CFL::cfl_speed(W, dx, dy, dz);
    assert(approx(got, expected));
    CONFIG::cfl_choice = previous;
}

void DRAGON_Test::verify_cfl_dispatch_max_3D() {
    PrimitiveState W = make_state(2.0, 3.0, -4.0, 5.0, 5.0);
    double dx = 0.25, dy = 0.10, dz = 0.50;

    int previous = CONFIG::cfl_choice;
    CONFIG::cfl_choice = CFL_MAX;

    double expected = CFL::cfl_max_speed(W, dx, dy, dz);
    double got = CFL::cfl_speed(W, dx, dy, dz);
    assert(approx(got, expected));
    CONFIG::cfl_choice = previous;
}

void DRAGON_Test::verify_cfl_dispatch_pow_3D() {
    PrimitiveState W = make_state(2.0, 3.0, -4.0, 5.0, 5.0);
    double dx = 0.25, dy = 0.10, dz = 0.50;

    int previous = CONFIG::cfl_choice;
    CONFIG::cfl_choice = 2;

    double expected = CFL::cfl_pow_speed(W, 2.0, dx, dy, dz);
    double got = CFL::cfl_speed(W, dx, dy, dz);
    assert(approx(got, expected));
    CONFIG::cfl_choice = previous;
}



//MARK: 1D Grid
void DRAGON_Test::verify_cfl_time_1D_uniform() {
    CONFIG::cfl_choice = CFL_ADD;
    Grid1D g(4, 0.5, 1);
    PrimitiveState W = make_state(1.0, 2.0, 0.0, 0.0, 1.0);
    for(int i = -1; i < 4 + 1; i++)  g[i] = W;

    double expected = CFL_coeff * g.dx / (2 + sqrt(_gamma));
    double got = CFL::cfl_time(g);
    assert(approx(got, expected));
}

void DRAGON_Test::verify_cfl_time_1D_uses_fastest_cell() {
    CONFIG::cfl_choice = CFL_ADD;
    Grid1D g(5, 0.25, 1);
    PrimitiveState slow = make_state(1.0, 0.1, 0.0, 0.0, 1.0);
    PrimitiveState fast = make_state(1.0, 9.0, 0.0, 0.0, 1.0);
    for(int i = -1; i < 5 + 1; i++) g[i] = slow;
    g[3] = fast;

    double expected = CFL_coeff * g.dx / (9 + sqrt(_gamma));
    double got = CFL::cfl_time(g);
    assert(approx(got, expected));
}

void DRAGON_Test::verify_cfl_time_1D_ignores_ghost_cells() {
    CONFIG::cfl_choice = CFL_ADD;
    Grid1D g(5, 0.25, 1);
    PrimitiveState slow = make_state(1.0, 0.1, 0.0, 0.0, 1.0);
    for(int i = -1; i < 5 + 1; i++) g[i] = slow;
    g[-1] = make_state(1.0, 1000.0, 0.0, 0.0, 1.0);
    g[5] = make_state(1.0, 1000.0, 0.0, 0.0, 1.0);

    double expected = CFL_coeff * g.dx / (0.1 + sqrt(_gamma));
    double got = CFL::cfl_time(g);
    assert(approx(got, expected));
}

//MARK: 2D Grid
void DRAGON_Test::verify_cfl_time_2D_uniform() {
    CONFIG::cfl_choice = CFL_ADD;
    Grid2D g(3,4, 0.5, 0.25, 1);
    PrimitiveState W = make_state(1.0, 2.0, 3.0, 4.0, 1.0);
    for(int i = -1; i < 3 + 1; i++){
        for(int j = -1; j < 4 + 1; j++){
            g[i,j] = W;
        }
    }

    double expected = CFL_coeff / ( (2 + sqrt(_gamma)) / g.dx  + (3 + sqrt(_gamma)) / g.dy );
    double got = CFL::cfl_time(g);
    assert(approx(got, expected));
}

void DRAGON_Test::verify_cfl_time_2D_visits_last_cell() {
    CONFIG::cfl_choice = CFL_ADD;
    PrimitiveState slow = make_state(1.0, 0.1, 0.0, 0.0, 1.0);
    PrimitiveState fast = make_state(1.0, 9.0, 0.0, 0.0, 1.0);
    {
        Grid2D g(3,4, 0.5,0.25, 1);
        for(int i = -1; i < 3 + 1; i++){
            for(int j = -1; j < 4 + 1; j++){
                g[i,j] = slow;
            }
        }
        g[2,3] = fast;
        
        double expected = CFL_coeff / ( (9 + sqrt(_gamma)) / g.dx  + (sqrt(_gamma)) / g.dy );
        double got = CFL::cfl_time(g);
        assert(approx(got, expected));
    }
    {
        Grid2D g(4,3, 0.5,0.25, 1);
        for(int i = -1; i < 4 + 1; i++){
            for(int j = -1; j < 3 + 1; j++){
                g[i,j] = slow;
            }
        }
        g[3,2] = fast;
        
        double expected = CFL_coeff / ( (9 + sqrt(_gamma)) / g.dx  + (sqrt(_gamma)) / g.dy );
        double got = CFL::cfl_time(g);
        assert(approx(got, expected));
    }
    
    
}

void DRAGON_Test::verify_cfl_time_2D_ignores_ghost_cells() {
    CONFIG::cfl_choice = CFL_ADD;
    Grid2D g(3,4, 0.5,0.25, 1);
    PrimitiveState slow = make_state(1.0, 0.1, 0.0, 0.0, 1.0);
    PrimitiveState fast = make_state(1.0, 100.0, 0.0, 0.0, 1.0);
    for(int i = -1; i < 3 + 1; i++){
        for(int j = -1; j < 4 + 1; j++){
            g[i,j] = (i<0 || i>=3 || j<0 || j>= 4) ? fast : slow;
        }
    }

    double expected = CFL_coeff / ( (0.1 + sqrt(_gamma)) / g.dx  + (sqrt(_gamma)) / g.dy );
    double got = CFL::cfl_time(g);
    assert(approx(got, expected));
}

//MARK: 3D Grid
void DRAGON_Test::verify_cfl_time_3D_uniform() {
    CONFIG::cfl_choice = CFL_ADD;
    Grid3D g(3,4,5, 0.5, 0.25,5.0, 1);
    PrimitiveState W = make_state(1.0, 2.0, 3.0, 4.0, 1.0);
    for(int i = -1; i < 3 + 1; i++){
        for(int j = -1; j < 4 + 1; j++){
            for(int k = -1; k < 5 + 1; k++){
                g[i,j,k] = W;
            }
        }
    }

    double expected = CFL_coeff / ( (2 + sqrt(_gamma)) / g.dx  + (3 + sqrt(_gamma)) / g.dy + (4 + sqrt(_gamma)) / g.dz);
    double got = CFL::cfl_time(g);
    assert(approx(got, expected));
}

void DRAGON_Test::verify_cfl_time_3D_uses_fastest_cell() {
    CONFIG::cfl_choice = CFL_ADD;
    Grid3D g(3,4,5, 0.5, 0.25,5.0, 1);
    PrimitiveState slow = make_state(1.0, 0.1, 0.0, 0.0, 1.0);
    PrimitiveState fast = make_state(1.0, 9.0, 0.0, 0.0, 1.0);
    for(int i = -1; i < 3 + 1; i++){
        for(int j = -1; j < 4 + 1; j++){
            for(int k = -1; k < 5 + 1; k++){
                g[i,j,k] = slow;
            }
        }
    }
    g[2,3,4] = fast;

    double expected = CFL_coeff / ( (9 + sqrt(_gamma)) / g.dx  + (sqrt(_gamma)) / g.dy + (sqrt(_gamma)) / g.dz );
    double got = CFL::cfl_time(g);
    assert(approx(got, expected));
}

void DRAGON_Test::verify_cfl_time_3D_ignores_ghost_cells() {
    CONFIG::cfl_choice = CFL_ADD;
    Grid3D g(3,4,5, 0.5, 0.25,5.0, 1);
    PrimitiveState slow = make_state(1.0, 0.1, 0.0, 0.0, 1.0);
    PrimitiveState fast = make_state(1.0, 100.0, 0.0, 0.0, 1.0);
    for(int i = -1; i < 3 + 1; i++){
        for(int j = -1; j < 4 + 1; j++){
            for(int k = -1; k < 5 + 1; k++){
                g[i,j,k] = (i<0 || i>=3 || j<0 || j>= 4 || k<0 || k>=5) ? fast : slow;
            }
        }
    }

    double expected = CFL_coeff / ( (0.1 + sqrt(_gamma)) / g.dx  + (sqrt(_gamma)) / g.dy + (sqrt(_gamma)) / g.dz);
    double got = CFL::cfl_time(g);
    assert(approx(got, expected));
}
//MARK: MHD
void DRAGON_Test::verify_cfl_mhd_speed_3D() {
#ifdef MHD
    PrimitiveState W = make_state(2.0, 3.0, -4.0, 5.0, 5.0);
    W.B = {2.0, -1.0, 3.0};
    double dx = 0.25, dy = 0.10, dz = 0.50;
    double a = W.c_fast();

    double sx = (std::fabs(W.v.x) + a) / dx;
    double sy = (std::fabs(W.v.y) + a) / dy;
    double sz = (std::fabs(W.v.z) + a) / dz;

    assert(approx(CFL::cfl_add_speed(W, dx, dy, dz), sx + sy + sz));
    assert(approx(CFL::cfl_max_speed(W, dx, dy, dz), std::max(sx, std::max(sy, sz))));
    assert(approx(CFL::cfl_pow_speed(W, 2, dx, dy, dz), sqrt(sx*sx + sy*sy + sz*sz)));
#endif
}

void DRAGON_Test::verify_cfl_mhd_time_3D() {
#ifdef MHD
    int previous = CONFIG::cfl_choice;
    CONFIG::cfl_choice = CFL_ADD;

    Grid3D g(3,4,5, 0.5, 0.25, 5.0, 1);
    PrimitiveState W = make_state(1.0, 2.0, 3.0, 4.0, 1.0);
    W.B = {1.0, -2.0, 0.5};
    for(int i = -1; i < 3 + 1; i++){
        for(int j = -1; j < 4 + 1; j++){
            for(int k = -1; k < 5 + 1; k++){
                g[i,j,k] = W;
            }
        }
    }

    double a = W.c_fast();
    double speed = (std::fabs(W.v.x) + a) / g.dx + (std::fabs(W.v.y) + a) / g.dy + (std::fabs(W.v.z) + a) / g.dz;
    assert(approx(CFL::cfl_time(g), CFL_coeff / speed));
    CONFIG::cfl_choice = previous;
#endif
}
