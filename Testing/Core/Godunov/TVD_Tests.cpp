//
//  TVD_Tests.cpp
//  DRAGON/Testing/Core/Godunov
//
//  Created by Bobbie Markwick on 18/06/2026.
//

#include "Testing.hpp"
#include "TVD.hpp"
#include "Config.h"
#include <iostream>

using namespace DRAGON_Test;

namespace {
    PrimitiveState make_tvd_state(double rho, double vx, double vy, double vz, double p) {
        PrimitiveState W = make_state(rho, vx, vy, vz, p);
#ifdef MHD
        W.B = {0.1 * rho, 0.2 * vx, 0.3 * vy};
#endif
        return W;
    }
}

void DRAGON_Test::verify_tvd(bool output) {
    if(output) std::cout << "TVD: \n";
    if(output) std::cout << "- Minmod: ";
    verify_tvd_scalar_minmod();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "- MC: ";
    verify_tvd_scalar_mc();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "- van Leer: ";
    verify_tvd_scalar_vanleer();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "- Superbee: ";
    verify_tvd_scalar_superbee();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "- van Albada: ";
    verify_tvd_scalar_vanalbada();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "- Vector limiters: ";
    verify_tvd_vec3_limiters();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "- Primitive limiters: ";
    verify_tvd_primitive_limiters();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "- Runtime dispatch: ";
    verify_tvd_limiter_dispatch();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "- MUSCL: ";
    verify_tvd_muscl_constant_state();
    verify_tvd_muscl_zero_dt_spatial_reconstruction();
    verify_tvd_muscl_nonzero_dt_predictor();
    verify_tvd_muscl_falls_back_when_unphysical();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "All TVD Tests Passed.\n\n";
}

//MARK: Scalar Limiters
void DRAGON_Test::verify_tvd_scalar_minmod() {
    assert(approx(TVD::minmod(2.0, 3.0), 2.0));
    assert(approx(TVD::minmod(3.0, 2.0), 2.0));
    assert(approx(TVD::minmod(-2.0, -3.0), -2.0));
    assert(approx(TVD::minmod(-2.0, 3.0), 0.0));
}

void DRAGON_Test::verify_tvd_scalar_mc() {
    assert(approx(TVD::MC(1.0, 3.0), 2.0));
    assert(approx(TVD::MC(3.0, 1.0), 2.0));
    assert(approx(TVD::MC(-1.0, -3.0), -2.0));
    assert(approx(TVD::MC(1.0, -3.0), 0.0));
}

void DRAGON_Test::verify_tvd_scalar_vanleer() {
    assert(approx(TVD::vanLeer(2.0, 4.0), 8.0 / 3.0));
    assert(approx(TVD::vanLeer(-2.0, -4.0), -8.0 / 3.0));
    assert(approx(TVD::vanLeer(2.0, -4.0), 0.0));
}

void DRAGON_Test::verify_tvd_scalar_superbee() {
    assert(approx(TVD::superbee(1.0, 4.0), 2.0));
    assert(approx(TVD::superbee(4.0, 1.0), 2.0));
    assert(approx(TVD::superbee(-1.0, -4.0), -2.0));
    assert(approx(TVD::superbee(1.0, -4.0), 0.0));
}

void DRAGON_Test::verify_tvd_scalar_vanalbada() {
    assert(approx(TVD::vanAlbada(2.0, 4.0), 12.0 / 5.0));
    assert(approx(TVD::vanAlbada(-2.0, -4.0), -12.0 / 5.0));
    assert(approx(TVD::vanAlbada(2.0, -4.0), 0.0));
}

//MARK: Vector/State Limiters
void DRAGON_Test::verify_tvd_vec3_limiters() {
    vec3 a = {2.0, -2.0, 2.0};
    vec3 b = {4.0, -4.0, -4.0};

    expect_close(TVD::minmod(a, b), {2.0, -2.0, 0.0});
    expect_close(TVD::MC(a, b), {3.0, -3.0, 0.0});
    expect_close(TVD::vanLeer(a, b), {8.0 / 3.0, -8.0 / 3.0, 0.0});
    expect_close(TVD::superbee(a, b), {4.0, -4.0, 0.0});
    expect_close(TVD::vanAlbada(a, b), {12.0 / 5.0, -12.0 / 5.0, 0.0});
}

void DRAGON_Test::verify_tvd_primitive_limiters() {
    PrimitiveState a = make_tvd_state(2.0, 4.0, -2.0, 3.0, 6.0);
    PrimitiveState b = make_tvd_state(4.0, 2.0, -5.0, -3.0, 3.0);

    PrimitiveState got = TVD::minmod(a, b);
    expect_close(got, make_tvd_state(2.0, 2.0, -2.0, 0.0, 3.0));

    got = TVD::MC(a, b);
    expect_close(got, make_tvd_state(3.0, 3.0, -3.5, 0.0, 4.5));

    got = TVD::vanLeer(a, b);
    expect_close(got, make_tvd_state(8.0 / 3.0, 8.0 / 3.0, -20.0 / 7.0, 0.0, 4.0));

    got = TVD::superbee(a, b);
    expect_close(got, make_tvd_state(4.0, 4.0, -4.0, 0.0, 6.0));

    got = TVD::vanAlbada(a, b);
    expect_close(got, make_tvd_state(2.4, 2.4, -70.0 / 29.0, 0.0, 3.6));
}

//MARK: Dispatch
void DRAGON_Test::verify_tvd_limiter_dispatch() {
    PrimitiveState a = make_tvd_state(2.0, 4.0, -2.0, 3.0, 6.0);
    PrimitiveState b = make_tvd_state(4.0, 2.0, -5.0, -3.0, 3.0);
    
    int prev = CONFIG::limiter_choice;

    CONFIG::limiter_choice = LIMITER_MINMOD;
    expect_close(TVD::limit(a, b), TVD::minmod(a, b));
    CONFIG::limiter_choice = LIMITER_MC;
    expect_close(TVD::limit(a, b), TVD::MC(a, b));
    CONFIG::limiter_choice = LIMITER_VANLEER;
    expect_close(TVD::limit(a, b), TVD::vanLeer(a, b));
    CONFIG::limiter_choice = LIMITER_SUPERBEE;
    expect_close(TVD::limit(a, b), TVD::superbee(a, b));
    CONFIG::limiter_choice = LIMITER_VANALBADA;
    expect_close(TVD::limit(a, b), TVD::vanAlbada(a, b));
    
    CONFIG::limiter_choice = prev; //Be a good citizen
}

//MARK: MUSCL
void DRAGON_Test::verify_tvd_muscl_constant_state() {
    int prev = CONFIG::limiter_choice;
    CONFIG::limiter_choice = LIMITER_MINMOD;
    PrimitiveState W = make_tvd_state(1.0, 0.2, -0.3, 0.4, 2.0);
    PrimitiveState L, R;

    TVD::MUSCL(W, L, W, R, W, 0.25);

    expect_close(L, W);
    expect_close(R, W);
    
    CONFIG::limiter_choice = prev; //Be a good citizen
}

void DRAGON_Test::verify_tvd_muscl_zero_dt_spatial_reconstruction() {
    int prev = CONFIG::limiter_choice;

    CONFIG::limiter_choice = LIMITER_MINMOD;
    PrimitiveState wL = make_state(1.0, 1.0, 2.0, 3.0, 4.0);
    PrimitiveState wC = make_state(2.0, 2.0, 3.0, 4.0, 5.0);
    PrimitiveState wR = make_state(3.0, 3.0, 4.0, 5.0, 6.0);
    PrimitiveState L, R;

    TVD::MUSCL(wL, L, wC, R, wR, 0.0);

    expect_close(L, make_state(1.5, 1.5, 2.5, 3.5, 4.5));
    expect_close(R, make_state(2.5, 2.5, 3.5, 4.5, 5.5));
    
    CONFIG::limiter_choice = prev; //Be a good citizen
}

void DRAGON_Test::verify_tvd_muscl_nonzero_dt_predictor() {
    int prev = CONFIG::limiter_choice;

    CONFIG::limiter_choice = LIMITER_MINMOD;
    PrimitiveState wL = make_state(1.0, 0.0, 0.0, 0.0, 1.0);
    PrimitiveState wC = make_state(2.0, 0.0, 0.0, 0.0, 2.0);
    PrimitiveState wR = make_state(3.0, 0.0, 0.0, 0.0, 3.0);
    PrimitiveState L, R;
    double dt_dL = 0.1;

    TVD::MUSCL(wL, L, wC, R, wR, dt_dL);

    PrimitiveState spatialL = make_state(1.5, 0.0, 0.0, 0.0, 1.5);
    PrimitiveState spatialR = make_state(2.5, 0.0, 0.0, 0.0, 2.5);
    ConservativeState UL(spatialL), UR(spatialR);
    ConservativeState correction = (UR.flux() - UL.flux()) * (0.5 * dt_dL);
    PrimitiveState expectedL = UL - correction;
    PrimitiveState expectedR = UR - correction;

    expect_close(L, expectedL, 1e-12, 1e-12);
    expect_close(R, expectedR, 1e-12, 1e-12);
    
    CONFIG::limiter_choice = prev; //Be a good citizen
}

void DRAGON_Test::verify_tvd_muscl_falls_back_when_unphysical() {
    int prev = CONFIG::limiter_choice;

    CONFIG::limiter_choice = LIMITER_SUPERBEE;
    PrimitiveState wL = make_state(1.0, 0.0, 0.0, 0.0, 1.0);
    PrimitiveState wC = make_state(1.0, 0.0, 0.0, 0.0, 1.0);
    PrimitiveState wR = make_state(10.0, 0.0, 0.0, 0.0, 1.0);
    PrimitiveState L, R;

    TVD::MUSCL(wL, L, wC, R, wR, 0.0);

    expect_close(L, wC);
    expect_close(R, wC);
    
    CONFIG::limiter_choice = prev; //Be a good citizen

}
