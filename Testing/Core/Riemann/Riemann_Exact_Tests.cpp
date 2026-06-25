//
//  Riemann_Exact_Tests.cpp
//  DRAGON/Testing/Core/Riemann
//
//  Created by Bobbie Markwick on 12/05/2026.
//

#include "Testing.hpp"
#include "Riemann.hpp"
#include <iostream>

using namespace DRAGON_Test;


//MARK: Helpers

namespace RiemannTestHelpers{
PrimitiveState mirrored_state(PrimitiveState W) {
    W.v.x *= -1.0;
#ifdef MHD
    W.B.x *= -1.0;
#endif
    return W;
}

ConservativeState mirrored_flux(ConservativeState F) {
    F.rho *= -1.0;
    F.p.y *= -1.0;
    F.p.z *= -1.0;
    F.E *= -1.0;
#ifdef MHD
    F.B.y *= -1.0;
    F.B.z *= -1.0;
#endif
    return F;
}
}
using namespace RiemannTestHelpers;


//MARK: Solver components
void DRAGON_Test::verify_riemann_f() {
    PrimitiveState W = make_state(1.0, 0.0, 0.0, 0.0, 1.0);
    assert(approx(Riemann::f(W.p, W), 0.0));
    
    double f_shock = Riemann::f(2.0, W);
    double f_rarefaction = Riemann::f(0.5, W);
    
    assert(f_shock > 0.0);
    assert(f_rarefaction < 0.0);
    assert(std::isfinite(f_shock));
    assert(std::isfinite(f_rarefaction));
}


void DRAGON_Test::verify_sample_mirror_restores_state() {
    PrimitiveState L = make_state(1.0,   0.0, 0.0, 0.0, 1.0);
    PrimitiveState R = make_state(0.125, 0.0, 0.0, 0.0, 0.1);

    RiemannSolution S = Riemann(L, R).exact();

    PrimitiveState A = S.sample(-0.5);
    PrimitiveState B = S.sample(-0.5);
    PrimitiveState C = S.sample( 0.5);
    PrimitiveState D = S.sample( 0.5);

    expect_close(A, B, 1e-10, 1e-10);
    expect_close(C, D, 1e-10, 1e-10);

    expect_close(S.wL, L);
    expect_close(S.wR, R);
}

//MARK: Test Problems
void DRAGON_Test::verify_exact_equal_state() {
    PrimitiveState W = make_state(1.0, 0.75, 0.2, -0.1, 1.0);
    ConservativeState expected = ConservativeState(W).flux(W.v);

    RiemannSolution S = Riemann(W, W).exact();
    expect_close(S.wL, W);
    expect_close(S.wR, W);
    expect_close(S.sL, W, 1e-10, 1e-10);
    expect_close(S.sR, W, 1e-10, 1e-10);
    expect_close(S.flux(), expected, 1e-10, 1e-10);
}
void DRAGON_Test::verify_exact_stationary_contact() {
    PrimitiveState L = make_state(1.0, 0.0, 0.0, 0.0, 1.0);
    PrimitiveState R = make_state(2.0, 0.0, 0.0, 0.0, 1.0);
    ConservativeState expected;
    expected.rho = 0.0;
    expected.p.x  = 1.0;
    expected.p.y  = 0.0;
    expected.p.z  = 0.0;
    expected.E   = 0.0;

    expect_close(Riemann(L,R).exact().flux(), expected, 1e-10, 1e-10);
}

void DRAGON_Test::verify_exact_supersonic_upwind() {
    {
        PrimitiveState L = make_state(1.0,   10.0, 0.0, 0.0, 1.0);
        PrimitiveState R = make_state(0.125, 10.0, 0.0, 0.0, 0.1);
        ConservativeState expected = ConservativeState(L).flux(L.v);
        expect_close(Riemann(L,R).exact().flux(), expected, 1e-10, 1e-10);
    }
    {
        PrimitiveState L = make_state(1.0,   -10.0, 0.0, 0.0, 1.0);
        PrimitiveState R = make_state(0.125, -10.0, 0.0, 0.0, 0.1);
        ConservativeState expected = ConservativeState(R).flux(R.v);
        expect_close(Riemann(L,R).exact().flux(), expected, 1e-10, 1e-10);
    }
}

void DRAGON_Test::verify_exact_sod() {
    PrimitiveState L = make_state(1.0,   0.0, 0.0, 0.0, 1.0);
    PrimitiveState R = make_state(0.125, 0.0, 0.0, 0.0, 0.1);

    RiemannSolution S = Riemann(L, R).exact();

    assert(approx(S.sL.p,   0.29394518766601785, 1e-10, 1e-10));
    assert(approx(S.sR.p,   0.29394518766601785, 1e-10, 1e-10));
    assert(approx(S.sL.v.x,  0.8411948521688083, 1e-10, 1e-10));
    assert(approx(S.sR.v.x,  0.8411948521688083, 1e-10, 1e-10));
    assert(approx(S.sL.rho, 0.4796890587209175, 1e-10, 1e-10));
    assert(approx(S.sR.rho, 0.22980574931194703, 1e-10, 1e-10));
    

    assert(approx(S.sL.v.y, L.v.y));
    assert(approx(S.sL.v.z, L.v.z));
    assert(approx(S.sR.v.y, R.v.y));
    assert(approx(S.sR.v.z, R.v.z));
    
    PrimitiveState W0 = S.sample(0.0);

    assert(approx(W0.rho, S.sL.rho, 1e-10, 1e-10));
    assert(approx(W0.v.x,  S.sL.v.x,  1e-10, 1e-10));
    assert(approx(W0.p,   S.sL.p,   1e-10, 1e-10));

    expect_close(S.sample(-100.0), L, 1e-10, 1e-10);
    expect_close(S.sample( 100.0), R, 1e-10, 1e-10);

    PrimitiveState fan = S.sample(-0.5);
    assert(fan.rho < L.rho && fan.rho > S.sL.rho);
    assert(fan.p   < L.p   && fan.p   > S.sL.p);
    assert(fan.v.x > L.v.x && fan.v.x < S.sL.v.x);
    expect_finite(fan);
}

void DRAGON_Test::verify_exact_supersonic_upwind_transverse() {
    {
        PrimitiveState L = make_state(1.0, 10.0, 2.0, -3.0, 1.0);
        PrimitiveState R = make_state(0.125, 10.0, -8.0, 9.0, 0.1);
        ConservativeState expected = ConservativeState(L).flux(L.v);
        expect_close(Riemann(L,R).exact().flux(), expected, 1e-10, 1e-10);
    }

    {
        PrimitiveState L = make_state(1.0, -10.0, 2.0, -3.0, 1.0);
        PrimitiveState R = make_state(0.125, -10.0, -8.0, 9.0, 0.1);
        ConservativeState expected = ConservativeState(R).flux(R.v);
        expect_close(Riemann(L,R).exact().flux(), expected, 1e-10, 1e-10);
    }
}

void DRAGON_Test::verify_riemann_exact_symmetry() {
    PrimitiveState L = make_state(1.0, 0.35, 0.2, -0.1, 1.0);
    PrimitiveState R = make_state(0.7, -0.45, -0.3, 0.4, 0.8);

    ConservativeState F = Riemann(L, R).exact().flux();
    ConservativeState reflected = Riemann(mirrored_state(R), mirrored_state(L)).exact().flux();

    expect_close(reflected, mirrored_flux(F));
}
