//
//  Core/Riemann/Roe_Tests.cpp
//  DRAGON/Testing
//
//  Created by Bobbie Markwick on 12/05/2026.
//

#include "Testing.hpp"
#include "Riemann.hpp"
#include <iostream>

using namespace DRAGON_Test;

namespace RiemannTestHelpers{
PrimitiveState mirrored_state(PrimitiveState W);
ConservativeState mirrored_flux(ConservativeState F);
}
using namespace RiemannTestHelpers;


//MARK: Test Problems
void DRAGON_Test::verify_roe_equal_state() {
    PrimitiveState W = make_state(1.0, 0.75, 0.2, -0.1, 1.0);

    Riemann problem(W, W);
    ConservativeState expected = ConservativeState(W).flux(W.v);

    expect_close(problem.Roe(), expected, 1e-10, 1e-10);
}
void DRAGON_Test::verify_roe_stationary_contact() {
    PrimitiveState L = make_state(1.0, 0.0, 0.0, 0.0, 1.0);
    PrimitiveState R = make_state(2.0, 0.0, 0.0, 0.0, 1.0);
    ConservativeState expected;
    expected.rho = 0.0;
    expected.p.x  = 1.0;
    expected.p.y  = 0.0;
    expected.p.z  = 0.0;
    expected.E   = 0.0;
    
    expect_close(Riemann(L,R).Roe(), expected, 1e-10, 1e-10);
}

void DRAGON_Test::verify_roe_supersonic_upwind() {
    {
        PrimitiveState L = make_state(1.0,   10.0, 0.0, 0.0, 1.0);
        PrimitiveState R = make_state(0.125, 10.0, 0.0, 0.0, 0.1);
        ConservativeState expected = ConservativeState(L).flux(L.v);
        expect_close(Riemann(L,R).Roe(), expected, 1e-10, 1e-10);
    }
    {
        PrimitiveState L = make_state(1.0,   -10.0, 0.0, 0.0, 1.0);
        PrimitiveState R = make_state(0.125, -10.0, 0.0, 0.0, 0.1);
        ConservativeState expected = ConservativeState(R).flux(R.v);
        expect_close(Riemann(L,R).Roe(), expected, 1e-10, 1e-10);
    }
}
void DRAGON_Test::verify_roe_supersonic_upwind_transverse() {
    {
        PrimitiveState L = make_state(1.0, 10.0, 2.0, -3.0, 1.0);
        PrimitiveState R = make_state(0.125, 10.0, -8.0, 9.0, 0.1);
        ConservativeState expected = ConservativeState(L).flux(L.v);
        expect_close(Riemann(L,R).Roe(),          expected, 1e-10, 1e-10);
    }
    {
        PrimitiveState L = make_state(1.0, -10.0, 2.0, -3.0, 1.0);
        PrimitiveState R = make_state(0.125, -10.0, -8.0, 9.0, 0.1);
        ConservativeState expected = ConservativeState(R).flux(R.v);
        expect_close(Riemann(L,R).Roe(),          expected, 1e-10, 1e-10);
    }
}
void DRAGON_Test::verify_roe_symmetry() {
    PrimitiveState L = make_state(1.0, 0.35, 0.2, -0.1, 1.0);
    PrimitiveState R = make_state(0.7, -0.45, -0.3, 0.4, 0.8);

    ConservativeState F = Riemann(L, R).Roe();
    ConservativeState reflected = Riemann(mirrored_state(R), mirrored_state(L)).Roe();
    expect_close(reflected, mirrored_flux(F));
}

void DRAGON_Test::verify_roe_entropy_fix_rarefactions() {
#ifdef Harten_Hyman
    {
        PrimitiveState L = make_state(1.0, -2.0, 0.2, -0.1, 0.4);
        PrimitiveState R = make_state(1.0,  0.0, -0.3, 0.4, 0.4);
        ConservativeState F = Riemann(L, R).Roe();
        expect_finite(F);
        assert((ConservativeState(L) - F*0.01).isPhysical());
        assert((ConservativeState(R) + F*0.01).isPhysical());
    }

    {
        PrimitiveState L = make_state(1.0,  0.0, 0.2, -0.1, 0.4);
        PrimitiveState R = make_state(1.0,  2.0, -0.3, 0.4, 0.4);
        ConservativeState F = Riemann(L, R).Roe();
        expect_finite(F);
        assert((ConservativeState(L) - F*0.01).isPhysical());
        assert((ConservativeState(R) + F*0.01).isPhysical());
    }
#endif
}
