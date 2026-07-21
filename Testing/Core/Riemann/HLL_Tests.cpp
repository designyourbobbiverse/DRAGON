//
//  HLL_Tests.cpp
//  DRAGON/Testing/Core/Riemann
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


//MARK: Manual Wave Speeds
void DRAGON_Test::verify_hll_manual_wave_speeds() {
    PrimitiveState L = make_state(1.0,   0.0, 0.0, 0.0, 1.0);
    PrimitiveState R = make_state(0.125, 0.0, 0.0, 0.0, 0.1);

    Riemann problem(L, R);

    ConservativeState UL(L), UR(R);
    ConservativeState FL = UL.flux(L.v);
    ConservativeState FR = UR.flux(R.v);

    double SL = -1.0;
    double SR =  2.0;

    ConservativeState expected = (SR * FL - SL * FR + SL * SR * (UR - UL)) / (SR - SL);

    expect_close(problem.HLL(SL, SR), expected, 1e-12, 1e-12);

    expect_close(problem.HLL( 0.1, 2.0), FL);
    expect_close(problem.HLL(-2.0,-0.1), FR);
}
void DRAGON_Test::verify_hllc_manual_wave_speeds() {
    PrimitiveState L = make_state(1.0,   0.0, 0.0, 0.0, 1.0);
    PrimitiveState R = make_state(0.125, 0.0, 0.0, 0.0, 0.1);

    Riemann problem(L, R);

    ConservativeState FL = ConservativeState(L).flux(L.v);
    ConservativeState FR = ConservativeState(R).flux(R.v);

    expect_close(problem.HLLC( 0.1, 2.0), FL);
    expect_close(problem.HLLC(-2.0,-0.1), FR);

    expect_finite(problem.HLLC(-1.0, 2.0));
}


//MARK: Hydro Test Problems
void DRAGON_Test::verify_hll_equal_state() {
    PrimitiveState W = make_state(1.0, 0.75, 0.2, -0.1, 1.0);
    ConservativeState expected = ConservativeState(W).flux(W.v);

    expect_close(Riemann(W,W).HLL(), expected, 1e-10, 1e-10);
    expect_close(Riemann(W,W).HLLC(), expected, 1e-10, 1e-10);
#ifdef MHD
    expect_close(Riemann(W,W).HLLD(), expected, 1e-10, 1e-10);
#endif
    expect_close(Riemann(W,W).HLLE(), expected, 1e-10, 1e-10);
}
void DRAGON_Test::verify_hll_stationary_contact() {
    PrimitiveState L = make_state(1.0, 0.0, 0.0, 0.0, 1.0);
    PrimitiveState R = make_state(2.0, 0.0, 0.0, 0.0, 1.0);
    ConservativeState expected;
    expected.rho = 0.0;
    expected.mom.x  = 1.0;
    expected.mom.y  = 0.0;
    expected.mom.z  = 0.0;
    expected.E   = 0.0;

    expect_close(Riemann(L,R).HLLC(),expected, 1e-10, 1e-10);
#ifdef MHD
    expect_close(Riemann(L,R).HLLD(), expected, 1e-10, 1e-10);
#endif
    // 2-wave HLL may diffuse contacts in evolution, but this exact stationary case should still be finite.
    expect_finite(Riemann(L,R).HLL());
    expect_finite(Riemann(L,R).HLLE());
}

void DRAGON_Test::verify_hll_supersonic_upwind() {
    {
        PrimitiveState L = make_state(1.0,   10.0, 0.0, 0.0, 1.0);
        PrimitiveState R = make_state(0.125, 10.0, 0.0, 0.0, 0.1);
        ConservativeState expected = ConservativeState(L).flux(L.v);
        expect_close(Riemann(L,R).HLL(), expected, 1e-10, 1e-10);
        expect_close(Riemann(L,R).HLLC(), expected, 1e-10, 1e-10);
    #ifdef MHD
        expect_close(Riemann(L,R).HLLD(), expected, 1e-10, 1e-10);
    #endif
        expect_close(Riemann(L,R).HLLE(), expected, 1e-10, 1e-10);
    }
    {
        PrimitiveState L = make_state(1.0,   -10.0, 0.0, 0.0, 1.0);
        PrimitiveState R = make_state(0.125, -10.0, 0.0, 0.0, 0.1);
        ConservativeState expected = ConservativeState(R).flux(R.v);
        expect_close(Riemann(L,R).HLL(), expected, 1e-10, 1e-10);
        expect_close(Riemann(L,R).HLLC(), expected, 1e-10, 1e-10);
    #ifdef MHD
        expect_close(Riemann(L,R).HLLD(), expected, 1e-10, 1e-10);
    #endif
        expect_close(Riemann(L,R).HLLE(), expected, 1e-10, 1e-10);
    }
}
void DRAGON_Test::verify_hll_supersonic_upwind_transverse() {
    {
        PrimitiveState L = make_state(1.0, 10.0, 2.0, -3.0, 1.0);
        PrimitiveState R = make_state(0.125, 10.0, -8.0, 9.0, 0.1);

        ConservativeState expected = ConservativeState(L).flux(L.v);

        expect_close(Riemann(L,R).exact().flux(), expected, 1e-10, 1e-10);
        expect_close(Riemann(L,R).HLL(),          expected, 1e-10, 1e-10);
        expect_close(Riemann(L,R).HLLC(),         expected, 1e-10, 1e-10);
    #ifdef MHD
        expect_close(Riemann(L,R).HLLD(), expected, 1e-10, 1e-10);
    #endif
        expect_close(Riemann(L,R).HLLE(), expected, 1e-10, 1e-10);
    }

    {
        PrimitiveState L = make_state(1.0, -10.0, 2.0, -3.0, 1.0);
        PrimitiveState R = make_state(0.125, -10.0, -8.0, 9.0, 0.1);

        ConservativeState expected = ConservativeState(R).flux(R.v);

        expect_close(Riemann(L,R).exact().flux(), expected, 1e-10, 1e-10);
        expect_close(Riemann(L,R).HLL(),          expected, 1e-10, 1e-10);
        expect_close(Riemann(L,R).HLLC(),         expected, 1e-10, 1e-10);
    #ifdef MHD
        expect_close(Riemann(L,R).HLLD(), expected, 1e-10, 1e-10);
    #endif
        expect_close(Riemann(L,R).HLLE(), expected, 1e-10, 1e-10);
    }
}
void DRAGON_Test::verify_hll_symmetry() {
    PrimitiveState L = make_state(1.0, 0.35, 0.2, -0.1, 1.0);
    PrimitiveState R = make_state(0.7, -0.45, -0.3, 0.4, 0.8);

    ConservativeState F = Riemann(L, R).HLL();
    ConservativeState reflected = Riemann(mirrored_state(R), mirrored_state(L)).HLL();
    expect_close(reflected, mirrored_flux(F));
    
    F = Riemann(L, R).HLLC();
    reflected = Riemann(mirrored_state(R), mirrored_state(L)).HLLC();
    expect_close(reflected, mirrored_flux(F));
    
    F = Riemann(L, R).HLLE();
    reflected = Riemann(mirrored_state(R), mirrored_state(L)).HLLE();
    expect_close(reflected, mirrored_flux(F));
    
#ifdef MHD
    L = make_mhd_state(1.0, 0.35, 0.2, -0.1, 1.0, 0.4, -0.3, 0.2);
    R = make_mhd_state(0.7, -0.45, -0.3, 0.4, 0.8, 0.4, 0.1, -0.5);
    
    F = Riemann(L, R).HLL();
    reflected = Riemann(mirrored_state(R), mirrored_state(L)).HLL();
    expect_close(reflected, mirrored_flux(F));
    
    F = Riemann(L, R).HLLD();
    reflected = Riemann(mirrored_state(R), mirrored_state(L)).HLLD();
    expect_close(reflected, mirrored_flux(F));
    
    F = Riemann(L, R).HLLE();
    reflected = Riemann(mirrored_state(R), mirrored_state(L)).HLLE();
    expect_close(reflected, mirrored_flux(F));
#endif
}

#ifdef MHD


//MARK: MHD Test Problems (HLLD)

void DRAGON_Test::verify_hlld_equal_state_nonzero_b() {
    PrimitiveState W = make_mhd_state(1.0, 0.75, 0.2, -0.1, 1.0, 0.4, -0.3, 0.2);
    ConservativeState expected = ConservativeState(W).flux(W.v);

    Riemann problem(W, W);
    ConservativeState F = problem.HLLD();

    expect_close(F, expected, 1e-10, 1e-10);
    expect_close(problem.L, W);
    expect_close(problem.R, W);
    expect_finite(F);
}

void DRAGON_Test::verify_hlld_supersonic_exterior_regions() {
    {
        PrimitiveState L = make_mhd_state(1.0,   10.0, 0.5, -0.3, 1.0, 0.4, -0.3, 0.2);
        PrimitiveState R = make_mhd_state(0.125, 10.0, -0.2, 0.4, 0.1, 0.4, 0.1, -0.5);
        ConservativeState expected = ConservativeState(L).flux(L.v);

        expect_close(Riemann(L,R).HLLD(), expected, 1e-10, 1e-10);
    }
    {
        PrimitiveState L = make_mhd_state(1.0,   -10.0, 0.5, -0.3, 1.0, 0.4, -0.3, 0.2);
        PrimitiveState R = make_mhd_state(0.125, -10.0, -0.2, 0.4, 0.1, 0.4, 0.1, -0.5);
        ConservativeState expected = ConservativeState(R).flux(R.v);

        expect_close(Riemann(L,R).HLLD(), expected, 1e-10, 1e-10);
    }
}

void DRAGON_Test::verify_hlld_averages_normal_field() {
    PrimitiveState L = make_mhd_state(1.0, 0.2, 0.5, -0.3, 1.0, 0.8, -0.3, 0.2);
    PrimitiveState R = make_mhd_state(0.7, -0.4, -0.2, 0.4, 0.8, -0.2, 0.1, -0.5);
    double expectedBx = 0.5 * (L.B.x + R.B.x);

    Riemann problem(L, R);
    ConservativeState F = problem.HLLD();

    assert(approx(problem.L.B.x, expectedBx));
    assert(approx(problem.R.B.x, expectedBx));
    assert(approx(F.B.x, 0.0));
    expect_finite(F);
}

void DRAGON_Test::verify_hlld_explicit_normal_field() {
    PrimitiveState L = make_mhd_state(1.0, 0.2, 0.5, -0.3, 1.0, 0.25, -0.3, 0.2);
    PrimitiveState R = make_mhd_state(0.7, -0.4, -0.2, 0.4, 0.8, 0.25, 0.1, -0.5);

    Riemann problem(L, R);
    ConservativeState F = problem.HLLD();

    assert(approx(problem.L.B.x, 0.25));
    assert(approx(problem.R.B.x, 0.25));
    assert(approx(F.B.x, 0.0));
    expect_finite(F);
}

void DRAGON_Test::verify_hlld_star_regions_finite() {
    PrimitiveState cases[][2] = {
        { make_mhd_state(1.0,  -0.7, 0.6, -0.2, 1.0,  0.4,  0.3, -0.2),
          make_mhd_state(0.8,   0.5, -0.4, 0.3, 0.7,  0.4, -0.2,  0.5) },
        { make_mhd_state(1.0,  -0.7, 0.6, -0.2, 1.0, -0.4,  0.3, -0.2),
          make_mhd_state(0.8,   0.5, -0.4, 0.3, 0.7, -0.4, -0.2,  0.5) },
        { make_mhd_state(1.0,  -1.2, 0.1, 0.5, 0.6,  0.0,  0.7, -0.4),
          make_mhd_state(1.2,   1.0, -0.3, 0.2, 0.9,  0.0, -0.5,  0.3) },
        { make_mhd_state(5.99924, 19.5975, 0.5, -0.2, 460.894,  1.0,  0.2, -0.6),
          make_mhd_state(5.99242,-6.19633, -0.4, 0.3, 46.0950,  1.0, -0.3,  0.4) }
    };

    for(auto& pair : cases) {
        ConservativeState F = Riemann(pair[0], pair[1]).HLLD();

        expect_finite(F);
        assert((ConservativeState(pair[0]) - F*0.001).isPhysical());
        assert((ConservativeState(pair[1]) + F*0.001).isPhysical());
    }
}



#endif
