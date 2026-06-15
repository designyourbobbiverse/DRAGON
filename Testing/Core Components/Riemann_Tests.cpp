//
//  Riemann_Tests.cpp
//  DRAGON/Testing
//
//  Created by Bobbie Markwick on 12/05/2026.
//

#include "Testing.hpp"
#include "Riemann.hpp"
#include <iostream>

using namespace DRAGON_Test;

void DRAGON_Test::verify_riemann(bool output){
    if(output) std::cout << "Exact Riemann Solver: \n";
    // Riemann construction / solution setup
    if(output) std::cout << "- Construction: ";
    verify_riemann_constructor();
    verify_riemann_solution_constructor();
    if(output) std::cout << "Passed\n";
    // Exact solver Tests
    if(output) std::cout << "- Equation: ";
    verify_riemann_f();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "- Sampler: ";
    verify_sample_mirror_restores_state();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "- Selected Problems: ";
    verify_exact_equal_state();
    verify_exact_stationary_contact();
    verify_exact_supersonic_upwind();
    verify_exact_supersonic_upwind_transverse();
    verify_exact_sod();
    if(output) std::cout << "All Tests Passed\n";

    // HLL Tests
    if(output) std::cout << "HLL/C Riemann Solver: ";
    verify_hll_equal_state();
    verify_hll_stationary_contact();
    verify_hll_supersonic_upwind();
    verify_hll_supersonic_upwind_transverse();
    verify_hll_manual_wave_speeds();
    verify_hllc_manual_wave_speeds();
    if(output) std::cout << "All Tests Passed\n";


    //Roe Tests
    if(output) std::cout << "Roe Riemann Solver: ";
    verify_roe_equal_state();
    verify_roe_stationary_contact();
    verify_roe_supersonic_upwind();
    verify_roe_supersonic_upwind_transverse();
    if(output) std::cout << "All Tests Passed\n";

    //Finiteness
    //if(output) std::cout << "Solver Finiteness: ";
    verify_approximate_solver_finiteness();
    //if(output) std::cout << "Passed\n";
    if(output) std::cout << "All Riemann Solver Tests Passed.\n\n";
}




//MARK: Constructors
void DRAGON_Test::verify_riemann_constructor(){
    PrimitiveState L = make_state(1.0, 2.0, 3.0, 4.0, 5.0);
    PrimitiveState R = make_state(6.0, 7.0, 8.0, 9.0, 10.0);
    Riemann problem(L, R);
    expect_close(problem.L, L);
    expect_close(problem.R, R);
}
void DRAGON_Test::verify_riemann_solution_constructor(){
    PrimitiveState L = make_state(1.0, 2.0, 3.0, 4.0, 5.0);
    PrimitiveState R = make_state(6.0, 7.0, 8.0, 9.0, 10.0);
    RiemannSolution S(Riemann(L, R));
    //Initial States Copied
    expect_close(S.wL, L);
    expect_close(S.wR, R);
    //Transverse Velocities
    assert(approx(S.sL.v.y, L.v.y));
    assert(approx(S.sL.v.z, L.v.z));
    assert(approx(S.sR.v.y, R.v.y));
    assert(approx(S.sR.v.z, R.v.z));
}

//MARK: Exact Solver Tests
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

void DRAGON_Test::verify_exact_equal_state() {
    PrimitiveState W = make_state(1.0, 0.75, 0.2, -0.1, 1.0);
    ConservativeState expected = ConservativeState(W).flux(W.v.x);

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
        ConservativeState expected = ConservativeState(L).flux(L.v.x);
        expect_close(Riemann(L,R).exact().flux(), expected, 1e-10, 1e-10);
    }
    {
        PrimitiveState L = make_state(1.0,   -10.0, 0.0, 0.0, 1.0);
        PrimitiveState R = make_state(0.125, -10.0, 0.0, 0.0, 0.1);
        ConservativeState expected = ConservativeState(R).flux(R.v.x);
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

void DRAGON_Test::verify_exact_supersonic_upwind_transverse() {
    {
        PrimitiveState L = make_state(1.0, 10.0, 2.0, -3.0, 1.0);
        PrimitiveState R = make_state(0.125, 10.0, -8.0, 9.0, 0.1);
        ConservativeState expected = ConservativeState(L).flux(L.v.x);
        expect_close(Riemann(L,R).exact().flux(), expected, 1e-10, 1e-10);
    }

    {
        PrimitiveState L = make_state(1.0, -10.0, 2.0, -3.0, 1.0);
        PrimitiveState R = make_state(0.125, -10.0, -8.0, 9.0, 0.1);
        ConservativeState expected = ConservativeState(R).flux(R.v.x);
        expect_close(Riemann(L,R).exact().flux(), expected, 1e-10, 1e-10);
    }
}
void DRAGON_Test::verify_hll_supersonic_upwind_transverse() {
    {
        PrimitiveState L = make_state(1.0, 10.0, 2.0, -3.0, 1.0);
        PrimitiveState R = make_state(0.125, 10.0, -8.0, 9.0, 0.1);

        ConservativeState expected = ConservativeState(L).flux(L.v.x);

        expect_close(Riemann(L,R).exact().flux(), expected, 1e-10, 1e-10);
        expect_close(Riemann(L,R).HLL(),          expected, 1e-10, 1e-10);
        expect_close(Riemann(L,R).HLLC(),         expected, 1e-10, 1e-10);
        expect_close(Riemann(L,R).Roe(),          expected, 1e-10, 1e-10);
    }

    {
        PrimitiveState L = make_state(1.0, -10.0, 2.0, -3.0, 1.0);
        PrimitiveState R = make_state(0.125, -10.0, -8.0, 9.0, 0.1);

        ConservativeState expected = ConservativeState(R).flux(R.v.x);

        expect_close(Riemann(L,R).exact().flux(), expected, 1e-10, 1e-10);
        expect_close(Riemann(L,R).HLL(),          expected, 1e-10, 1e-10);
        expect_close(Riemann(L,R).HLLC(),         expected, 1e-10, 1e-10);
        expect_close(Riemann(L,R).Roe(),          expected, 1e-10, 1e-10);
    }
}

//MARK: HLL/HLLC Tests
void DRAGON_Test::verify_hll_equal_state() {
    PrimitiveState W = make_state(1.0, 0.75, 0.2, -0.1, 1.0);
    ConservativeState expected = ConservativeState(W).flux(W.v.x);

    expect_close(Riemann(W,W).HLL(), expected, 1e-10, 1e-10);
    expect_close(Riemann(W,W).HLLC(), expected, 1e-10, 1e-10);
}
void DRAGON_Test::verify_hll_stationary_contact() {
    PrimitiveState L = make_state(1.0, 0.0, 0.0, 0.0, 1.0);
    PrimitiveState R = make_state(2.0, 0.0, 0.0, 0.0, 1.0);
    ConservativeState expected;
    expected.rho = 0.0;
    expected.p.x  = 1.0;
    expected.p.y  = 0.0;
    expected.p.z  = 0.0;
    expected.E   = 0.0;

    expect_close(Riemann(L,R).HLLC(),expected, 1e-10, 1e-10);
    // HLL may diffuse contacts in evolution, but this exact stationary case should still be finite.
    expect_finite(Riemann(L,R).HLL());
}

void DRAGON_Test::verify_hll_supersonic_upwind() {
    {
        PrimitiveState L = make_state(1.0,   10.0, 0.0, 0.0, 1.0);
        PrimitiveState R = make_state(0.125, 10.0, 0.0, 0.0, 0.1);
        ConservativeState expected = ConservativeState(L).flux(L.v.x);
        expect_close(Riemann(L,R).HLL(), expected, 1e-10, 1e-10);
        expect_close(Riemann(L,R).HLLC(), expected, 1e-10, 1e-10);
    }
    {
        PrimitiveState L = make_state(1.0,   -10.0, 0.0, 0.0, 1.0);
        PrimitiveState R = make_state(0.125, -10.0, 0.0, 0.0, 0.1);
        ConservativeState expected = ConservativeState(R).flux(R.v.x);
        expect_close(Riemann(L,R).HLL(), expected, 1e-10, 1e-10);
        expect_close(Riemann(L,R).HLLC(), expected, 1e-10, 1e-10);
    }
}
void DRAGON_Test::verify_hll_manual_wave_speeds() {
    PrimitiveState L = make_state(1.0,   0.0, 0.0, 0.0, 1.0);
    PrimitiveState R = make_state(0.125, 0.0, 0.0, 0.0, 0.1);

    Riemann problem(L, R);

    ConservativeState UL(L), UR(R);
    ConservativeState FL = UL.flux(L.v.x);
    ConservativeState FR = UR.flux(R.v.x);

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

    ConservativeState FL = ConservativeState(L).flux(L.v.x);
    ConservativeState FR = ConservativeState(R).flux(R.v.x);

    expect_close(problem.HLLC( 0.1, 2.0), FL);
    expect_close(problem.HLLC(-2.0,-0.1), FR);

    expect_finite(problem.HLLC(-1.0, 2.0));
}



//MARK: Roe Tests
void DRAGON_Test::verify_roe_equal_state() {
    PrimitiveState W = make_state(1.0, 0.75, 0.2, -0.1, 1.0);

    Riemann problem(W, W);
    ConservativeState expected = ConservativeState(W).flux(W.v.x);

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
        ConservativeState expected = ConservativeState(L).flux(L.v.x);
        expect_close(Riemann(L,R).Roe(), expected, 1e-10, 1e-10);
    }
    {
        PrimitiveState L = make_state(1.0,   -10.0, 0.0, 0.0, 1.0);
        PrimitiveState R = make_state(0.125, -10.0, 0.0, 0.0, 0.1);
        ConservativeState expected = ConservativeState(R).flux(R.v.x);
        expect_close(Riemann(L,R).Roe(), expected, 1e-10, 1e-10);
    }
}
void DRAGON_Test::verify_roe_supersonic_upwind_transverse() {
    {
        PrimitiveState L = make_state(1.0, 10.0, 2.0, -3.0, 1.0);
        PrimitiveState R = make_state(0.125, 10.0, -8.0, 9.0, 0.1);
        ConservativeState expected = ConservativeState(L).flux(L.v.x);
        expect_close(Riemann(L,R).Roe(),          expected, 1e-10, 1e-10);
    }
    {
        PrimitiveState L = make_state(1.0, -10.0, 2.0, -3.0, 1.0);
        PrimitiveState R = make_state(0.125, -10.0, -8.0, 9.0, 0.1);
        ConservativeState expected = ConservativeState(R).flux(R.v.x);
        expect_close(Riemann(L,R).Roe(),          expected, 1e-10, 1e-10);
    }
}

//MARK: Verify Finiteness
void DRAGON_Test::expect_finite(const ConservativeState& U) {
    assert(std::isfinite(U.rho));
    assert(std::isfinite(U.p.x));
    assert(std::isfinite(U.p.y));
    assert(std::isfinite(U.p.z));
    assert(std::isfinite(U.E));
}
void DRAGON_Test::expect_finite(const PrimitiveState& W) {
    assert(std::isfinite(W.rho));
    assert(std::isfinite(W.v.x));
    assert(std::isfinite(W.v.y));
    assert(std::isfinite(W.v.z));
    assert(std::isfinite(W.p));
}

void DRAGON_Test::verify_approximate_solver_finiteness() {
    PrimitiveState cases[][2] = {
        { make_state(1.0,   0.0, 0.0, 0.0, 1.0),  make_state(0.125, 0.0, 0.0, 0.0, 0.1)  }, // Sod-like
        { make_state(0.125, 0.0, 0.0, 0.0, 0.1),  make_state(1.0,   0.0, 0.0, 0.0, 1.0)  }, // reverse Sod
        { make_state(1.0,  -2.0, 0.0, 0.0, 0.4),  make_state(1.0,   2.0, 0.0, 0.0, 0.4)  }, // rarefaction-ish
        { make_state(5.99924, 19.5975, 0.0, 0.0, 460.894),
          make_state(5.99242,-6.19633, 0.0, 0.0, 46.0950) }, // strong shock-ish
        { make_state(1.0, 0.0, 0.3, -0.2, 1.0),
          make_state(0.5, 0.0, -0.4, 0.1, 0.5) } // transverse velocities
    };

    for (auto& pair : cases) {
        Riemann problem(pair[0], pair[1]);

        expect_finite(problem.exact().flux());
        expect_finite(problem.HLL());
        expect_finite(problem.HLLC());
        expect_finite(problem.Roe());
    }
}
