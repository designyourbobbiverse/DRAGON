//
//  Riemann_Tests.cpp
//  DRAGON/Testing/Core/Riemann
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
    verify_approximate_solver_finiteness();
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
    verify_riemann_exact_symmetry();
    if(output) std::cout << "All Tests Passed\n";
    
    // HLL Tests
    if(output) std::cout << "HLL/C/E Riemann Solvers: ";
    verify_hll_equal_state();
    verify_hll_stationary_contact();
    verify_hll_supersonic_upwind();
    verify_hll_supersonic_upwind_transverse();
    verify_hll_symmetry();
    verify_hll_manual_wave_speeds();
    verify_hllc_manual_wave_speeds();
    if(output) std::cout << "All Tests Passed\n";
#ifdef MHD
    if(output) std::cout << "HLLD Riemann Solver: ";
    verify_hlld_equal_state_nonzero_b();
    verify_hlld_supersonic_exterior_regions();
    verify_hlld_averages_normal_field();
    verify_hlld_explicit_normal_field();
    verify_hlld_star_regions_finite();
    if(output) std::cout << "All Tests Passed\n";
#endif


    //Roe Tests
    if(output) std::cout << "Roe Riemann Solver: ";
    verify_roe_equal_state();
    verify_roe_stationary_contact();
    verify_roe_supersonic_upwind();
    verify_roe_supersonic_upwind_transverse();
    verify_roe_symmetry();
    verify_roe_entropy_fix_rarefactions();
    if(output) std::cout << "All Tests Passed\n";

    
    if(output) std::cout << "Riemann Dispatch: ";
    verify_riemann_flux_dispatch();
    verify_riemann_flux_safety_check();
    verify_riemann_flux_dimension_wrappers();
#ifdef MHD
    verify_riemann_flux_dimension_wrappers_mhd();
#endif
    if(output) std::cout << "Passed\n";


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

//MARK: Flux Dispatch
void DRAGON_Test::verify_riemann_flux_dispatch(){
    PrimitiveState L = make_state(1.0, 0.25, 0.2, -0.1, 1.0);
    PrimitiveState R = make_state(0.8, -0.15, -0.3, 0.4, 0.7);

    CONFIG::riemann_choice = RIEMANN_EXACT;
    expect_close(Riemann(L,R).flux(), Riemann(L,R).exact().flux(), 1e-10, 1e-10);

    CONFIG::riemann_choice = RIEMANN_HLL;
    expect_close(Riemann(L,R).flux(), Riemann(L,R).HLL(), 1e-10, 1e-10);

    CONFIG::riemann_choice = RIEMANN_HLLE;
    expect_close(Riemann(L,R).flux(), Riemann(L,R).HLLE(), 1e-10, 1e-10);

    CONFIG::riemann_choice = RIEMANN_HLLC;
    expect_close(Riemann(L,R).flux(), Riemann(L,R).HLLC(), 1e-10, 1e-10);

    CONFIG::riemann_choice = RIEMANN_ROE;
    expect_close(Riemann(L,R).flux(), Riemann(L,R).Roe(), 1e-10, 1e-10);

#ifdef MHD
    CONFIG::riemann_choice = RIEMANN_HLLD;
    expect_close(Riemann(L,R).flux(), Riemann(L,R).HLLD(), 1e-10, 1e-10);

    CONFIG::riemann_choice = RIEMANN_HLLX;
    expect_close(Riemann(L,R).flux(), Riemann(L,R).HLLD(), 1e-10, 1e-10);
#else
    CONFIG::riemann_choice = RIEMANN_HLLX;
    expect_close(Riemann(L,R).flux(), Riemann(L,R).HLLC(), 1e-10, 1e-10);
#endif

    CONFIG::riemann_choice = RIEMANN_EXACT;
}

void DRAGON_Test::verify_riemann_flux_safety_check(){
    PrimitiveState W = make_state(1.0, 0.2, -0.3, 0.4, 1.0);

    int previous = CONFIG::riemann_choice;
    CONFIG::riemann_choice = RIEMANN_HLL;
    ConservativeState unchecked = Riemann(W,W).flux(0.0);
    ConservativeState checked = Riemann(W,W).flux(0.01);

    expect_close(checked, unchecked, 1e-10, 1e-10);
    assert((ConservativeState(W) - checked*0.01).isPhysical());
    assert((ConservativeState(W) + checked*0.01).isPhysical());

#ifdef RIEMANN_VERIFY_FALLBACK
    PrimitiveState L = make_state(1.0, -3.0, 0.2, -0.1, 0.4);
    PrimitiveState R = make_state(1.0,  3.0, -0.3, 0.4, 0.4);
    CONFIG::riemann_choice = RIEMANN_HLLC;

    ConservativeState selected = Riemann(L, R).flux(0.0);
    ConservativeState fallback = Riemann(L, R).flux(1.0);

    bool selectedPhysical = (ConservativeState(L) - selected).isPhysical()
                         && (ConservativeState(R) + selected).isPhysical();
    if(!selectedPhysical) {
        ConservativeState hlle = Riemann(L, R).HLLE();
        ConservativeState exact = Riemann(L, R).exact().flux();
        assert(fallback == hlle || fallback == exact);
    }
    expect_finite(fallback);
#endif

    CONFIG::riemann_choice = previous;
}
//MARK: Dimension Wrappers
void DRAGON_Test::verify_riemann_flux_dimension_wrappers(){
    PrimitiveState L = make_state(1.0, 0.6, -0.2, 0.3, 1.0);
    PrimitiveState R = make_state(0.7, -0.4, 0.5, -0.1, 0.8);
    
    int previous =  CONFIG::riemann_choice;
    CONFIG::riemann_choice = RIEMANN_HLL;

    Riemann xProblem(L, R);
    expect_close(xProblem.flux_X(), Riemann(L,R).flux(), 1e-10, 1e-10);
    expect_close(xProblem.L, L);
    expect_close(xProblem.R, R);

    Riemann yProblem(L, R);
    ConservativeState expectedY = Riemann(L.swappedXY(), R.swappedXY()).flux().swappedXY();
    expect_close(yProblem.flux_Y(), expectedY, 1e-10, 1e-10);
    expect_close(yProblem.L, L);
    expect_close(yProblem.R, R);

    Riemann zProblem(L, R);
    ConservativeState expectedZ = Riemann(L.swappedXZ(), R.swappedXZ()).flux().swappedXZ();
    expect_close(zProblem.flux_Z(), expectedZ, 1e-10, 1e-10);
    expect_close(zProblem.L, L);
    expect_close(zProblem.R, R);

    CONFIG::riemann_choice = previous;
}

#ifdef MHD
void DRAGON_Test::verify_riemann_flux_dimension_wrappers_mhd(){
    PrimitiveState L = make_mhd_state(1.0, 0.6, -0.2, 0.3, 1.0, 0.4, -0.3, 0.2);
    PrimitiveState R = make_mhd_state(0.7, -0.4, 0.5, -0.1, 0.8, 0.4, -0.3, 0.2);

    int previous = CONFIG::riemann_choice;
    int choices[] = { RIEMANN_HLL, RIEMANN_HLLE, RIEMANN_HLLD, RIEMANN_HLLX };

    for(int choice : choices) {
        CONFIG::riemann_choice = choice;

        Riemann xProblem(L, R);
        expect_close(xProblem.flux_X(), Riemann(L,R).flux(), 1e-10, 1e-10);
        expect_close(xProblem.L, L);
        expect_close(xProblem.R, R);

        Riemann yProblem(L, R);
        ConservativeState expectedY = Riemann(L.swappedXY(), R.swappedXY()).flux().swappedXY();
        expect_close(yProblem.flux_Y(), expectedY, 1e-10, 1e-10);
        expect_close(yProblem.L, L);
        expect_close(yProblem.R, R);

        Riemann zProblem(L, R);
        ConservativeState expectedZ = Riemann(L.swappedXZ(), R.swappedXZ()).flux().swappedXZ();
        expect_close(zProblem.flux_Z(), expectedZ, 1e-10, 1e-10);
        expect_close(zProblem.L, L);
        expect_close(zProblem.R, R);
    }

    CONFIG::riemann_choice = previous;
}
#endif


//MARK: Verify Finiteness
void DRAGON_Test::expect_finite(const ConservativeState& U) {
    assert(std::isfinite(U.rho));
    assert(std::isfinite(U.p.x));
    assert(std::isfinite(U.p.y));
    assert(std::isfinite(U.p.z));
    assert(std::isfinite(U.E));
#ifdef MHD
    assert(std::isfinite(U.B.x));
    assert(std::isfinite(U.B.y));
    assert(std::isfinite(U.B.z));
#endif
}
void DRAGON_Test::expect_finite(const PrimitiveState& W) {
    assert(std::isfinite(W.rho));
    assert(std::isfinite(W.v.x));
    assert(std::isfinite(W.v.y));
    assert(std::isfinite(W.v.z));
    assert(std::isfinite(W.p));
#ifdef MHD
    assert(std::isfinite(W.B.x));
    assert(std::isfinite(W.B.y));
    assert(std::isfinite(W.B.z));
#endif
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
        expect_finite(problem.HLLE());
        expect_finite(problem.HLLC());
        expect_finite(problem.Roe());
    }

#ifdef MHD
    PrimitiveState mhdCases[][2] = {
        { make_mhd_state(1.0,   0.0, 0.2, -0.1, 1.0,  0.4,  0.3, -0.2),
          make_mhd_state(0.125, 0.0, -0.3, 0.4, 0.1, -0.1, -0.2,  0.5) },
        { make_mhd_state(1.0,  -2.0, 0.4, 0.1, 0.4,  0.2,  0.5, -0.3),
          make_mhd_state(1.0,   2.0, -0.2, 0.3, 0.4,  0.2, -0.4,  0.1) },
        { make_mhd_state(5.99924, 19.5975, 0.5, -0.2, 460.894,  1.0,  0.2, -0.6),
          make_mhd_state(5.99242,-6.19633, -0.4, 0.3, 46.0950, -0.5, -0.3,  0.4) }
    };

    for (auto& pair : mhdCases) {
        Riemann problem(pair[0], pair[1]);

        expect_finite(problem.HLL());
        expect_finite(problem.HLLE());
        expect_finite(problem.HLLD());
    }
#endif
}
