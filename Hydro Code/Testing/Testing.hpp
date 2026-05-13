//
//  Testing.hpp
//  Hydro Code
//
//  Created by Bobbie Markwick on 12/05/2026.
//

#ifndef Testing_hpp
#define Testing_hpp

#include <stdio.h>
#include <cassert>
#include "FluidElement.hpp"

namespace DRAGON_Test{


//MARK: Numerical Close-ness
bool approx(double a, double b, double rel = 1e-12, double abs = 1e-12);
void expect_close(const ConservativeState& a, const ConservativeState& b,  double rel = 1e-12, double abs = 1e-12);
void expect_close(const PrimitiveState& a, const PrimitiveState& b, double rel = 1e-12, double abs = 1e-12);
void expect_finite(const ConservativeState& U);
void expect_finite(const PrimitiveState& W);

//MARK: Fluid Element
void verify_fluid_element();
//Data Structure & Coversion
void verify_constructors();
void verify_conversion();
void verify_flux();
void verify_enthalpy();
void verify_swaps_P();//Primitive
void verify_swaps_C();//Conservative
//Arithmetic
void verify_add();
void verify_sub();
void verify_mult();
void verify_div();



//MARK: Riemann Solvers
void verify_riemann();
// Riemann construction / solution setup
void verify_riemann_constructor();
void verify_riemann_solution_constructor();
// Exact solver Tests
void verify_exact_equal_state();
void verify_exact_stationary_contact();
void verify_exact_supersonic_upwind();
void verify_exact_sod();
void verify_riemann_f();
void verify_sample_mirror_restores_state();
// HLL Tests
void verify_hll_equal_state();
void verify_hll_stationary_contact();
void verify_hll_supersonic_upwind();
void verify_hll_manual_wave_speeds();
void verify_hllc_manual_wave_speeds();
//Roe Tests
void verify_roe_equal_state();
void verify_roe_stationary_contact();
void verify_roe_supersonic_upwind();
//Finiteness
void verify_approximate_solver_finiteness();
}


#endif /* Testing_hpp */
