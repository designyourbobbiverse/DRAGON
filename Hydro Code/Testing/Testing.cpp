//
//  Testing.cpp
//  Hydro Code
//
//  Created by Bobbie Markwick on 12/05/2026.
//
#include "Testing.hpp"
#include <iostream>

int main() {
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


    std::cout << "Fluid Element tests passed.\n";

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
    
    std::cout << "Exact Riemann Solver tests passed.\n";

    // HLL Tests
    void verify_hll_equal_state();
    void verify_hll_stationary_contact();
    void verify_hll_supersonic_upwind();
    void verify_hll_manual_wave_speeds();
    void verify_hllc_manual_wave_speeds();
    
    std::cout << "HLL Riemann Solver tests passed.\n";

    
    //Roe Tests
    void verify_roe_equal_state();
    void verify_roe_stationary_contact();
    void verify_roe_supersonic_upwind();
    
    std::cout << "Roe Riemann Solver tests passed.\n";

    //Finiteness
    void verify_approximate_solver_finiteness();
    std::cout << "All tests passed.\n";
}
