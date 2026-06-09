//
//  Testing.cpp
//  DRAGON/Testing
//
//  Created by Bobbie Markwick on 12/05/2026.
//
#include "Testing.hpp"
#include <iostream>

using namespace DRAGON_Test;

int main() {
    //MARK: Fluid Element
    verify_fluid_element();
    std::cout << "Fluid Arithmetic tests passed.\n\n";

    //MARK: Riemann Solvers
    // Riemann construction / solution setup
    verify_riemann_constructor();
    verify_riemann_solution_constructor();
    // Exact solver Tests
    verify_riemann_f();
    verify_sample_mirror_restores_state();
    verify_exact_equal_state();
    verify_exact_stationary_contact();
    verify_exact_supersonic_upwind();
    verify_exact_sod();

    std::cout << "Exact Riemann Solver tests passed.\n";

    // HLL Tests
    verify_hll_equal_state();
    verify_hll_stationary_contact();
    verify_hll_supersonic_upwind();
    verify_hll_manual_wave_speeds();
    verify_hllc_manual_wave_speeds();

    std::cout << "HLL/C Riemann Solver tests passed.\n";


    //Roe Tests
    verify_roe_equal_state();
    verify_roe_stationary_contact();
    verify_roe_supersonic_upwind();

    std::cout << "Roe Riemann Solver tests passed.\n";

    //Finiteness
    verify_approximate_solver_finiteness();
    std::cout << "All Riemann Solver tests passed.\n\n";

    //MARK: Grid Structure
    verify_grid();
    std::cout << "Grid Structure tests passed.\n\n";

    
    //MARK: Boundary Conditions
    verify_boundary_constructors();
    //Outflow
    verify_boundary_outflow();
    std::cout << "Outflow Boundary tests passed.\n";
    //Reflective
    verify_boundary_reflective();
    std::cout << "Reflective Boundary tests passed.\n";
    //Periodic
    verify_boundary_periodic();
    std::cout << "Periodic Boundary tests passed.\n";
    //Fixed
    verify_boundary_fixed();
    std::cout << "Fixed State Boundary tests passed.\n";
    //Composition
    verify_boundary_composition();
    std::cout << "Boundary Composition tests passed.\n";

    std::cout << "All Boundary Setup tests passed.\n\n";

    
    
    std::cout << "All tests passed.\n";

    
    
}
