//
//  Testing.hpp
//  DRAGON/Testing
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
void verify_fluid_element(bool output = true);
PrimitiveState make_state(double rho, double vx, double vy, double vz, double p);
//Data Structure & Coversion
void verify_constructors();
void verify_conversion();
void verify_flux();
void verify_flux_add();
void verify_enthalpy();
void verify_swaps_P();//Primitive
void verify_swaps_C();//Conservative
//Arithmetic
void verify_add();
void verify_sub();
void verify_mult();
void verify_div();



//MARK: Riemann Solvers
void verify_riemann(bool output = true);
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
void verify_exact_supersonic_upwind_transverse();
// HLL Tests
void verify_hll_equal_state();
void verify_hll_stationary_contact();
void verify_hll_supersonic_upwind();
void verify_hll_manual_wave_speeds();
void verify_hllc_manual_wave_speeds();
void verify_hll_supersonic_upwind_transverse();
//Roe Tests
void verify_roe_equal_state();
void verify_roe_stationary_contact();
void verify_roe_supersonic_upwind();
void verify_roe_supersonic_upwind_transverse();
//Finiteness
void verify_approximate_solver_finiteness();


//MARK: Grid Structure
void verify_grid(bool output = true);
void verify_grid1D();
void verify_grid2D();
void verify_grid3D();
PrimitiveState make_tagged_state(double tag);



//MARK: Boundary Conditions
void verify_boundary(bool output = true);
//Constructors
void verify_boundary_constructors();
//Outflow
void verify_boundary_outflow();
void verify_boundary_outflow_1D();
void verify_boundary_outflow_2D();
void verify_boundary_outflow_3D();
void verify_boundary_outflow_1D_gated();
void verify_boundary_outflow_2D_gated();
void verify_boundary_outflow_3D_gated_X();
void verify_boundary_outflow_3D_gated_Y();
void verify_boundary_outflow_3D_gated_Z();
//Reflective
void verify_boundary_reflective();
void verify_boundary_reflective_1D();
void verify_boundary_reflective_2D();
void verify_boundary_reflective_3D();
//Periodic
void verify_boundary_periodic();
void verify_boundary_periodic_1D();
void verify_boundary_periodic_2D();
void verify_boundary_periodic_3D();
//Fixed
void verify_boundary_fixed();
void verify_boundary_fixed_1D();
void verify_boundary_fixed_2D();
void verify_boundary_fixed_3D();
//Composition
void verify_boundary_composition();
void verify_boundary_set_missing_faces_outflow_2D();
void verify_boundary_composition_3D();
void verify_boundary_composition_order();


//MARK: CFL
void verify_cfl(bool output=true);
//Max Speed
void verify_cfl_max_speed_1D();
void verify_cfl_max_speed_2D();
void verify_cfl_max_speed_3D();
//Add Speeds
void verify_cfl_add_speed_1D();
void verify_cfl_add_speed_2D();
void verify_cfl_add_speed_3D();
//L2 norm
void verify_cfl_pow_speed_1D();
void verify_cfl_pow_speed_2D();
void verify_cfl_pow_speed_3D();
//Dispatch
void verify_cfl_dispatch_add_3D();
void verify_cfl_dispatch_max_3D();
void verify_cfl_dispatch_pow_3D();
//1D
void verify_cfl_time_1D_uniform();
void verify_cfl_time_1D_uses_fastest_cell();
void verify_cfl_time_1D_ignores_ghost_cells();
//2D
void verify_cfl_time_2D_uniform();
void verify_cfl_time_2D_visits_last_cell();
void verify_cfl_time_2D_ignores_ghost_cells();
//3D
void verify_cfl_time_3D_uniform();
void verify_cfl_time_3D_uses_fastest_cell();
void verify_cfl_time_3D_ignores_ghost_cells();

//MARK: Godunov Scheme
//1D
void verify_godunov_1D(bool output = true);
void verify_god_uniform_stationary_1D();
void verify_god_uniform_moving_1D();
void verify_god_periodic_conservation_1D();
void verify_god_dt0_1D();
//2D
void verify_godunov_2D_Split(bool output = true);
void verify_godunov_2D_Unsplit(bool output = true);
void verify_god_uniform_stationary_2D(bool split = true);
void verify_god_uniform_moving_2D(bool split = true);
void verify_god_periodic_conservation_2D(bool split = true);
void verify_god_dt0_2D(bool split = true);
void verify_2D_X_match_1D(bool split = true);
void verify_2D_Y_match_1D(bool split = true);
//3D Split
void verify_godunov_3D_Split(bool output = true);
void verify_godunov_3D_Unsplit(bool output = true);
void verify_god_uniform_stationary_3D(bool split = true);
void verify_god_uniform_moving_3D(bool split = true);
void verify_god_periodic_conservation_3D(bool split = true);
void verify_god_dt0_3D(bool split = true);
void verify_3D_X_match_1D(bool split = true);
void verify_3D_Y_match_1D(bool split = true);
void verify_3D_Z_match_1D(bool split = true);

}
#endif /* Testing_hpp */
