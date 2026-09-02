//
//  Testing.hpp
//  DRAGON/Testing
//
//  Created by Bobbie Markwick on 12/05/2026.
//

#ifndef Testing_hpp
#define Testing_hpp

#include <cassert>

#include "Config.h"
#include "FluidElement/FluidElement.hpp"
#include "Hydro/ExtendedArray/ExtendedArray.hpp"
using namespace DRAGON;

namespace DRAGON_Test{



//MARK: Numerical Close-ness
bool approx(double a, double b, double rel = 1e-15, double abs = 1e-18);
void expect_close(const vec3& a, const vec3& b,  double rel = 1e-15, double abs = 1e-18);
void expect_close(const ConservativeState& a, const ConservativeState& b,  double rel = 1e-15, double abs = 1e-18);
void expect_close(const PrimitiveState& a, const PrimitiveState& b, double rel = 1e-15, double abs = 1e-18);
void expect_finite(const ConservativeState& U);
void expect_finite(const PrimitiveState& W);


//MARK: Fluid Element
void verify_fluid_element(bool output = true);
PrimitiveState make_state(double rho, double vx, double vy, double vz, double p);
#ifdef MHD
PrimitiveState make_mhd_state(double rho, double vx, double vy, double vz, double p, double Bx, double By, double Bz);
#endif
//Data Structure & Coversion
void verify_constructors();
void verify_conversion();
void verify_flux();
void verify_flux_add();
void verify_enthalpy();
void verify_wavespeeds();
void verify_physicality();
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
void verify_riemann_exact_symmetry();
// Flux dispatch and wrappers
void verify_riemann_flux_dispatch();
void verify_riemann_flux_safety_check();
void verify_riemann_flux_dimension_wrappers();
void verify_riemann_flux_dimension_wrappers_mhd();
// HLL Tests
void verify_hll_equal_state();
void verify_hll_stationary_contact();
void verify_hll_supersonic_upwind();
void verify_hll_symmetry();
void verify_hll_manual_wave_speeds();
void verify_hllc_manual_wave_speeds();
void verify_hll_supersonic_upwind_transverse();
void verify_hlld_equal_state_nonzero_b();
void verify_hlld_supersonic_exterior_regions();
void verify_hlld_averages_normal_field();
void verify_hlld_explicit_normal_field();
void verify_hlld_star_regions_finite();
//Roe Tests
void verify_roe_equal_state();
void verify_roe_symmetry();
void verify_roe_stationary_contact();
void verify_roe_supersonic_upwind();
void verify_roe_supersonic_upwind_transverse();
void verify_roe_entropy_fix_rarefactions();
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
void verify_boundary_constructors(bool output = true);
//Outflow
void verify_boundary_outflow(bool output = true);
    void verify_boundary_outflow_1D();
    void verify_boundary_outflow_2D();
    void verify_boundary_outflow_3D();
    void verify_boundary_outflow_1D_gated();
    void verify_boundary_outflow_2D_gated();
    void verify_boundary_outflow_3D_gated_X();
    void verify_boundary_outflow_3D_gated_Y();
    void verify_boundary_outflow_3D_gated_Z();
//Reflective
void verify_boundary_reflective(bool output = true);
    void verify_boundary_reflective_1D();
    void verify_boundary_reflective_2D();
    void verify_boundary_reflective_3D();
//Periodic
void verify_boundary_periodic(bool output = true);
    void verify_boundary_periodic_1D();
    void verify_boundary_periodic_2D();
    void verify_boundary_periodic_3D();
//Fixed
void verify_boundary_fixed(bool output = true);
    void verify_boundary_fixed_1D();
    void verify_boundary_fixed_2D();
    void verify_boundary_fixed_3D();
//Ignore
void verify_boundary_ignore(bool output = true);
    void verify_boundary_ignore_1D();
    void verify_boundary_ignore_2D();
    void verify_boundary_ignore_3D();
    void verify_boundary_ignore_blocks_implicit_outflow();
//Composition
void verify_boundary_composition(bool output = true);
    void verify_boundary_set_missing_faces_outflow_2D();
    void verify_boundary_composition_3D();
    void verify_boundary_composition_order();
    void verify_boundary_composition_overlap_order();


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
void verify_cfl_mhd_speed_3D();
void verify_cfl_mhd_time_3D();
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


//MARK: TVD
void verify_tvd(bool output = true);
void verify_tvd_scalar_minmod();
void verify_tvd_scalar_mc();
void verify_tvd_scalar_vanleer();
void verify_tvd_scalar_superbee();
void verify_tvd_scalar_vanalbada();
void verify_tvd_vec3_limiters();
void verify_tvd_primitive_limiters();
void verify_tvd_limiter_dispatch();
void verify_tvd_muscl_constant_state();
void verify_tvd_muscl_zero_dt_spatial_reconstruction();
void verify_tvd_muscl_nonzero_dt_predictor();
void verify_tvd_muscl_falls_back_when_unphysical();


//MARK: Godunov Scheme
//1D
void verify_godunov_1D(bool output = true);
    void verify_god_uniform_stationary_1D();
    void verify_god_uniform_moving_1D();
    void verify_god_periodic_conservation_1D();
    void verify_god_dist_grid_1D();
    void verify_god_dt0_1D();
//2D
void verify_godunov_2D_Split(bool output = true);
void verify_godunov_2D_Unsplit(bool output = true);
    void verify_god_uniform_stationary_2D(bool split);
    void verify_god_uniform_moving_2D(bool split);
    void verify_god_periodic_conservation_2D(bool split);
    void verify_god_dist_grid_2D();
    void verify_god_dt0_2D(bool split);
    void verify_2D_X_match_1D(bool split);
    void verify_2D_Y_match_1D(bool split);
    void verify_ctu_diagonal_contact_2D();
    void verify_ctu_blast_2D();
//3D
void verify_godunov_3D_Split(bool output = true);
void verify_godunov_3D_Unsplit(bool output = true);
    void verify_god_uniform_stationary_3D(bool split);
    void verify_god_uniform_moving_3D(bool split);
    void verify_god_periodic_conservation_3D(bool split);
    void verify_god_dist_grid_3D();
    void verify_god_dt0_3D(bool split);
    void verify_3D_X_match_1D(bool split);
    void verify_3D_Y_match_1D(bool split);
    void verify_3D_Z_match_1D(bool split);
    void verify_ctu_diagonal_contact_3D();
    void verify_ctu_blast_3D();

#ifdef MHD
//MARK: Godunov Scheme (MHD)
//1D
void verify_godunov_1D_MHD(bool output = true);
    void verify_god_uniform_stationary_1D_MHD();
    void verify_god_uniform_moving_1D_MHD();
    void verify_god_periodic_conservation_1D_MHD();
//2D
void verify_godunov_2D_MHD(bool output = true);
    void verify_god_uniform_stationary_2D_MHD();
    void verify_god_uniform_moving_2D_MHD();
    void verify_god_periodic_conservation_2D_MHD();
    void verify_ct_stationary_2D();
    void verify_ct_uniform_Bz0_2D();
void verify_god_dist_grid_2D_MHD();
void verify_godunov_2D_MHD_Match(bool output = true);
    void verify_2D_X_match_1D_MHD();
    void verify_2D_Y_match_1D_MHD();

//3D
void verify_godunov_3D_MHD(bool output = true);
    void verify_god_uniform_stationary_3D_MHD();
    void verify_god_uniform_moving_3D_MHD();
    void verify_god_periodic_conservation_3D_MHD();
    void verify_ct_stationary_3D();
    void verify_ct_uniform_Bz0_3D();
void verify_god_dist_grid_3D_MHD();
void verify_godunov_3D_MHD_Match(bool output = true);
    void verify_3D_X_match_1D_MHD();
    void verify_3D_Y_match_1D_MHD();
    void verify_3D_Z_match_1D_MHD();


//MARK: Constrained Transport
//2D
void verify_ct_2D(bool output = true);
//Faraday Update
    void verify_ct_divergence_2D();
        void assert_divergenceless(const ExtendedArray2D<vec3>& B, double dx, double dy);
    void verify_ct_stokes_theorem_2D();
    void verify_ct_uniform_E_2D();
    //Face and Body Field Machinery
    void verify_ct_body_fields_2D();
    void verify_ct_copy_face_fields_2D();
//3D
void verify_ct_3D(bool output = true);
//Faraday Update
    void verify_ct_divergence_3D();
        void assert_divergenceless(const ExtendedArray3D<vec3>& B, double dx, double dy, double dz);
    void verify_ct_stokes_theorem_3D();
    void verify_ct_uniform_E_3D();
//Face and Body Field Machinery
    void verify_ct_body_fields_3D();
    void verify_ct_copy_face_fields_3D();


#endif

//MARK: Passive Scalars
void verify_passive_scalars(bool output = true);
//Key Management
void verify_passive_scalar_add();
void verify_passive_scalar_remove();
void verify_passive_scalar_no_op();
void verify_passive_scalar_bad_lookup();
//1D/2D/3D add/remove
void verify_passive_scalar_add_1D();
void verify_passive_scalar_remove_1D();
void verify_passive_scalar_add_2D();
void verify_passive_scalar_remove_2D();
void verify_passive_scalar_add_3D();
void verify_passive_scalar_remove_3D();
//Advection
void verify_passive_scalar_zero_flux_1D();
void verify_passive_scalar_zero_flux_2D();
void verify_passive_scalar_zero_flux_3D();
void verify_passive_scalar_uniform_1D();
void verify_passive_scalar_uniform_2D();
void verify_passive_scalar_uniform_3D();
void verify_passive_scalar_conserve_1D();
void verify_passive_scalar_conserve_2D();
void verify_passive_scalar_conserve_3D();
void verify_passive_scalar_dir_1D();
void verify_passive_scalar_dir_2D();
void verify_passive_scalar_dir_3D();

//MARK: DRAGONHOARD
void verify_HOARD(bool output = true);
void verify_IO1D();
void verify_IO2D();
void verify_IO3D();
void verify_IO_dim_assert();
void verify_IO_size_assert();


//MARK: DRAGONWING
void verify_WING(bool output = true);
//Memory Management
void verify_WING_alloc_prim_dimensions();
void verify_WING_alloc_flux_dimensions();
void verify_WING_alloc_vec3_dimensions();
void verify_WING_array_reuse();
void verify_WING_arrayguard_get();
void verify_WING_double_release();
void verify_WING_scope_release();
void verify_WING_array_purge_survival();
//Atomics
void verify_WING_fbcount_no_op();
void verify_WING_fbcount_reset();
void verify_WING_fbcount_single_no_throw();
void verify_WING_fbcount_single_over_throw();
void verify_WING_fbcount_cumulative_throw();
//Multithreading
void verify_WING_restart_message();
void verify_WING_launchParallel_calls_advance();
void verify_WING_launchParallel_catches();
void verify_WING_checkpoints();
void verify_WING_checkpoints_restart();
}


#endif /* Testing_hpp */
