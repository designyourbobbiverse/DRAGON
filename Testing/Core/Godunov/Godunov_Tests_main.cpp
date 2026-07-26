//
//  Godunov_Tests_main.cpp
//  DRAGON/Testing/Core/Godunov
//
//  Created by Bobbie Markwick on 09/06/2026.
//

#include "Testing.hpp"
#include "Boundary.hpp"
#include "Grid.hpp"
#include "CFL.hpp"
#include <iostream>

using namespace DRAGON_Test;

//MARK: 1D
void DRAGON_Test::verify_godunov_1D(bool output){
    int prev = CONFIG::riemann_choice;
    CONFIG::riemann_choice = RIEMANN_EXACT;
    
    if(output) std::cout<<"Godunov Scheme (1D): \n";
    if(output) std::cout<<"- Zero Time: ";
    verify_god_dt0_1D();
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- Uniform Flows: ";
    verify_god_uniform_stationary_1D();
    verify_god_uniform_moving_1D();
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- Periodic Conservation: ";
    verify_god_periodic_conservation_1D();
    if(output) std::cout<<"Passed\n";
    #ifdef MHD
    verify_godunov_1D_MHD();
    #endif
    if(output) std::cout<<"- Parallelization: ";
    verify_god_dist_grid_1D();
    if(output) std::cout<<"Passed\n";

    CONFIG::riemann_choice = prev;
}

//MARK: 2D Split
void DRAGON_Test::verify_godunov_2D_Split(bool output){
    int prev = CONFIG::riemann_choice;
    CONFIG::riemann_choice = RIEMANN_EXACT;

    if(output) std::cout<<"Split Scheme (2D): \n";
    if(output) std::cout<<"- Zero Time: ";
    verify_god_dt0_2D(true);
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- Uniform Flows: ";
    verify_god_uniform_stationary_2D(true);
    verify_god_uniform_moving_2D(true);
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- Periodic Conservation: ";
    verify_god_periodic_conservation_2D(true);
    if(output) std::cout<<"Passed\n";
    #ifndef DIMENSION_UNSPLIT
    if(output) std::cout<<"- Parallelization: ";
    verify_god_dist_grid_2D();
    if(output) std::cout<<"Passed\n";
    #endif
    if(output) std::cout<<"- 1D Match (X): ";
    verify_2D_X_match_1D(true);
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- 1D Match (Y): ";
    verify_2D_Y_match_1D(true);
    if(output) std::cout<<"Passed\n";
    
    CONFIG::riemann_choice = prev;
}

//MARK: 2D Unsplit
void DRAGON_Test::verify_godunov_2D_Unsplit(bool output){
    int prev = CONFIG::riemann_choice;
    CONFIG::riemann_choice = RIEMANN_EXACT;

    if(output) std::cout<<"Unsplit Scheme (2D): \n";

    if(output) std::cout<<"- Zero Time: ";
    verify_god_dt0_2D(false);
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- Uniform Flows: ";
    verify_god_uniform_stationary_2D(false);
    
    
    verify_god_uniform_moving_2D(false);
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- Periodic Conservation: ";
    verify_god_periodic_conservation_2D(false);
    if(output) std::cout<<"Passed\n";
    #ifdef MHD
    verify_godunov_2D_MHD();
    #endif

    #ifdef DIMENSION_UNSPLIT
    if(output) std::cout<<"- Parallelization: ";
    verify_god_dist_grid_2D();
    #ifdef MHD
    verify_god_dist_grid_2D_MHD();
    #endif
    if(output) std::cout<<"Passed\n";
    #endif
    
    #ifdef CTU
    if(output) std::cout<<"- Diagonal Contact: ";
    verify_ctu_diagonal_contact_2D();
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- Blast: ";
    verify_ctu_blast_2D();
    if(output) std::cout<<"Passed\n";
    #endif
    
    if(output) std::cout<<"- 1D Match (X): ";
    verify_2D_X_match_1D(false);
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- 1D Match (Y): ";
    verify_2D_Y_match_1D(false);
    if(output) std::cout<<"Passed\n";
#ifdef MHD
    if(output) std::cout<<"- 1D Match (MHD-X): ";
    verify_2D_X_match_1D_MHD();
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- 1D Match (MHD-Y): ";
    verify_2D_Y_match_1D_MHD();
    if(output) std::cout<<"Passed\n";
#endif
    CONFIG::riemann_choice = prev;
}

//MARK: 3D Split
void DRAGON_Test::verify_godunov_3D_Split(bool output){
    int prev = CONFIG::riemann_choice;
    CONFIG::riemann_choice = RIEMANN_EXACT;

    if(output) std::cout<<"Split Scheme (3D): \n";
    if(output) std::cout<<"- Zero Time: ";
    verify_god_dt0_3D(true);
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- Uniform Flows: ";
    verify_god_uniform_stationary_3D(true);
    verify_god_uniform_moving_3D(true);
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- Periodic Conservation: ";
    verify_god_periodic_conservation_3D(true);
    if(output) std::cout<<"Passed\n";
    #ifndef DIMENSION_UNSPLIT
    if(output) std::cout<<"- Parallelization: ";
    verify_god_dist_grid_3D();
    if(output) std::cout<<"Passed\n";
    #endif
    
    if(output) std::cout<<"- 1D Match (X): ";
    verify_3D_X_match_1D(true);
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- 1D Match (Y): ";
    verify_3D_Y_match_1D(true);
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- 1D Match (Z): ";
    verify_3D_Z_match_1D(true);
    if(output) std::cout<<"Passed\n";
    
    CONFIG::riemann_choice = prev;
}
//MARK: 3D Unsplit
void DRAGON_Test::verify_godunov_3D_Unsplit(bool output){
    int prev = CONFIG::riemann_choice;
    CONFIG::riemann_choice = RIEMANN_HLLC;

    if(output) std::cout<<"Unsplit Scheme (3D): \n";
    if(output) std::cout<<"- Zero Time: ";
    verify_god_dt0_3D(false);
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- Uniform Flows: ";
    verify_god_uniform_stationary_3D(false);
    verify_god_uniform_moving_3D(false);
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- Periodic Conservation: ";
    verify_god_periodic_conservation_3D(false);
    if(output) std::cout<<"Passed\n";

    #ifdef MHD
    verify_godunov_3D_MHD();
    #endif

    #ifdef DIMENSION_UNSPLIT
    if(output) std::cout<<"- Parallelization: ";
    verify_god_dist_grid_3D();
    #ifdef MHD
    verify_god_dist_grid_3D_MHD();
    #endif
    if(output) std::cout<<"Passed\n";
    #endif
   
    #ifdef CTU
    if(output) std::cout<<"- Diagonal Contact: ";
    verify_ctu_diagonal_contact_3D();
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- Blast: ";
    verify_ctu_blast_3D();
    if(output) std::cout<<"Passed\n";
    #endif
    
    if(output) std::cout<<"- 1D Match (X): ";
    verify_3D_X_match_1D(false);
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- 1D Match (Y): ";
    verify_3D_Y_match_1D(false);
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- 1D Match (Z): ";
    verify_3D_Z_match_1D(false);
    if(output) std::cout<<"Passed\n";
    #ifdef MHD
    if(output) std::cout<<"- 1D Match (MHD-X): ";
    verify_3D_X_match_1D_MHD();
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- 1D Match (MHD-Y): ";
    verify_3D_Y_match_1D_MHD();
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- 1D Match (MHD-Z): ";
    verify_3D_Z_match_1D_MHD();
    if(output) std::cout<<"Passed\n";
    #endif
    

    CONFIG::riemann_choice = prev;
}
