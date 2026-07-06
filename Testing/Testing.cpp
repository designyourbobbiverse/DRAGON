//
//  Testing.cpp
//  DRAGON/Testing
//
//  Created by Bobbie Markwick on 12/05/2026.
//
#include "Testing.hpp"
#include <iostream>

using namespace DRAGON_Test;


//Testing TODOs:
//AMRGrid
//DRAGONWING

int main() {
    verify_fluid_element();
    verify_riemann();
    verify_grid();
    verify_boundary();
    verify_tvd();
    verify_cfl();
    
    verify_IO();

    verify_godunov_1D();
    verify_godunov_2D_Split();
    verify_godunov_2D_Unsplit();
    verify_godunov_3D_Split();
    verify_godunov_3D_Unsplit();
    std::cout << "All Godunov Scheme tests passed.\n\n";
    
    #ifdef MHD
    verify_ct_2D();
    verify_ct_3D();
    std::cout << "All Constrained Transport Tests passed.\n\n";
    #endif
    
    std::cout << "All tests passed.\n";

    
    
}
