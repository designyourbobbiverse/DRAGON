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
//CFL with MHD
//HLLD-specific tests with MHD
//Add magnetic field checks to boundary tests

int main() {
    
    verify_fluid_element();
    verify_riemann();
    verify_grid();
    verify_boundary();
    verify_tvd();
    verify_cfl();

    verify_godunov_1D();
    verify_godunov_2D_Split();
    verify_godunov_2D_Unsplit();
    verify_godunov_3D_Split();
    verify_godunov_3D_Unsplit();
    std::cout << "All Godunov Scheme tests passed.\n\n";

    
    std::cout << "All tests passed.\n";

    
    
}
