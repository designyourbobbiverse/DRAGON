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
    verify_fluid_element();
    verify_riemann();
    verify_grid();
    verify_boundary();
    verify_cfl();
    
    verify_godunov_1D();
    verify_godunov_2D_Split();
    verify_godunov_3D_Split();
    std::cout << "All tests passed.\n";

    
    
}
