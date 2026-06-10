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
    
    std::cout << "All tests passed.\n";

    
    
}
