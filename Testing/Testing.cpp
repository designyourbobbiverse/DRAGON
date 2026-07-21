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
    
    auto start = std::chrono::system_clock::now();

    verify_fluid_element();
    
    verify_grid();
    verify_boundary();
   // verify_IO();
    
    verify_riemann();
    verify_tvd();
    verify_cfl();

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
    
    double clock_time = std::chrono::duration<double>(std::chrono::system_clock::now() - start).count();
    int m = floor((clock_time)/60.0);
    double s = round((clock_time - m*60)*100)/100.0;
    std::cout << "Time: "<< m <<"m " << s <<"s \n";


    
    
}
