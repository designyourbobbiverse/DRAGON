//
//  Testing.hpp
//  Hydro Code
//
//  Created by Bobbie Markwick on 12/05/2026.
//

#ifndef Testing_hpp
#define Testing_hpp

#include <stdio.h>
#include <cassert>
#include "FluidElement.hpp"

namespace DRAGON_Test{

void verify_fluid_element();

//Numerical Close-ness
bool approx(double a, double b, double rel = 1e-12, double abs = 1e-12);
void expect_close(const ConservativeState& a, const ConservativeState& b,  double rel = 1e-12, double abs = 1e-12);
void expect_close(const PrimitiveState& a, const PrimitiveState& b, double rel = 1e-12, double abs = 1e-12);

//Fluid Element Data Structure
void verify_constructors();
void verify_conversion();
void verify_flux();
void verify_enthalpy();
void verify_swaps_P();//Primitive
void verify_swaps_C();//Conservative
//Fluid Element Arithmetic
void verify_add();
void verify_sub();
void verify_mult();
void verify_div();


}


#endif /* Testing_hpp */
