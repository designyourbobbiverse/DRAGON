//
//  TVD.hpp
//  DRAGON
//
//  Created by Bobbie Markwick on 10/06/2026.
//
#include "FluidElement.hpp"

namespace TVD{

    void MUSCL(const PrimitiveState& wL, PrimitiveState& _L, const PrimitiveState& wC, PrimitiveState& _R, const PrimitiveState& wR, double dt_dL);

//Limiter Functions
    double minmod(double a, double b);
    PrimitiveState minmod(const PrimitiveState& a, const PrimitiveState& b);

}

