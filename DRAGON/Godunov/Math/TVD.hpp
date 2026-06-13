//
//  TVD.hpp
//  DRAGON
//
//  Created by Bobbie Markwick on 10/06/2026.
//
#include "FluidElement.hpp"

namespace TVD{

    void MUSCL(const PrimitiveState& wL, PrimitiveState& _L, const PrimitiveState& wC, PrimitiveState& _R, const PrimitiveState& wR, double dt_dL);

    PrimitiveState limit(const PrimitiveState& a, const PrimitiveState& b);

    //Limiter Functions
    double minmod(double a, double b);
    PrimitiveState minmod(const PrimitiveState& a, const PrimitiveState& b);

    double MC(double a, double b);
    PrimitiveState MC(const PrimitiveState& a, const PrimitiveState& b);

    double vanLeer(double a, double b);
    PrimitiveState vanLeer(const PrimitiveState& a, const PrimitiveState& b);

    double superbee(double a, double b);
    PrimitiveState superbee(const PrimitiveState& a, const PrimitiveState& b);

    double vanAlbada(double a, double b);
    PrimitiveState vanAlbada(const PrimitiveState& a, const PrimitiveState& b);

}

