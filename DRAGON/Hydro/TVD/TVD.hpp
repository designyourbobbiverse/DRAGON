//
//  TVD.hpp
//  DRAGON/Hydro/TVD
//
//  Created by Bobbie Markwick on 10/06/2026.
//
#include "FluidElement.hpp"

namespace TVD{
    //MUSCL Reconstruction
    void MUSCL(const PrimitiveState& wL, PrimitiveState& _L, const PrimitiveState& wC, PrimitiveState& _R, const PrimitiveState& wR, double dt_dL);

    //Apply the limiter chosen in Config.h
    PrimitiveState limit(const PrimitiveState& a, const PrimitiveState& b);

    //Limiter Functions
    double minmod(double a, double b);
    vec3 minmod(const vec3& a, const vec3& b);
    PrimitiveState minmod(const PrimitiveState& a, const PrimitiveState& b);

    double MC(double a, double b);
    vec3 MC(const vec3& a, const vec3& b);
    PrimitiveState MC(const PrimitiveState& a, const PrimitiveState& b);

    double vanLeer(double a, double b);
    vec3 vanLeer(const vec3& a, const vec3& b);
    PrimitiveState vanLeer(const PrimitiveState& a, const PrimitiveState& b);

    double superbee(double a, double b);
    vec3 superbee(const vec3& a, const vec3& b);
    PrimitiveState superbee(const PrimitiveState& a, const PrimitiveState& b);

    double vanAlbada(double a, double b);
    vec3 vanAlbada(const vec3& a, const vec3& b);
    PrimitiveState vanAlbada(const PrimitiveState& a, const PrimitiveState& b);

}
