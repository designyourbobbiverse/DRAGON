//
//  Sources.hpp
//  DRAGON/Source
//  User-Facing Header file
//
//  Created by Bobbie Markwick on 22/09/2026.
//

#ifndef UniformGravity_hpp
#define UniformGravity_hpp


#include "Source/Source.hpp"


namespace DRAGON::Source{

class UniformGravity: public MomentumSource {
public:
    vec3 g; //Acceleration due to gravity
    UniformGravity(vec3 g_): g(g_) {}
    UniformGravity(double gx, double gy, double gz): g({gx, gy, gz}) {}

    vec3 force(const PrimitiveState& w, double t) const override;
};

}


#endif
