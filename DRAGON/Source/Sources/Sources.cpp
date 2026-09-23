//
//  Sources.cpp
//  DRAGON/Source/Sources
//
//  Created by Bobbie Markwick on 04/09/2026.
//

#include "Sources.hpp"

using namespace DRAGON;
using namespace DRAGON::Source;


vec3 UniformGravity::force(const PrimitiveState& w, double t) const {
    return w.rho * g;
}
