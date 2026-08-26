//
//  Hydro/ExtendedArray/ArrayTypes.hpp
//  DRAGON/Hydro/ExtendedArray
//
//  Created by Bobbie Markwick on 31/07/2026.
//

#ifndef Array_Types_hpp
#define Array_Types_hpp

#include "Hydro/ExtendedArray/ExtendedArray.hpp" //For ExtendedArray
#include "FluidElement/FluidElement.hpp"         //For PrimitiveState, ConservativeState, vec3

namespace DRAGON{
typedef ExtendedArray1D<PrimitiveState> FluidArray1D;
typedef ExtendedArray2D<PrimitiveState> FluidArray2D;
typedef ExtendedArray3D<PrimitiveState> FluidArray3D;

typedef ExtendedArray2D<ConservativeState> FluxArray2D;
typedef ExtendedArray3D<ConservativeState> FluxArray3D;

typedef ExtendedArray2D<vec3> MagneticArray2D;
typedef ExtendedArray3D<vec3> MagneticArray3D;
}

#endif
