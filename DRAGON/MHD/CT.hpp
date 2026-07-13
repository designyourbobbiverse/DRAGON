//
//  CT.hpp
//  DRAGON/MHD
//
//  Created by Bobbie Markwick on 24/06/2026.
//
#include "Grid.hpp"
#include "Config.h"

typedef ExtendedArray2D<PrimitiveState> FluidArray2D;
typedef ExtendedArray2D<ConservativeState> FluxArray2D;
typedef ExtendedArray3D<PrimitiveState> FluidArray3D;
typedef ExtendedArray3D<ConservativeState> FluxArray3D;

#ifdef MHD
typedef ExtendedArray2D<vec3> MagneticArray2D;
typedef ExtendedArray3D<vec3> MagneticArray3D;

namespace CT {
    void computeBodyFields(const MagneticArray2D& B, FluidArray2D& w);
    void computeFaceFields(const MagneticArray2D& _A, MagneticArray2D& _B, double dx, double dy);
    void copyFaceFields_X( FluidArray2D& _L,const MagneticArray2D& _B, FluidArray2D& _R);
    void copyFaceFields_Y( FluidArray2D& _L,const MagneticArray2D& _B, FluidArray2D& _R);
    void computeElectric(MagneticArray2D& E, const FluxArray2D& F_X,const FluxArray2D& F_Y, int ghosts=0);
    void updatePotential(MagneticArray2D& _A, const MagneticArray2D& E, double dt, int ghosts=0);

    void computeBodyFields(const MagneticArray3D& B, FluidArray3D& w);
    void computeFaceFields(const MagneticArray3D& _A, MagneticArray3D& _B, double dx, double dy, double dz);
    void copyFaceFields_X( FluidArray3D& _L,const MagneticArray3D& _B, FluidArray3D& _R);
    void copyFaceFields_Y( FluidArray3D& _L,const MagneticArray3D& _B, FluidArray3D& _R);
    void copyFaceFields_Z( FluidArray3D& _L,const MagneticArray3D& _B, FluidArray3D& _R);
    void computeElectric(MagneticArray3D& E, const FluxArray3D& F_X,const FluxArray3D& F_Y,const FluxArray3D& F_Z, int ghosts=0);
    void updatePotential(MagneticArray3D& _A, const MagneticArray3D& E, double dt, int ghosts=0);

}
#endif
