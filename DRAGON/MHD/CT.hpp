//
//  CT.hpp
//  DRAGON/MHD
//
//  Created by Bobbie Markwick on 24/06/2026.
//
#ifndef __CT_hpp__
#define __CT_hpp__

#include "ArrayTypes.hpp"
#include "Config.h"

#ifdef MHD
namespace DRAGON::CT {

//Electric Fields
//2D
void computeElectric(MagneticArray2D& E, const FluxArray2D& F_X,const FluxArray2D& F_Y, int ghosts=0);
void upwindElectric(MagneticArray2D& E, const FluxArray2D& F_X,const FluxArray2D& F_Y, const MagneticArray2D& Eref, int ghosts=0);
void bodyElectric(const FluidArray2D& w, MagneticArray2D& E, int ghosts=0); //1001001 SOS (iykyk)
//3D
void computeElectric(MagneticArray3D& E, const FluxArray3D& F_X,const FluxArray3D& F_Y,const FluxArray3D& F_Z, int ghosts=0);
void upwindElectric(MagneticArray3D& E, const FluxArray3D& F_X,const FluxArray3D& F_Y,const FluxArray3D& F_Z, const MagneticArray3D& Eref, int ghosts=0);
void bodyElectric(const FluidArray3D& w, MagneticArray3D& E, int ghosts=0); //1001001 SOS (iykyk)

//Faraday: Advance B by -curl(E)*dt
//2D
void Faraday(const MagneticArray2D& E, MagneticArray2D& _B, double dt_dx, double dt_dy, int ghosts = 0);
//3D
void Faraday(const MagneticArray3D& E, MagneticArray3D& _B, double dt_dx, double dt_dy, double dt_dz, int ghosts = 0);


//Face Fields -> Body fields or Interface States
//2D
void computeBodyFields(const MagneticArray2D& B, FluidArray2D& w);
void copyFaceFields_X( FluidArray2D& _L,const MagneticArray2D& _B, FluidArray2D& _R);
void copyFaceFields_Y( FluidArray2D& _L,const MagneticArray2D& _B, FluidArray2D& _R);
//3D
void computeBodyFields(const MagneticArray3D& B, FluidArray3D& w);
void copyFaceFields_X( FluidArray3D& _L,const MagneticArray3D& _B, FluidArray3D& _R);
void copyFaceFields_Y( FluidArray3D& _L,const MagneticArray3D& _B, FluidArray3D& _R);
void copyFaceFields_Z( FluidArray3D& _L,const MagneticArray3D& _B, FluidArray3D& _R);

}
#endif
#endif
