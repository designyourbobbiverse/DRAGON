
//
//  B_fields.cpp
//  DRAGON/MHD
//
//  Created by Bobbie Markwick on 16/07/2026.
//  CT Implementaiton based in part on https://arxiv.org/abs/0710.0424
//      Evans and Hawley (1988). https://doi.org/10.1086/166684
//      Gardiner and Stone (2005 & 2008). https://arxiv.org/abs/astro-ph/0501557   https://arxiv.org/abs/0712.2634


#include "MHD/CT.hpp"
#include "Hydro/Grid.hpp" //Implements Grid::initialize_B_fields()

#include "Config.h"
#include "Constants.h" //For _1_8pi
using namespace DRAGON;

#ifdef MHD

//MARK: Faraday's Law
//Advances B by -curl(E)*dt
void CT::Faraday(const MagneticArray2D& E, MagneticArray2D& _B, double dt_dx, double dt_dy, int g){
    
    const int nx = _B.getSizeX()-1, ny = _B.getSizeY()-1;
    
    for (int i=-g; i<=nx+g; i++) {
        for (int j=-g; j<ny+g; j++) {
            _B[i,j].x -= (E[i,j+1].z - E[i,j].z) * dt_dy;
        }
    }
    for (int i=-g; i<nx+g; i++) {
        for (int j=-g; j<=ny+g; j++) {
            _B[i,j].y -= (E[i,j].z - E[i+1,j].z) * dt_dx;
        }
    }
    for (int i=-g; i<nx+g; i++) {
        for (int j=-g; j<ny+g; j++) {
            _B[i,j].z -= (E[i+1,j].y - E[i,j].y) * dt_dx -  (E[i,j+1].x - E[i,j].x) * dt_dy;
        }
    }
}
void CT::Faraday(const MagneticArray3D& E, MagneticArray3D& _B, double dt_dx, double dt_dy, double dt_dz, int g){
    
    const int nx = _B.getSizeX()-1, ny = _B.getSizeY()-1, nz = _B.getSizeZ()-1;
    
    for (int i=-g; i<=nx+g; i++) {
        for (int j=-g; j<ny+g; j++) {
            for (int k=-g; k<nz+g; k++) {
                _B[i,j,k].x -= (E[i,j+1,k].z - E[i,j,k].z) * dt_dy - (E[i,j,k+1].y - E[i,j,k].y) * dt_dz;
            }
        }
    }
    for (int i=-g; i<nx+g; i++) {
        for (int j=-g; j<=ny+g; j++) {
            for (int k=-g; k<nz+g; k++) {
                _B[i,j,k].y -= (E[i,j,k+1].x - E[i,j,k].x) * dt_dz - (E[i+1,j,k].z - E[i,j,k].z) * dt_dx;
            }
        }
    }
    for (int i=-g; i<nx+g; i++) {
        for (int j=-g; j<ny+g; j++) {
            for (int k=-g; k<=nz+g; k++) {
                _B[i,j,k].z -= (E[i+1,j,k].y - E[i,j,k].y) * dt_dx -  (E[i,j+1,k].x - E[i,j,k].x) * dt_dy;
            }
        }
    }
}
//MARK: Body Fields
//Body fields are an average of the adjacent face states
namespace DRAGON::CT{
void computeBodyField(const MagneticArray2D& B, FluidArray2D& w, int i, int j, bool consv_E){
    const vec3 _B =  {(B[i,j].x + B[i+1,j].x)*0.5, (B[i,j].y + B[i,j+1].y)*0.5, B[i,j].z};
    if (consv_E) {
        ConservativeState U(w[i,j]);//Update the conservative element to keep quantities conserved
        U.B = _B;
        w[i,j] = U;
    } else { //Just update the magnetic field, keeps thermal energy where it was (numerically safer, but less physical)
        w[i,j].B = _B;
    }
}
void computeBodyField(const MagneticArray3D& B, FluidArray3D& w, int i, int j, int k, bool consv_E){
    const vec3 _B = (B[i,j,k] + vec3{B[i+1,j,k].x, B[i,j+1,k].y, B[i,j,k+1].z}) * 0.5;
    if (consv_E) {
        ConservativeState U(w[i,j,k]);//Update the conservative element to keep quantities conserved
        U.B = _B;
        w[i,j,k] = U;
    } else { //Just update the magnetic field, keeps thermal energy where it was (numerically safer, but less physical)
        w[i,j,k].B = _B;
    }
}

//Low-beta CT update has a tendancy to drive pressure below zero, which is unphysical
//One way to avoid this is to keep thermal energy (instead of total energy) fixed during the CT update.
//This is a numerical-stability vs physical-realism tradeoff, controlled in Config.h
bool shouldProtectThermal(const PrimitiveState& w){
    #if CT_ENERGY_CONSV == CHOOSE_RUNTIME || defined(TESTMODE)
    switch (Config::CT_energy_choice) {
        case CT_CONSV_TOTAL_E: return false;
        case CT_CONSV_THERMAL: return true;
        case CT_CONSV_BETA_GATED: {
            double B2 = (w.B * w.B);
            if (B2 == 0) return false;
            double beta = w.p / (_1_8pi * w.B*w.B);
            return beta <= Config::ct_energy_beta;
        }
        default: return false;
    }
    #elif CT_ENERGY_CONSV == CT_CONSV_TOTAL_E
    return false;
    #elif CT_ENERGY_CONSV == CT_CONSV_THERMAL
    return true;
    #elif CT_CONSV_BETA_GATED
    double B2 = (w.B * w.B);
    if (B2 == 0) return false;
    double beta = 8*M_PI*w.p / (w.B*w.B);
    return beta <= Config::ct_energy_beta;
    #endif
}
}

void Grid2D::initialize_B_fields(){
    const int nx = w.getSizeX(), ny = w.getSizeY(), ng = w.getGhosts();
    boundary.apply(*this); //Apply boundary conditions before initializing B
    
    for (int i=-ng; i<nx+ng; i++) {
        for (int j=-ng; j<ny+ng; j++) {
            CT::computeBodyField(B,w,i,j,false); //When initializing, always preserve the user's initial thermal pressure
        }
    }
}

void Grid2D::computeBodyAveragedFields(const MagneticArray2D& B){ CT::computeBodyFields(B, w); }
void CT::computeBodyFields(const MagneticArray2D& B, FluidArray2D& w){
    const int nx = w.getSizeX(), ny = w.getSizeY(), ng = w.getGhosts();
    for (int i=-ng; i<nx+ng; i++) {
        for (int j=-ng; j<ny+ng; j++) {
            computeBodyField(B, w, i, j, !shouldProtectThermal(w[i,j]));
        }
    }
}

void Grid3D::initialize_B_fields(){
    const int nx = w.getSizeX(), ny = w.getSizeY(), nz = w.getSizeZ(), ng = w.getGhosts();
    boundary.apply(*this); //Apply boundary conditions before initializing B
    
    for (int i=-ng; i<nx+ng; i++) {
        for (int j=-ng; j<ny+ng; j++) {
            for (int k=-ng; k<nz+ng; k++) {
                CT::computeBodyField(B,w,i,j,k,false); //When initializing, always preserve the user's initial thermal pressure
            }
        }
    }
}

void Grid3D::computeBodyAveragedFields(const MagneticArray3D& B){ CT::computeBodyFields(B, w); }
void CT::computeBodyFields(const MagneticArray3D& B, FluidArray3D& w){
    const int nx = w.getSizeX(), ny = w.getSizeY(), nz = w.getSizeZ(), ng = w.getGhosts();
    for (int i=-ng; i<nx+ng; i++) {
        for (int j=-ng; j<ny+ng; j++) {
            for (int k=-ng; k<nz+ng; k++) {
                computeBodyField(B, w, i, j, k, !shouldProtectThermal(w[i,j,k]));
            }
        }
    }
}



//MARK: Face-Centred Field Copy
void CT::copyFaceFields_X( FluidArray2D& _L,const MagneticArray2D& _B, FluidArray2D& _R){
    const int nx = _L.getSizeX(), ny = _L.getSizeY(), g = _B.getGhosts();
    for (int i=-g; i<nx+g; i++) {
        for (int j=-g; j<ny+g; j++) {
            _L[i,j].B.x = _B[i,j].x;
            _R[i,j].B.x = _B[i+1,j].x;
        }
    }
}
void CT::copyFaceFields_Y( FluidArray2D& _L,const MagneticArray2D& _B, FluidArray2D& _R){
    const int nx = _L.getSizeX(), ny = _L.getSizeY(), g = _B.getGhosts();
    for (int i=-g; i<nx+g; i++) {
        for (int j=-g; j<ny+g; j++) {
            _L[i,j].B.y = _B[i,j].y;
            _R[i,j].B.y = _B[i,j+1].y;
        }
    }
}
void CT::copyFaceFields_X( FluidArray3D& _L,const MagneticArray3D& _B, FluidArray3D& _R){
    const int nx = _L.getSizeX(), ny = _L.getSizeY(), nz = _L.getSizeZ(), g = _B.getGhosts();
    for (int i=-g; i<nx+g; i++) {
        for (int j=-g; j<ny+g; j++) {
            for (int k=-g; k<nz+g; k++) {
                _L[i,j,k].B.x = _B[i,j,k].x;
                _R[i,j,k].B.x = _B[i+1,j,k].x;
            }
        }
    }
}
void CT::copyFaceFields_Y( FluidArray3D& _L,const MagneticArray3D& _B, FluidArray3D& _R){
    const int nx = _L.getSizeX(), ny = _L.getSizeY(), nz = _L.getSizeZ(), g = _B.getGhosts();
    for (int i=-g; i<nx+g; i++) {
        for (int j=-g; j<ny+g; j++) {
            for (int k=-g; k<nz+g; k++) {
                _L[i,j,k].B.y = _B[i,j,k].y;
                _R[i,j,k].B.y = _B[i,j+1,k].y;
            }
        }
    }
}
void CT::copyFaceFields_Z( FluidArray3D& _L,const MagneticArray3D& _B, FluidArray3D& _R){
    const int nx = _L.getSizeX(), ny = _L.getSizeY(), nz = _L.getSizeZ(), g = _B.getGhosts();
    for (int i=-g; i<nx+g; i++) {
        for (int j=-g; j<ny+g; j++) {
            for (int k=-g; k<nz+g; k++) {
                _L[i,j,k].B.z = _B[i,j,k].z;
                _R[i,j,k].B.z = _B[i,j,k+1].z;
            }
        }
    }
}

#endif


