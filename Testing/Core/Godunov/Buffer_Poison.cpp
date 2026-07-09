//
//  Buffer_Poison.cpp
//  DRAGON/Testing/Core/Godunov
//
//  Created by Bobbie Markwick on 9/07/2026.
//

#include "GridBuffers.hpp"
#include "FluidElement.hpp"
#include "Config.h"
#include <cmath>


void GridBuffers1D::poison(){
    for(int n=0;n<2;n++) {
        int g = prim[n]->getGhosts();
        for(int i=-g; i<g+prim[n]->getSize(); i++){
            (*prim[n])[i].rho = NAN;
            (*prim[n])[i].p = NAN;
            (*prim[n])[i].v = {NAN,NAN,NAN};
        }
    }
}
void GridBuffers2D::poison(){
    for(int n=0;n<4;n++) {
        int g = prim[n]->getGhosts();
        for(int i=-g; i<g+prim[n]->getSizeX(); i++){
            for(int j=-g; j<g+prim[n]->getSizeY(); j++){
                (*prim[n])[i,j].rho = NAN;
                (*prim[n])[i,j].p = NAN;
                (*prim[n])[i,j].v = {NAN,NAN,NAN};
            }
        }
    }
    for(int n=0;n<2;n++) {
        int g = flux[n]->getGhosts();
        for(int i=-g; i<g+flux[n]->getSizeX(); i++){
            for(int j=-g; j<g+flux[n]->getSizeY(); j++){
                (*flux[n])[i,j].rho = NAN;
                (*flux[n])[i,j].E = NAN;
                (*flux[n])[i,j].p = {NAN,NAN,NAN};
            }
        }
    }
    #ifdef MHD
    for(int n=0;n<2;n++) {
        int g = mag[n]->getGhosts();
        for(int i=-g; i<g+mag[n]->getSizeX(); i++){
            for(int j=-g; j<g+mag[n]->getSizeY(); j++){
                (*mag[n])[i,j] = {NAN,NAN,NAN};
            }
        }
    }
    #endif
}
void GridBuffers3D::poison(){
    for(int n=0;n<6;n++) {
        int g = prim[n]->getGhosts();
        for(int i=-g; i<g+prim[n]->getSizeX(); i++){
            for(int j=-g; j<g+prim[n]->getSizeY(); j++){
                for(int k=-g; k<g+prim[n]->getSizeZ(); k++){
                    (*prim[n])[i,j,k].rho = NAN;
                    (*prim[n])[i,j,k].p = NAN;
                    (*prim[n])[i,j,k].v = {NAN,NAN,NAN};
                }
            }
        }
    }
    for(int n=0;n<3;n++) {
        int g = flux[n]->getGhosts();
        for(int i=-g; i<g+flux[n]->getSizeX(); i++){
            for(int j=-g; j<g+flux[n]->getSizeY(); j++){
                for(int k=-g; k<g+flux[n]->getSizeZ(); k++){
                    (*flux[n])[i,j,k].rho = NAN;
                    (*flux[n])[i,j,k].E = NAN;
                    (*flux[n])[i,j,k].p = {NAN,NAN,NAN};
                }
            }
        }
    }
    #ifdef CTU
    for(int n=0;n<2;n++) {
        int g = ctuprim[n]->getGhosts();
        for(int i=-g; i<g+ctuprim[n]->getSizeX(); i++){
            for(int j=-g; j<g+ctuprim[n]->getSizeY(); j++){
                for(int k=-g; k<g+ctuprim[n]->getSizeZ(); k++){
                    (*ctuprim[n])[i,j,k].rho = NAN;
                    (*ctuprim[n])[i,j,k].p = NAN;
                    (*ctuprim[n])[i,j,k].v = {NAN,NAN,NAN};
                }
            }
        }
    }
    for(int n=0;n<6;n++) {
        int g = ctuflux[n]->getGhosts();
        for(int i=-g; i<g+ctuflux[n]->getSizeX(); i++){
            for(int j=-g; j<g+ctuflux[n]->getSizeY(); j++){
                for(int k=-g; k<g+ctuflux[n]->getSizeZ(); k++){
                    (*ctuflux[n])[i,j,k].rho = NAN;
                    (*ctuflux[n])[i,j,k].E = NAN;
                    (*ctuflux[n])[i,j,k].p = {NAN,NAN,NAN};
                }
            }
        }
    }
    #endif
    #ifdef MHD
    for(int n=0;n<2;n++) {
        int g = mag[n]->getGhosts();
        for(int i=-g; i<g+mag[n]->getSizeX(); i++){
            for(int j=-g; j<g+mag[n]->getSizeY(); j++){
                for(int k=-g; k<g+mag[n]->getSizeZ(); k++){
                    (*mag[n])[i,j,k] = {NAN,NAN,NAN};
                }
            }
        }
    }
    #endif
}
