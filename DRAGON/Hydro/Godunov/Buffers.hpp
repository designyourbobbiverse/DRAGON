//
//  Buffers.hpp
//  DRAGON
//
//  Created by Bobbie Markwick on 09/07/2026.
//

#ifndef BUFFERS_HPP
#define BUFFERS_HPP

#include "ExtendedArray.hpp"
#include "FluidElement.hpp"
#include <cmath>
#pragma #once


struct Buffers2D{
    ExtendedArray2D<PrimitiveState>* prim[4];
    ExtendedArray2D<ConservativeState>* flux[2];
    #ifdef MHD
    ExtendedArray2D<vec3>* mag[2];
    #endif
    
    Buffers2D(int nx, int ny, int g){
        for(int i=0;i<4;i++) {
            prim[i] = new ExtendedArray2D<PrimitiveState>(nx,ny,g);
        }
        for(int i=0;i<2;i++) {
            flux[i] = new ExtendedArray2D<ConservativeState>(nx,ny,g);
        }
        #ifdef MHD
        for(int i=0;i<2;i++) {
            mag[i] = new ExtendedArray2D<vec3>(nx+1,ny+1,g);
        }
        #endif
    }
    
    void poison(){
        for(int n=0;n<4;n++) {
            int g = prim[n]->getGhosts();
            for(int i=-g; i<g+prim[n]->getSizeX(); i++){
                for(int j=0; j<g+prim[n]->getSizeY(); j++){
                    (*prim[n])[i,j].rho = NAN;
                    (*prim[n])[i,j].p = NAN;
                    (*prim[n])[i,j].v = {NAN,NAN,NAN};
                }
            }
        }
        for(int n=0;n<2;n++) {
            int g = flux[n]->getGhosts();
            for(int i=-g; i<g+flux[n]->getSizeX(); i++){
                for(int j=0; j<g+flux[n]->getSizeY(); j++){
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
                for(int j=0; j<g+mag[n]->getSizeY(); j++){
                    (*mag[n])[i,j] = {NAN,NAN,NAN};
                }
            }
        }
        #endif
    }
    
};

struct Buffers3D{
    ExtendedArray3D<PrimitiveState>* prim[6];
    ExtendedArray3D<ConservativeState>* flux[3];
    #ifdef CTU
    ExtendedArray3D<PrimitiveState>* ctuprim[2];
    ExtendedArray3D<ConservativeState>* ctuflux[6];
    #endif
    #ifdef MHD
    ExtendedArray3D<vec3>* mag[2];
    #endif
    
    Buffers3D(int nx, int ny, int nz, int g){
        for(int i=0;i<6;i++) {
            prim[i] = new ExtendedArray3D<PrimitiveState>(nx,ny,nz,g);
        }
        for(int i=0;i<3;i++) {
            flux[i] = new ExtendedArray3D<ConservativeState>(nx,ny,nz,g);
        }
        #ifdef CTU
        for(int i=0;i<2;i++) {
            ctuprim[i] = new ExtendedArray3D<PrimitiveState>(nx,ny,nz,g);
        }
        for(int i=0;i<6;i++) {
            ctuflux[i] = new ExtendedArray3D<ConservativeState>(nx,ny,nz,g);
        }
        #endif
        #ifdef MHD
        for(int i=0;i<2;i++) {
            mag[i] = new ExtendedArray3D<vec3>(nx+1,ny+1,nz+1,g);
        }
        #endif
    }
    
    void poison(){
        for(int n=0;n<6;n++) {
            int g = prim[n]->getGhosts();
            for(int i=-g; i<g+prim[n]->getSizeX(); i++){
                for(int j=0; j<g+prim[n]->getSizeY(); j++){
                    for(int k=0; k<g+prim[n]->getSizeZ(); k++){
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
                for(int j=0; j<g+flux[n]->getSizeY(); j++){
                    for(int k=0; k<g+flux[n]->getSizeZ(); k++){
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
                for(int j=0; j<g+ctuprim[n]->getSizeY(); j++){
                    for(int k=0; k<g+ctuprim[n]->getSizeZ(); k++){
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
                for(int j=0; j<g+ctuflux[n]->getSizeY(); j++){
                    for(int k=0; k<g+ctuflux[n]->getSizeZ(); k++){
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
                for(int j=0; j<g+mag[n]->getSizeY(); j++){
                    for(int k=0; k<g+mag[n]->getSizeZ(); k++){
                        (*mag[n])[i,j,k] = {NAN,NAN,NAN};
                    }
                }
            }
        }
        #endif
    }
    
};


#endif
