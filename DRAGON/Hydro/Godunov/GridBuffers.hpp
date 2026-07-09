//
//  Buffers.hpp
//  DRAGON
//
//  Created by Bobbie Markwick on 09/07/2026.
//

#ifndef BUFFERS_HPP
#define BUFFERS_HPP

#include "ExtendedArray.hpp"
#include "Config.h"
#pragma #once

struct PrimitiveState;
struct ConservativeState;
struct vec3;

struct GridBuffers2D{
    ExtendedArray2D<PrimitiveState>* prim[4];
    ExtendedArray2D<ConservativeState>* flux[2];
    #ifdef MHD
    ExtendedArray2D<vec3>* mag[2];
    #endif
    
    GridBuffers2D(int nx, int ny, int g){
        for(int i=0;i<4;i++)  prim[i] = new ExtendedArray2D<PrimitiveState>(nx,ny,g);
        for(int i=0;i<2;i++)  flux[i] = new ExtendedArray2D<ConservativeState>(nx,ny,g);
        #ifdef MHD
        for(int i=0;i<2;i++)  mag[i] = new ExtendedArray2D<vec3>(nx+1,ny+1,g);
        #endif
    }
    ~GridBuffers2D(){
        for(int i=0;i<4;i++)  delete prim[i];
        for(int i=0;i<2;i++)  delete flux[i];
        #ifdef MHD
        for(int i=0;i<2;i++)  delete mag[i];
        #endif
    }
    #ifdef TESTMODE
    void poison();
    #endif
    
};

struct GridBuffers3D{
    ExtendedArray3D<PrimitiveState>* prim[6];
    ExtendedArray3D<ConservativeState>* flux[3];
    #ifdef CTU
    ExtendedArray3D<PrimitiveState>* ctuprim[2];
    ExtendedArray3D<ConservativeState>* ctuflux[6];
    #endif
    #ifdef MHD
    ExtendedArray3D<vec3>* mag[2];
    #endif
    
    GridBuffers3D(int nx, int ny, int nz, int g){
        for(int i=0;i<6;i++)  prim[i] = new ExtendedArray3D<PrimitiveState>(nx,ny,nz,g);
        for(int i=0;i<3;i++)  flux[i] = new ExtendedArray3D<ConservativeState>(nx,ny,nz,g);
        #ifdef CTU
        for(int i=0;i<2;i++)  ctuprim[i] = new ExtendedArray3D<PrimitiveState>(nx,ny,nz,g);
        for(int i=0;i<6;i++)  ctuflux[i] = new ExtendedArray3D<ConservativeState>(nx,ny,nz,g);
        #endif
        #ifdef MHD
        for(int i=0;i<2;i++)  mag[i] = new ExtendedArray3D<vec3>(nx+1,ny+1,nz+1,g);
        #endif
    }
    ~GridBuffers3D(){
        for(int i=0;i<6;i++)  delete prim[i];
        for(int i=0;i<3;i++)  delete flux[i];
        #ifdef CTU
        for(int i=0;i<2;i++)  delete ctuprim[i];
        for(int i=0;i<6;i++)  delete ctuflux[i];
        #endif
        #ifdef MHD
        for(int i=0;i<2;i++)  delete mag[i];
        #endif
    }
    #ifdef TESTMODE
    void poison();
    #endif    
};


#endif
