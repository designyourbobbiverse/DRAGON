//
//  BufferGrids.cpp
//  DRAGONWING
//
//  Created by Bobbie Markwick on 9/07/2026.
//

#include "DragonWing.hpp"
#include "DRAGONWING_Config.hpp"
#include "ExtendedArray.hpp"
#include "FluidElement.hpp"

#ifdef REUSE_AUX_GRIDS
#include <vector>

#ifndef MULTITHREAD_UNAVAILABLE
#include <mutex>
static std::mutex mutex;
#define MUTEX_LOCK_SCOPE std::lock_guard mutex_lock_guard(mutex);
#else
#define MUTEX_LOCK_SCOPE 
#endif

namespace {
template <class T> struct ArrayItem{ T* array; bool active; };

std::vector<ArrayItem<ExtendedArray1D<PrimitiveState>>> prim1D;
std::vector<ArrayItem<ExtendedArray2D<PrimitiveState>>> prim2D;
std::vector<ArrayItem<ExtendedArray3D<PrimitiveState>>> prim3D;

std::vector<ArrayItem<ExtendedArray2D<ConservativeState>>> flux2D;
std::vector<ArrayItem<ExtendedArray3D<ConservativeState>>> flux3D;

std::vector<ArrayItem<ExtendedArray2D<vec3>>> vec2D;
std::vector<ArrayItem<ExtendedArray3D<vec3>>> vec3D;
}

//MARK: Request 1D
ExtendedArray1D<PrimitiveState>* DRAGONWING::requestPrimitiveArray(int nx, int g){
    MUTEX_LOCK_SCOPE
    //Search for an existing available item
    for(auto& item : prim1D){
        if(item.active) continue; //Someone else is using this one
        auto arr = item.array;
        if(arr->getSize() == nx && arr->getGhosts() == g){
            item.active = true;
            return arr;
        }
    }
    //Create a new grid if needed
    auto arr = new ExtendedArray1D<PrimitiveState>(nx,g);
    prim1D.push_back({arr, true});
    return arr;
}

//MARK: Request 2D
template <class T>
ExtendedArray2D<T>* request2DArray(std::vector<ArrayItem<ExtendedArray2D<T>>> arrs, int nx, int ny, int g){
    MUTEX_LOCK_SCOPE
    //Search for an existing available item
    for(auto& item : arrs){
        if(item.active) continue; //Someone else is using this one
        auto arr = item.array;
        if(arr->getSizeX() == nx && arr->getSizeY() == ny && arr->getGhosts() == g){
            item.active = true;
            return arr;
        }
    }
    //Create a new grid if needed
    auto arr = new ExtendedArray2D<T>(nx,ny,g);
    arrs.push_back({arr, true});
    return arr;
}

ExtendedArray2D<PrimitiveState>* DRAGONWING::requestPrimitiveArray(int nx, int ny, int g){
    return request2DArray(prim2D, nx,ny,g);
}
ExtendedArray2D<ConservativeState>* DRAGONWING::requestFluxArray(int nx, int ny, int g){
    return request2DArray(flux2D, nx,ny,g);
}
ExtendedArray2D<vec3>* DRAGONWING::requestVec3Array(int nx, int ny, int g){
    return request2DArray(vec2D, nx,ny,g);
}
//MARK: Request 3D
template <class T>
ExtendedArray3D<T>* request3DArray(std::vector<ArrayItem<ExtendedArray3D<T>>> arrs, int nx, int ny, int nz, int g){
    MUTEX_LOCK_SCOPE
    //Search for an existing available item
    for(auto& item : arrs){
        if(item.active) continue; //Someone else is using this one
        auto arr = item.array;
        if(arr->getSizeX() == nx && arr->getSizeY() == ny && arr->getSizeZ() == nz && arr->getGhosts() == g){
            item.active = true;
            return arr;
        }
    }
    //Create a new grid if needed
    auto arr = new ExtendedArray3D<T>(nx,ny,nz,g);
    arrs.push_back({arr, true});
    return arr;
}

ExtendedArray3D<PrimitiveState>* DRAGONWING::requestPrimitiveArray(int nx, int ny, int nz, int g){
    return request3DArray(prim3D, nx,ny,nz,g);
}
ExtendedArray3D<ConservativeState>* DRAGONWING::requestFluxArray(int nx, int ny, int nz, int g){
    return request3DArray(flux3D, nx,ny,nz,g);
}
ExtendedArray3D<vec3>* DRAGONWING::requestVec3Array(int nx, int ny, int nz, int g){
    return request3DArray(vec3D, nx,ny,nz,g);
}

//MARK: Item Release
template <class T>
void releaseArray(std::vector<ArrayItem<T>> arrs, T* arr){
    MUTEX_LOCK_SCOPE
    //Search for an existing available item
    for(auto& item : arrs){
        if(item.array != arr) continue;
        item.active = false;
        return;
    }
    //Not found in the vector, delete to prevent leak
    delete arr;
}

void DRAGONWING::releaseArray(ExtendedArray1D<PrimitiveState>* arr){ releaseArray(prim1D, arr); }
void DRAGONWING::releaseArray(ExtendedArray2D<PrimitiveState>* arr){ releaseArray(prim2D, arr); }
void DRAGONWING::releaseArray(ExtendedArray3D<PrimitiveState>* arr){ releaseArray(prim3D, arr); }

void DRAGONWING::releaseArray(ExtendedArray2D<ConservativeState>* arr){ releaseArray(flux2D, arr); }
void DRAGONWING::releaseArray(ExtendedArray3D<ConservativeState>* arr){ releaseArray(flux3D, arr); }

void DRAGONWING::releaseArray(ExtendedArray2D<vec3>* arr){ releaseArray(vec2D, arr); }
void DRAGONWING::releaseArray(ExtendedArray3D<vec3>* arr){ releaseArray(vec3D, arr); }

//MARK: Cleanup
template <class T>
void purgeBuffers(std::vector<ArrayItem<T>> array){
    for(auto& item : array){
        if(!item.active)  delete item.array;
    }
}
void DRAGONWING::purgeAllBuffers(){
    MUTEX_LOCK_SCOPE
    purgeBuffers(prim1D);
    purgeBuffers(prim2D);
    purgeBuffers(prim3D);
    purgeBuffers(flux2D);
    purgeBuffers(flux3D);
    
    prim1D.clear();
    prim2D.clear();
    prim3D.clear();
    flux2D.clear();
    flux3D.clear();
    vec2D.clear();
    vec3D.clear();
}
#else
//MARK: No-Reuse (wrap new/delete)

ExtendedArray1D<PrimitiveState>* DRAGONWING::requestPrimitiveArray(int nx, int g){
    return new ExtendedArray1D<PrimitiveState>(nx,g);
}

ExtendedArray2D<PrimitiveState>* DRAGONWING::requestPrimitiveArray(int nx, int ny, int g){
    return new ExtendedArray2D<PrimitiveState>(nx,ny,g);
}
ExtendedArray2D<ConservativeState>* DRAGONWING::requestFluxArray(int nx, int ny, int g){
    return new ExtendedArray2D<ConservativeState>(nx,ny,g);
}
ExtendedArray2D<vec3>* DRAGONWING::requestVec3Array(int nx, int ny, int g){
    return new ExtendedArray2D<vec3>(nx,ny,g);
}
ExtendedArray3D<PrimitiveState>* DRAGONWING::requestPrimitiveArray(int nx, int ny, int nz, int g){
    return new ExtendedArray3D<PrimitiveState>(nx,ny,nz,g);
}
ExtendedArray3D<ConservativeState>* DRAGONWING::requestFluxArray(int nx, int ny, int nz, int g){
    return new ExtendedArray3D<ConservativeState>(nx,ny,nz,g);
}
ExtendedArray3D<vec3>* DRAGONWING::requestVec3Array(int nx, int ny, int nz, int g){
    return new ExtendedArray3D<vec3>(nx,ny,nz,g);
}

void DRAGONWING::releaseArray(ExtendedArray1D<PrimitiveState>* arr){ delete arr; }
void DRAGONWING::releaseArray(ExtendedArray2D<PrimitiveState>* arr){ delete arr; }
void DRAGONWING::releaseArray(ExtendedArray3D<PrimitiveState>* arr){ delete arr; }

void DRAGONWING::releaseArray(ExtendedArray2D<ConservativeState>* arr){ delete arr; }
void DRAGONWING::releaseArray(ExtendedArray3D<ConservativeState>* arr){ delete arr; }

void DRAGONWING::releaseArray(ExtendedArray2D<vec3>* arr){ delete arr; }
void DRAGONWING::releaseArray(ExtendedArray3D<vec3>* arr){ delete arr; }

void DRAGONWING::purgeAllBuffers(){ }
#endif
