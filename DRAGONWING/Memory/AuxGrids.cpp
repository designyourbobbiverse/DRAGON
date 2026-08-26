//
//  BufferGrids.cpp
//  DRAGONWING
//
//  Created by Bobbie Markwick on 9/07/2026.
//

#include "DragonWing.hpp"

#include "DRAGONWING_Config.hpp"
using namespace DRAGON;

#if defined(REUSE_AUX_GRIDS) && !defined(MULTITHREAD_UNAVAILABLE)
//Keep track of all the arrays that have been allocated, and which ones are in use
#include <vector>
template <class T> struct ArrayItem{ T* array; bool active; };
static std::vector<ArrayItem<ExtendedArray1D<PrimitiveState>>> prim1D;
static std::vector<ArrayItem<ExtendedArray2D<PrimitiveState>>> prim2D;
static std::vector<ArrayItem<ExtendedArray3D<PrimitiveState>>> prim3D;
static std::vector<ArrayItem<ExtendedArray1D<ConservativeState>>> flux1D;
static std::vector<ArrayItem<ExtendedArray2D<ConservativeState>>> flux2D;
static std::vector<ArrayItem<ExtendedArray3D<ConservativeState>>> flux3D;
static std::vector<ArrayItem<ExtendedArray1D<vec3>>> vec1D;
static std::vector<ArrayItem<ExtendedArray2D<vec3>>> vec2D;
static std::vector<ArrayItem<ExtendedArray3D<vec3>>> vec3D;

//Make sure access is threadsafe
#include <mutex>
static std::mutex pm;
static std::mutex cm;
static std::mutex vm;

//MARK: Request 1D
template <class T>
DRAGONWING::ArrayGuard<ExtendedArray1D<T>> requestArrays(std::vector<ArrayItem<ExtendedArray1D<T>>>& prealloc, std::mutex& m,  int N, int nx, int g){
    auto selected = std::vector<ExtendedArray1D<T>*>();

    //Search for an existing available item
    std::unique_lock lock(m);
    for (auto& item : prealloc) {
        if (item.active) continue; //Check that the array isn't currently in use
        auto arr = item.array;
        if (arr->getSize() == nx && arr->getGhosts() == g) { //Check that the array is of the correct size
            item.active = true;
            selected.push_back(arr);
            if (selected.size() == N) return DRAGONWING::ArrayGuard(selected); //If we have enough arrays, we are done
        }
    }
    lock.unlock(); //The next part doesn't require synchronization, so don't hold everyone else up
    
    //Create a new grids if needed
    auto new_arrs = std::vector<ExtendedArray1D<T>*>();
    while(selected.size() < N) {
        auto arr = new ExtendedArray1D<T>(nx,g);
        new_arrs.push_back(arr);
        selected.push_back(arr);
    }
    //Add the new arrays to the global list
    lock.lock(); //accessing prealloc means we need to relock the mutex
    for (auto& a : new_arrs) prealloc.push_back({a, true});
    lock.unlock();
    
    return DRAGONWING::ArrayGuard(selected);
}
//Call the above function with the appropriate buffer collection and mutex
DRAGONWING::ArrayGuard<ExtendedArray1D<PrimitiveState>> DRAGONWING::requestPrimitiveArrays(int N, int nx, int g){
    return requestArrays(prim1D, pm, N, nx, g);
}
DRAGONWING::ArrayGuard<ExtendedArray1D<ConservativeState>> DRAGONWING::requestFluxArrays(int N, int nx, int g){
    return requestArrays(flux1D, cm, N, nx, g);
}
DRAGONWING::ArrayGuard<ExtendedArray1D<vec3>> DRAGONWING::requestVec3Arrays(int N, int nx, int g){
    return requestArrays(vec1D, vm, N, nx, g);
}

//MARK: Request 2D
template <class T>
DRAGONWING::ArrayGuard<ExtendedArray2D<T>> requestArrays(std::vector<ArrayItem<ExtendedArray2D<T>>>& prealloc, std::mutex& m,  int N, int nx, int ny, int g){
    auto selected = std::vector<ExtendedArray2D<T>*>();

    //Search for an existing available item
    std::unique_lock lock(m);
    for (auto& item : prealloc) {
        if (item.active) continue; //Check that the array isn't currently in use
        auto arr = item.array;
        if (arr->getSizeX() == nx && arr->getSizeY() == ny && arr->getGhosts() == g) { //Check that the array is of the correct size
            item.active = true;
            selected.push_back(arr);
            if (selected.size() == N) return DRAGONWING::ArrayGuard(selected); //If we have enough arrays, we are done
        }
    }
    lock.unlock(); //The next part doesn't require synchronization, so don't hold everyone else up

    //Create a new grids if needed
    auto new_arrs = std::vector<ExtendedArray2D<T>*>();
    while(selected.size() < N) {
        auto arr = new ExtendedArray2D<T>(nx,ny,g);
        new_arrs.push_back(arr);
        selected.push_back(arr);
    }
    //Add the new arrays to the global list
    lock.lock(); //accessing prealloc means we need to relock the mutex
    for (auto& a : new_arrs) prealloc.push_back({a, true});
    lock.unlock();
    
    return DRAGONWING::ArrayGuard(selected);
}
//Call the above function with the appropriate buffer collection and mutex
DRAGONWING::ArrayGuard<ExtendedArray2D<PrimitiveState>> DRAGONWING::requestPrimitiveArrays(int N, int nx, int ny, int g){
    return requestArrays(prim2D, pm, N, nx, ny, g);
}

DRAGONWING::ArrayGuard<ExtendedArray2D<ConservativeState>> DRAGONWING::requestFluxArrays(int N, int nx, int ny, int g){
    return requestArrays(flux2D, cm, N, nx, ny, g);
}
DRAGONWING::ArrayGuard<ExtendedArray2D<vec3>> DRAGONWING::requestVec3Arrays(int N, int nx, int ny, int g){
    return requestArrays(vec2D, vm, N, nx, ny, g);
}

//MARK: Request 3D
template <class T>
DRAGONWING::ArrayGuard<ExtendedArray3D<T>> requestArrays(std::vector<ArrayItem<ExtendedArray3D<T>>>& prealloc, std::mutex& m,  int N, int nx, int ny, int nz, int g){
    auto selected = std::vector<ExtendedArray3D<T>*>();

    //Search for an existing available item
    std::unique_lock lock(m);
    for (auto& item : prealloc) {
        if (item.active) continue; //Check that the array isn't currently in use
        auto arr = item.array; //Check that the array is of the correct size
        if (arr->getSizeX() == nx && arr->getSizeY() == ny && arr->getSizeZ() == nz && arr->getGhosts() == g) {
            item.active = true;
            selected.push_back(arr);
            if (selected.size() == N) return DRAGONWING::ArrayGuard(selected); //If we have enough arrays, we are done
        }
    }
    lock.unlock(); //The next part doesn't require synchronization, so don't hold everyone else up
    
    //Create a new grids if needed
    auto new_arrs = std::vector<ExtendedArray3D<T>*>();
    while(selected.size() < N) {
        auto arr = new ExtendedArray3D<T>(nx,ny,nz,g);
        new_arrs.push_back(arr);
        selected.push_back(arr);
    }
    //Add the new arrays to the global list
    lock.lock(); //accessing prealloc means we need to relock the mutex
    for (auto& a : new_arrs) prealloc.push_back({a, true});
    lock.unlock();
    
    return DRAGONWING::ArrayGuard(selected);
}

//Call the above function with the appropriate buffer collection and mutex
DRAGONWING::ArrayGuard<ExtendedArray3D<PrimitiveState>> DRAGONWING::requestPrimitiveArrays(int N, int nx, int ny, int nz, int g){
    return requestArrays(prim3D, pm, N, nx, ny, nz, g);
}
DRAGONWING::ArrayGuard<ExtendedArray3D<ConservativeState>> DRAGONWING::requestFluxArrays(int N, int nx, int ny, int nz, int g){
    return requestArrays(flux3D, cm, N, nx, ny, nz, g);
}
DRAGONWING::ArrayGuard<ExtendedArray3D<vec3>> DRAGONWING::requestVec3Arrays(int N, int nx, int ny, int nz, int g){
    return requestArrays(vec3D, vm, N, nx, ny, nz, g);
}

//MARK: Item Release
template <class T>
void releaseArray(std::vector<ArrayItem<T>>& arrs, T* arr){
    //Search for an existing available item
    for (auto& item : arrs) {
        if (item.array != arr) continue;
        item.active = false; //If we found this array, mark it as available for use
        return;
    }
    //Not found in the vector (perhaps purgeAllBuffers() was called), delete to prevent a leak
    delete arr;
}

//Lock the appropriate mutex, then release the array (either back to the pool, or from memory if no longer assoicated with the pool)
void DRAGONWING::Internal::releaseArray(ExtendedArray1D<PrimitiveState>* arr){ std::lock_guard lock(pm); releaseArray(prim1D, arr); }
void DRAGONWING::Internal::releaseArray(ExtendedArray2D<PrimitiveState>* arr){ std::lock_guard lock(pm); releaseArray(prim2D, arr); }
void DRAGONWING::Internal::releaseArray(ExtendedArray3D<PrimitiveState>* arr){ std::lock_guard lock(pm); releaseArray(prim3D, arr); }

void DRAGONWING::Internal::releaseArray(ExtendedArray1D<ConservativeState>* arr){ std::lock_guard lock(cm); releaseArray(flux1D, arr); }
void DRAGONWING::Internal::releaseArray(ExtendedArray2D<ConservativeState>* arr){ std::lock_guard lock(cm); releaseArray(flux2D, arr); }
void DRAGONWING::Internal::releaseArray(ExtendedArray3D<ConservativeState>* arr){ std::lock_guard lock(cm); releaseArray(flux3D, arr); }

void DRAGONWING::Internal::releaseArray(ExtendedArray1D<vec3>* arr){ std::lock_guard lock(vm); releaseArray(vec1D, arr); }
void DRAGONWING::Internal::releaseArray(ExtendedArray2D<vec3>* arr){ std::lock_guard lock(vm); releaseArray(vec2D, arr); }
void DRAGONWING::Internal::releaseArray(ExtendedArray3D<vec3>* arr){ std::lock_guard lock(vm); releaseArray(vec3D, arr); }

//MARK: Cleanup
template <class T>
static void purgeBuffers(std::vector<ArrayItem<T>>& array){
    for (auto& item : array) {
        //If the item is inactive, we need to deallocate it from system memory
        //If the item is active, it will get deleted when releaseArray() is called
        if (!item.active)  delete item.array;
    }
}
void DRAGONWING::purgeAllBuffers(){
    { //Purge all of the buffers for each type, locking the appropraite mutex while doing so
        std::lock_guard lock(pm);
        purgeBuffers(prim1D);
        purgeBuffers(prim2D);
        purgeBuffers(prim3D);
        prim1D.clear();
        prim2D.clear();
        prim3D.clear();
    }
    {
        std::lock_guard lock(cm);
        purgeBuffers(flux1D);
        purgeBuffers(flux2D);
        purgeBuffers(flux3D);
        flux1D.clear();
        flux2D.clear();
        flux3D.clear();
    }
    {
        std::lock_guard lock(vm);
        purgeBuffers(vec1D);
        purgeBuffers(vec2D);
        purgeBuffers(vec3D);
        vec1D.clear();
        vec2D.clear();
        vec3D.clear();
    }
}
#else
//MARK: No-Reuse (wrap new/delete)
//Allocate N arrays
template <class T>
static DRAGONWING::ArrayGuard<ExtendedArray1D<T>> requestArrays(int N, int nx, int g){
    auto selected = std::vector<ExtendedArray1D<T>*>();
    selected.reserve(N);
    while(selected.size() < N) selected.push_back(new ExtendedArray1D<T>(nx,g));
    return DRAGONWING::ArrayGuard(selected);
}
DRAGONWING::ArrayGuard<ExtendedArray1D<PrimitiveState>> DRAGONWING::requestPrimitiveArrays(int N, int nx, int g){
    return requestArrays<PrimitiveState>(N, nx, g);
}
DRAGONWING::ArrayGuard<ExtendedArray1D<ConservativeState>> DRAGONWING::requestFluxArrays(int N, int nx, int g){
    return requestArrays<ConservativeState>(N, nx, g);
}
DRAGONWING::ArrayGuard<ExtendedArray1D<vec3>> DRAGONWING::requestVec3Arrays(int N, int nx, int g){
    return requestArrays<vec3>(N, nx, g);
}

template <class T>
static DRAGONWING::ArrayGuard<ExtendedArray2D<T>> requestArrays(int N, int nx, int ny, int g){
    auto selected = std::vector<ExtendedArray2D<T>*>();
    selected.reserve(N);
    while(selected.size() < N) selected.push_back(new ExtendedArray2D<T>(nx,ny,g));
    return DRAGONWING::ArrayGuard(selected);
}
DRAGONWING::ArrayGuard<ExtendedArray2D<PrimitiveState>> DRAGONWING::requestPrimitiveArrays(int N, int nx, int ny, int g){
    return requestArrays<PrimitiveState>(N, nx, ny, g);
}

DRAGONWING::ArrayGuard<ExtendedArray2D<ConservativeState>> DRAGONWING::requestFluxArrays(int N, int nx, int ny, int g){
    return requestArrays<ConservativeState>(N, nx, ny, g);
}
DRAGONWING::ArrayGuard<ExtendedArray2D<vec3>> DRAGONWING::requestVec3Arrays(int N, int nx, int ny, int g){
    return requestArrays<vec3>(N, nx, ny, g);
}


template <class T>
DRAGONWING::ArrayGuard<ExtendedArray3D<T>> requestArrays(int N, int nx, int ny, int nz, int g){
    auto selected = std::vector<ExtendedArray3D<T>*>();
    selected.reserve(N);
    while(selected.size() < N) selected.push_back(new ExtendedArray3D<T>(nx,ny,nz,g));
    return DRAGONWING::ArrayGuard(selected);
}
DRAGONWING::ArrayGuard<ExtendedArray3D<PrimitiveState>> DRAGONWING::requestPrimitiveArrays(int N, int nx, int ny, int nz, int g){
    return requestArrays<PrimitiveState>(N, nx, ny, nz, g);
}

DRAGONWING::ArrayGuard<ExtendedArray3D<ConservativeState>> DRAGONWING::requestFluxArrays(int N, int nx, int ny, int nz, int g){
    return requestArrays<ConservativeState>(N, nx, ny, nz, g);
}
DRAGONWING::ArrayGuard<ExtendedArray3D<vec3>> DRAGONWING::requestVec3Arrays(int N, int nx, int ny, int nz, int g){
    return requestArrays<vec3>(N, nx, ny, nz, g);
}

//Without a pool, release just means delete
void DRAGONWING::Internal::releaseArray(ExtendedArray1D<PrimitiveState>* arr){ delete arr; }
void DRAGONWING::Internal::releaseArray(ExtendedArray2D<PrimitiveState>* arr){ delete arr; }
void DRAGONWING::Internal::releaseArray(ExtendedArray3D<PrimitiveState>* arr){ delete arr; }

void DRAGONWING::Internal::releaseArray(ExtendedArray1D<ConservativeState>* arr){ delete arr; }
void DRAGONWING::Internal::releaseArray(ExtendedArray2D<ConservativeState>* arr){ delete arr; }
void DRAGONWING::Internal::releaseArray(ExtendedArray3D<ConservativeState>* arr){ delete arr; }

void DRAGONWING::Internal::releaseArray(ExtendedArray1D<vec3>* arr){ delete arr; }
void DRAGONWING::Internal::releaseArray(ExtendedArray2D<vec3>* arr){ delete arr; }
void DRAGONWING::Internal::releaseArray(ExtendedArray3D<vec3>* arr){ delete arr; }

void DRAGONWING::purgeAllBuffers(){ } //No pool to purge
#endif
