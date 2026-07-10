//
//  DragonWing.hpp
//  DRAGONWING
//
//  Created by Bobbie Markwick on 11/06/2026.
//

#ifndef DRAGON_WING_hpp
#define DRAGON_WING_hpp
#include <string>

class Grid;
struct PrimitiveState;
struct ConservativeState;
struct vec3;
template <class T> struct ExtendedArray1D;
template <class T> struct ExtendedArray2D;
template <class T> struct ExtendedArray3D;

namespace DRAGONWING{

//MARK: Multithreading
void* launchParallel(Grid* grid, double dt);

void initialize(int nthreads);
bool requestRestart(std::string msg = ""); //Returns true if genuinely in multithread mode
std::string restartMsg();
void reportCheckpoint1();
void reportCheckpoint2();
bool waitForRelease();
bool waitForCheckpoint1(); //Returns false iff someone requested a restart
bool waitForCheckpoint2(); //Returns false iff someone requested a restart

//MARK: Memory Management
ExtendedArray1D<PrimitiveState>* requestPrimitiveArray(int nx, int g);
ExtendedArray2D<PrimitiveState>* requestPrimitiveArray(int nx, int ny, int g);
ExtendedArray3D<PrimitiveState>* requestPrimitiveArray(int nx, int ny, int nz, int g);

ExtendedArray2D<ConservativeState>* requestFluxArray(int nx, int ny, int g);
ExtendedArray3D<ConservativeState>* requestFluxArray(int nx, int ny, int nz, int g);

ExtendedArray2D<vec3>* requestVec3Array(int nx, int ny, int g);
ExtendedArray3D<vec3>* requestVec3Array(int nx, int ny, int nz, int g);

void releaseArray(ExtendedArray1D<PrimitiveState>* arr);
void releaseArray(ExtendedArray2D<PrimitiveState>* arr);
void releaseArray(ExtendedArray3D<PrimitiveState>* arr);

void releaseArray(ExtendedArray1D<ConservativeState>* arr);
void releaseArray(ExtendedArray2D<ConservativeState>* arr);
void releaseArray(ExtendedArray3D<ConservativeState>* arr);

void releaseArray(ExtendedArray2D<vec3>* arr);
void releaseArray(ExtendedArray3D<vec3>* arr);

void purgeAllBuffers();

}


#endif
