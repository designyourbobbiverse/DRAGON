//
//  DragonWing.hpp
//  DRAGONWING
//
//  Created by Bobbie Markwick on 11/06/2026.
//

#ifndef DRAGON_WING_hpp
#define DRAGON_WING_hpp
#include <string>
#include <vector>
#include "ArrayGuard.hpp"
class Grid;


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
//Return Type: ArrayGuard
//Automatically releases the arrays when it goes out of scope. Can release early by calling guard.release()
//guard[i] or guard.get(i) returns a pointer to the ith array
ArrayGuard<ExtendedArray1D<PrimitiveState>> requestPrimitiveArrays(int N, int nx, int g);
ArrayGuard<ExtendedArray1D<ConservativeState>> requestFluxArrays(int N, int nx, int g);
ArrayGuard<ExtendedArray1D<vec3>> requestVec3Arrays(int N, int nx, int g);
ArrayGuard<ExtendedArray2D<PrimitiveState>> requestPrimitiveArrays(int N, int nx, int ny, int g);
ArrayGuard<ExtendedArray2D<ConservativeState>> requestFluxArrays(int N, int nx, int ny, int g);
ArrayGuard<ExtendedArray2D<vec3>> requestVec3Arrays(int N, int nx, int ny, int g);
ArrayGuard<ExtendedArray3D<PrimitiveState>> requestPrimitiveArrays(int N, int nx, int ny, int nz, int g);
ArrayGuard<ExtendedArray3D<ConservativeState>> requestFluxArrays(int N, int nx, int ny, int nz, int g);
ArrayGuard<ExtendedArray3D<vec3>> requestVec3Arrays(int N, int nx, int ny, int nz, int g);





void purgeAllBuffers();

}


#endif
