//
//  DragonWing.hpp
//  DRAGONWING
//
//  Created by Bobbie Markwick on 11/06/2026.
//

#ifndef DRAGON_WING_hpp
#define DRAGON_WING_hpp

#include "ArrayGuard.hpp"
#include "ThreadPool.hpp"

namespace DRAGONWING{

//MARK: Multithreading
//Multithreading is done via the ThreadPool class. Important functions:
    //ThreadPool(int n); //Creates a pool with n threads
    //void* launchParallel(Grid* grid, double dt); //Creates a thread which calls grid->advance(dt, false);
    //std::string restartMsg(); //returns (+clears) the error message if something requested a restart.
    //bool waitForCompletion(); //Waits for all threads to finish. Returns false iff anyone requested a restart

//Everything below automatically invokes the ThreadPool associated with the current thread. No need to pass the object around yourself.
bool requestRestart(std::string msg = ""); //Returns true if genuinely in multithread mode
void reportCheckpoint1();
void reportCheckpoint2();
bool waitForRelease(); //Returns false iff someone requested a restart
bool waitForCheckpoint1(); //Returns false iff someone requested a restart
bool waitForCheckpoint2(); //Returns false iff someone requested a restart

//MARK: Memory Management
//Parameters: N is the number of arrays to return
    //nx/ny/nz/g are the size parameters for the arrays themselves
//Return Type: ArrayGuard
    //Automatically releases the arrays when it goes out of scope. Can release early by calling guard.release()
    //guard[i] or guard.get(i) returns a pointer to the ith array
ArrayGuard<DRAGON::ExtendedArray1D<DRAGON::PrimitiveState>> requestPrimitiveArrays(int N, int nx, int g);
ArrayGuard<DRAGON::ExtendedArray1D<DRAGON::ConservativeState>> requestFluxArrays(int N, int nx, int g);
ArrayGuard<DRAGON::ExtendedArray1D<DRAGON::vec3>> requestVec3Arrays(int N, int nx, int g);
ArrayGuard<DRAGON::ExtendedArray2D<DRAGON::PrimitiveState>> requestPrimitiveArrays(int N, int nx, int ny, int g);
ArrayGuard<DRAGON::ExtendedArray2D<DRAGON::ConservativeState>> requestFluxArrays(int N, int nx, int ny, int g);
ArrayGuard<DRAGON::ExtendedArray2D<DRAGON::vec3>> requestVec3Arrays(int N, int nx, int ny, int g);
ArrayGuard<DRAGON::ExtendedArray3D<DRAGON::PrimitiveState>> requestPrimitiveArrays(int N, int nx, int ny, int nz, int g);
ArrayGuard<DRAGON::ExtendedArray3D<DRAGON::ConservativeState>> requestFluxArrays(int N, int nx, int ny, int nz, int g);
ArrayGuard<DRAGON::ExtendedArray3D<DRAGON::vec3>> requestVec3Arrays(int N, int nx, int ny, int nz, int g);

void purgeAllBuffers();

}


#endif
