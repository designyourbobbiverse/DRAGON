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
bool requestRestart(std::string msg = ""); //Request that all threads restart the timestep. Returns true iff in multithread mode
//If any thread requests (or has already requested) a step-restart, all waitFor functions will promptly return false. If on the other hand the waitFor condition is successfully reached, the waitFor functions will return true.
bool waitForRelease(); //Call this before entering memory-heavy phase 1 to avoid too many threads fighting for memory usage.
    //The exact number of threads permitted to be in phase 1 at a time is controlled in DRAGONWING_Config.hpp
void reportCheckpoint1(); //Report that this thread's subgrid has completed phase 1 by successfully computing (but not committing) its update
bool waitForCheckpoint1(); //Wait until all threads have reported checkpoint 1
void reportCheckpoint2(); //Report that this thread's subgrid has fully completed (and committed) its update
//To wait for checkpoint 2, call ThreadPool::waitForCompletion() from a thread outside the ThreadPool.

//MARK: Memory Management
//Parameters: N is the number of arrays to return
    //nx/ny/nz/g are the size parameters for the arrays themselves
    //Do not assume you will get an array of zeros: Arrays are recycled and not cleared between uses, so they may contain stale data.
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
