//
//  ThreadPool.hpp
//  DRAGONWING/Multithreading
//
//  Created by Bobbie Markwick on 31/07/2026.
//

#ifndef ThreadPool_hpp
#define ThreadPool_hpp

#include <string>   //For error messages
namespace DRAGON{ class Grid;}

#include "DRAGONWING_Config.hpp"
#ifndef MULTITHREAD_UNAVAILABLE
#include <thread>             //Execute tasks in parallel
#include <mutex>              //For synchronisation
#include <condition_variable> //For synchronisation
#include <deque>              //For thread tracking
#include <atomic>             //For Fallback Counting

namespace DRAGONWING{
struct ThreadArgs{ //The things a thread needs to know to do its work
    DRAGON::Grid* grid; //The subgrid this thread is supposed to run on
    double dt; //The time the subgird is supposed to advance by
};
}
#endif

namespace DRAGONWING {

class ThreadPool{
private:
    #ifndef MULTITHREAD_UNAVAILABLE
    std::deque<std::thread> threads; //All of the threads in the pool
    std::deque<DRAGONWING::ThreadArgs> args; //All of the arguments given to each thread
    std::mutex mutex; //A mutex to protect inter-thread communication
    std::condition_variable cv; //To wake up the threads when a checkpoint (or restart request) has been reached

    int active_phase_1 = 0; //How many threads are in memory-heavy phase 1. Max value = DRAGONWING::Config::phase_1_max_threads
    int reached_checkpoint_1 = 0; //How many threads have completed phase 1 (most of the calculation)
    int reached_checkpoint_2 = 0; //How many threads have completed phase 2 (committing the updates)
    #endif
    int nthreads; //Total number of threads in the pool

    //Restart requests
    bool abort_requested = false;
    std::string restart_msg = "";
    
    std::atomic<std::size_t> fallback_counter;

public:
    ThreadPool(int N): nthreads(N){} //Create a pool with room for n threads
    void* launchParallel(DRAGON::Grid* grid, double dt); //Execute grid->advance(dt) on a new thread in the pool
    
    //Error Handling
    void reportFallback(const int weight, const int threshold); //increase fallback_counter, throw if at/above threshold
    void requestRestart(std::string msg = ""); //Something went wrong, request a restart
    std::string restartMsg(); //The message reported by the first thread to call requestRestart()
    
    //Checkpoints
    void reportCheckpoint1(); //Report that this thread's subgrid has computed (but not committed) its updated state
    void reportCheckpoint2(); //Report that this thread's subgrid has fully completed (and committed) its update
    //The following waitFor___ functions return false iff a restart has been requested
    bool waitForRelease(); //If active_phase_1 is at its limit, wait until someone reaches checkpoint 1 before this thread entersphase 1
    bool waitForCheckpoint1(); //Wait until all threads report checkpoint 1, or until any thread requests a restart
    bool waitForCompletion(); //Wait until all threads finish their work, or until any thread requests a restart
};
}

#endif
