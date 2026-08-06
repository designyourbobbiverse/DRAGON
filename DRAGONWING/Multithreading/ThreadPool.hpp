//
//  ThreadPool.h
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
#include <deque>    //For thread tracking

namespace DRAGONWING{
struct ThreadArgs{
    DRAGON::Grid* grid;
    double dt;
};
}
#endif

namespace DRAGONWING {

class ThreadPool{
private:
    #ifndef MULTITHREAD_UNAVAILABLE
    std::deque<std::thread> threads;
    std::deque<DRAGONWING::ThreadArgs> args;
    std::mutex mutex;
    std::condition_variable cv;

    int active_phase_1 = 0;
    int reached_checkpoint_1 = 0;
    int reached_checkpoint_2 = 0;
    #endif
    int nthreads;

    bool abort_requested = false;
    std::string restart_msg = "";

public:
    ThreadPool(int nthreads): nthreads(nthreads) {}
    void* launchParallel(DRAGON::Grid* grid, double dt);
    
    void requestRestart(std::string msg = "");
    std::string restartMsg();
    
    void reportCheckpoint1();
    void reportCheckpoint2();
    bool waitForRelease();
    bool waitForCheckpoint1(); //Returns false iff someone requested a restart
    bool waitForCompletion(); //Returns false iff someone requested a restart
};
}

#endif
