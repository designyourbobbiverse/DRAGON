//
//  pthread_parallel.cpp
//  DRAGONWING
//
//  Created by Bobbie Markwick on 11/06/2026.
//

#include "DragonWing.hpp"
#include <thread>
#include <mutex>
#include <condition_variable>
#include <iostream>
#include <deque>
#include "Grid.hpp"


namespace {
    struct ThreadArgs{
        Grid* grid;
        double dt;
    };


    void* advance(void* arg) {
        ThreadArgs* args = static_cast<ThreadArgs*>(arg);
        args->grid->advance(args->dt,false);
        return nullptr;
    }

    std::deque<pthread_t> threads;
    std::deque<ThreadArgs> args;
    int nthreads = 0;

    std::string restart_msg;

    std::mutex mutex;
    std::condition_variable cv;
    int reached_checkpoint_1 = 0;
    int reached_checkpoint_2 = 0;
    bool abort_requested = false;

}

//MARK: Launching
void DRAGONWING::initialize(int _nthreads){
    nthreads = _nthreads;
    threads.clear();
    args.clear();
    reached_checkpoint_1 = 0;
    reached_checkpoint_2 = 0;
    abort_requested = false;
}


void* DRAGONWING::launchParallel(Grid* grid, double dt){
    if (nthreads == 0) {
        std::cerr << "Did not initialize DRAGONWING before attempting to launch parallel task";
        return nullptr;
    }
    
    args.push_back({grid,dt});
    pthread_t thread;

    int err = pthread_create(&thread, nullptr, advance, &args.back());
    if (err != 0) {
        std::cerr << "Failed to create thread\n";
        return nullptr;
    }
    
    threads.push_back(thread);
    
    return &threads.back();
}


//MARK: Checkpoint reporting
void DRAGONWING::reportCheckpoint1(){
    if(nthreads == 0) return; //Single thread mode

    std::unique_lock lock(mutex);
    bool done = (++reached_checkpoint_1 == nthreads);
    lock.unlock();
    if (done) cv.notify_all();
}
void DRAGONWING::reportCheckpoint2(){
    if(nthreads == 0) return; //Single thread mode

    std::unique_lock lock(mutex);
    bool done = (++reached_checkpoint_2 == nthreads);
    lock.unlock();
    if (done) cv.notify_all();
}

bool DRAGONWING::requestRestart(std::string msg){
    std::unique_lock lock(mutex);
    restart_msg = msg;
    if(nthreads == 0) return false; //Single thread mode
    abort_requested = true;
    lock.unlock();
    cv.notify_all();
    return true;
}

std::string DRAGONWING::restartMsg(){
    std::lock_guard lock(mutex);
    auto msg = restart_msg;
    restart_msg = "";
    return msg;
}




//MARK: Synchronization
bool DRAGONWING::waitForCheckpoint1(){
    if(nthreads == 0) return true; //Single thread mode

    std::unique_lock lock(mutex);
    //Wait until everybody else is done
    cv.wait(lock, [&] { return abort_requested || reached_checkpoint_1 == nthreads; });
    return !abort_requested;
}
bool DRAGONWING::waitForCheckpoint2(){
    if(nthreads == 0) return true; //Single thread mode
    
    //Wait until everybody else is done
    std::unique_lock lock(mutex);
    cv.wait(lock, [&] { return abort_requested || reached_checkpoint_2 == nthreads; });
    bool success = !abort_requested;
    lock.unlock();
    //Wait for all the threads to clean up
    for (pthread_t& thread : threads) pthread_join(thread, nullptr);
    return success;
}


