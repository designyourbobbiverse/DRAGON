//
//  ThreadPool.cpp
//  DRAGONWING/Multithreading
//
//  Created by Bobbie Markwick on 11/06/2026.
//

#include "DragonWing.hpp"

#include "DRAGONWING_Config.hpp"

#include "Hydro/Grid.hpp"   //For grid.advance()

#include <iostream>         //For std::cerr

//MARK: Threadpool access
//The thread pool associated with the current thread
#ifndef MULTITHREAD_UNAVAILABLE
thread_local DRAGONWING::ThreadPool* current_thread_pool;
#else
static DRAGONWING::ThreadPool* current_thread_pool;
#endif

//Call the corresponding function on this thread's threadpool, if one exists
void DRAGONWING::reportCheckpoint1(){
    if (!current_thread_pool) return; //Single thread mode
    current_thread_pool->reportCheckpoint1();
}
void DRAGONWING::reportCheckpoint2(){
    if (!current_thread_pool) return; //Single thread mode
    current_thread_pool->reportCheckpoint2();
}
bool DRAGONWING::requestRestart(std::string msg){
    if (!current_thread_pool) return false; //Single thread mode
    current_thread_pool->requestRestart(msg);
    return true;
}
bool DRAGONWING::waitForRelease(){
    if (!current_thread_pool) return true; //Single thread mode
    return current_thread_pool->waitForRelease();
}
bool DRAGONWING::waitForCheckpoint1(){
    if (!current_thread_pool) return true; //Single thread mode
    return current_thread_pool->waitForCheckpoint1();
}

//MARK: Thread launching
#ifndef MULTITHREAD_UNAVAILABLE
void* DRAGONWING::ThreadPool::launchParallel(DRAGON::Grid* grid, double dt){
    args.push_back({grid,dt}); //Set up the thread arguments (grid + dt)
    ThreadArgs* thread_args = &args.back();
    
    try {
        threads.emplace_back([thread_args, this]{//Run the code on a new thread
            current_thread_pool = this; //Give the thread access to the pool
            thread_args->grid->advance(thread_args->dt, false); //Do the actual work
        });
    } catch (const std::system_error&) {
        std::cerr << "Failed to create thread\n";
        args.pop_back();
        return nullptr;
    }
    return &threads.back();
}
#else
void* DRAGONWING::ThreadPool::launchParallel(Grid* grid, double dt){
    current_thread_pool = this;
    grid->advance(dt, false);
    current_thread_pool = nullptr;
    return nullptr;
}
#endif

//MARK: Error Handling
std::string DRAGONWING::ThreadPool::restartMsg(){
    #ifndef MULTITHREAD_UNAVAILABLE
    std::lock_guard lock(mutex);
    #endif
    //Return the current message + clear the field to make room for the next error
    auto msg = restart_msg;
    restart_msg = "";
    return msg;
}

#ifndef MULTITHREAD_UNAVAILABLE
//MARK: Checkpoint reporting
void DRAGONWING::ThreadPool::reportCheckpoint1(){
    std::unique_lock lock(mutex);
    active_phase_1--; //Allow someone else to start doing memory-intensive work
    ++reached_checkpoint_1;
    lock.unlock();
    cv.notify_all();
}
void DRAGONWING::ThreadPool::reportCheckpoint2(){
    std::unique_lock lock(mutex);
    bool done = (++reached_checkpoint_2 == nthreads);
    lock.unlock();
    if (done) cv.notify_all();
}
void DRAGONWING::ThreadPool::requestRestart(std::string msg){
    std::unique_lock lock(mutex);
    if (restart_msg.size() < 1) restart_msg = msg;
    abort_requested = true;
    lock.unlock();
    cv.notify_all();
}

//MARK: Synchronization
bool DRAGONWING::ThreadPool::waitForRelease(){
    std::unique_lock lock(mutex);
    if (Config::phase_1_max_threads <= 0) return !abort_requested; //nonpositive phase_1_max_threads => no limit on threads in phase 1
    //Wait until there is room to proceed
    cv.wait(lock, [&] { return abort_requested || (active_phase_1 < Config::phase_1_max_threads); });
    if (!abort_requested)  active_phase_1++;
    return !abort_requested;
}
bool DRAGONWING::ThreadPool::waitForCheckpoint1(){
    std::unique_lock lock(mutex);
    //Wait until everybody else is done
    cv.wait(lock, [&] { return abort_requested || reached_checkpoint_1 == nthreads; });
    return !abort_requested;
}

bool DRAGONWING::ThreadPool::waitForCompletion(){
    //Wait until everybody else is done
    std::unique_lock lock(mutex);
    cv.wait(lock, [&] { return abort_requested || reached_checkpoint_2 == nthreads; });
    bool success = !abort_requested;
    lock.unlock();
    //Wait for all the threads to clean up
    for (auto& thread : threads) {
        if (thread.joinable()) thread.join();
    }
    return success;
}

#else //MARK: Single Threaded Stubs
void DRAGONWING::ThreadPool::reportCheckpoint1(){ }
void DRAGONWING::ThreadPool::reportCheckpoint2(){ }
void DRAGONWING::ThreadPool::requestRestart(std::string msg){
    if (restart_msg.size() < 1) restart_msg = msg;
    abort_requested = true;
}
bool DRAGONWING::ThreadPool::waitForRelease(){ return !abort_requested; }
bool DRAGONWING::ThreadPool::waitForCheckpoint1(){ return !abort_requested; }
bool DRAGONWING::ThreadPool::waitForCompletion(){ return !abort_requested; }
#endif
