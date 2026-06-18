//
//  pthread_parallel.cpp
//  DRAGON
//
//  Created by Bobbie Markwick on 11/06/2026.
//

#include "DragonWing.hpp"
#include <pthread.h>
#include <iostream>
#include <deque>
#include "AMRGrid.hpp"


namespace DRARGONWING{
    struct ThreadArgs_Grid1D {
        AMRGrid1D* grid;
        double dt;
    };
    struct ThreadArgs_Grid2D {
        AMRGrid2D* grid;
        double dt;
    };
    struct ThreadArgs_Grid3D {
        AMRGrid3D* grid;
        double dt;
    };


    void* advance1D(void* arg) {
        ThreadArgs_Grid1D* args = static_cast<ThreadArgs_Grid1D*>(arg);
        args->grid->advance(args->dt,false);
        return nullptr;
    }
    void* advance2D(void* arg) {
        ThreadArgs_Grid2D* args = static_cast<ThreadArgs_Grid2D*>(arg);
        args->grid->advance(args->dt,false);
        return nullptr;
    }
    void* advance3D(void* arg) {
        ThreadArgs_Grid3D* args = static_cast<ThreadArgs_Grid3D*>(arg);
        args->grid->advance(args->dt,false);
        return nullptr;
    }

    std::deque<pthread_t> threads;
    std::deque<ThreadArgs_Grid1D> args1D;
    std::deque<ThreadArgs_Grid2D> args2D;
    std::deque<ThreadArgs_Grid3D> args3D;

}

void* DRARGONWING::launchParallel(AMRGrid1D* grid, double dt){
    args1D.push_back({grid,dt});
    pthread_t thread;

    int err = pthread_create(&thread, nullptr, advance1D, &args1D.back());
    if (err != 0) {
        std::cerr << "Failed to create thread\n";
        return nullptr;
    }
    
    threads.push_back(thread);
    
    return &threads.back();
}

void* DRARGONWING::launchParallel(AMRGrid2D* grid, double dt){
    args2D.push_back({grid,dt});
    pthread_t thread;

    int err = pthread_create(&thread, nullptr, advance2D, &args2D.back());
    if (err != 0) {
        std::cerr << "Failed to create thread\n";
        return nullptr;
    }
    
    threads.push_back(thread);
    
    return &threads.back();
}

void* DRARGONWING::launchParallel(AMRGrid3D* grid, double dt){
    args3D.push_back({grid,dt});
    pthread_t thread;

    int err = pthread_create(&thread, nullptr, advance3D, &args3D.back());
    if (err != 0) {
        std::cerr << "Failed to create thread\n";
        return nullptr;
    }
    
    threads.push_back(thread);
    
    return &threads.back();
}

void DRARGONWING::synchronize(void* thread){
    if(thread == nullptr) {
        for (pthread_t& thread : threads) pthread_join(thread, nullptr);
    } else {
        pthread_t _thread = *static_cast<pthread_t*>(thread);
        pthread_join(_thread, nullptr);
    }
}

