//
//  pthread_parallel.cpp
//  DRAGON
//
//  Created by Bobbie Markwick on 11/06/2026.
//

#include "Dragon_Wing.hpp"
#include <pthread.h>
#include <iostream>
#include <deque>
#include "Grid.hpp"

namespace DRARGONWING{
    struct ThreadArgs_Grid2D {
        Grid2D* grid;
        double dt;
    };

    void* advance2D(void* arg) {
        ThreadArgs_Grid2D* args = static_cast<ThreadArgs_Grid2D*>(arg);
        args->grid->advance(args->dt,false);
        return nullptr;
    }

    std::deque<pthread_t> threads;
    std::deque<ThreadArgs_Grid2D> args;

}

void* DRARGONWING::launchParallel(Grid2D* grid, double dt){
    args.push_back({grid,dt});
    pthread_t thread;

    int err = pthread_create(&thread, nullptr, advance2D, &args.back());
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

