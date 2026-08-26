//
//  FallbackCounter.cpp
//  DRAGONWING/Atomics
//
//  Created by Bobbie Markwick on 26/08/2026.
//

#include "DragonWing.hpp"

#include "DRAGONWING_Config.hpp"
#ifndef MULTITHREAD_UNAVAILABLE
#include <atomic>    //For Fallback Counting
#endif
#include <stdexcept> //For throwing exceptions

#ifdef MULTITHREAD_UNAVAILABLE
static int fallback_counter;
#else
static std::atomic<std::size_t> fallback_counter;
#endif

void DRAGONWING::reportFallback(int weight, int threshold){
    if(weight <= 0 || threshold < 0) return; //User doesn't want track this type of fallback
    
    #ifdef MULTITHREAD_UNAVAILABLE
    fallback_counter += weight;
    if(fallback_counter > threshold) throw std::runtime_error("Too many fallback behaviours reported");
    #else
    std::size_t prev = fallback_counter.fetch_add(weight, std::memory_order_relaxed);
    if(prev+weight > threshold) throw std::runtime_error("Too many fallback behaviours reported");
    #endif
}

void DRAGONWING::resetFallbacks(){ fallback_counter = 0; }

