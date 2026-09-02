//
//  FallbackCounter.cpp
//  DRAGONWING/Atomics
//
//  Created by Bobbie Markwick on 26/08/2026.
//

#include "DragonWing.hpp"

#include <atomic>    //For Fallback Counting

#include <stdexcept> //For throwing exceptions

static std::atomic<std::size_t> fallback_counter;

void DRAGONWING::reportFallback(int weight, int threshold){
    if(weight <= 0 || threshold < 0) return; //User doesn't want track this type of fallback
    
    std::size_t prev = fallback_counter.fetch_add(weight, std::memory_order_relaxed);
    if(prev+weight > threshold) throw std::runtime_error("Too many fallback behaviours reported");
}

void DRAGONWING::resetFallbacks(){ fallback_counter = 0; }

