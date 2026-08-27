//
//  WING_Atomic_Tests.cpp
//  DRAGON/Testing/Infrastructure/WING
//
//  Created by Bobbie Markwick on 09/06/2026.
//

#include "Testing.hpp"
#include "DragonWing.hpp"
using namespace DRAGON_Test;

//Things not tested:
//REUSE_AUX_GRIDS disabled
//Concurrency Stress Testing


void DRAGON_Test::verify_WING_fbcount_no_op(){
    for(int i=0; i<1000;i++){
        DRAGONWING::reportFallback(-1, -i);
        DRAGONWING::reportFallback(-1, 0);
        DRAGONWING::reportFallback(-1, 1);
        DRAGONWING::reportFallback(0, -i);
        DRAGONWING::reportFallback(0, 0);
        DRAGONWING::reportFallback(0, i);
    }
    DRAGONWING::resetFallbacks();
}
void DRAGON_Test::verify_WING_fbcount_reset(){
    DRAGONWING::resetFallbacks();
    DRAGONWING::reportFallback(1, 1);
    DRAGONWING::resetFallbacks();
    DRAGONWING::reportFallback(2, 2);
    DRAGONWING::resetFallbacks();
}
void DRAGON_Test::verify_WING_fbcount_single_no_throw(){
    for(int i=0; i<1000;i++){
        for(int j=0; j<i; j++){
            DRAGONWING::reportFallback(j, i);
            DRAGONWING::resetFallbacks();
        }
    }
}
void DRAGON_Test::verify_WING_fbcount_single_over_throw(){
    try{
        DRAGONWING::reportFallback(2, 1);
        assert(false); //This case is supposed to throw
    } catch(...) { }
    DRAGONWING::resetFallbacks();
    try{
        DRAGONWING::reportFallback(12,3);
        assert(false); //This case is supposed to throw
    } catch(...) { }
    DRAGONWING::resetFallbacks();
    try{
        DRAGONWING::reportFallback(21,12);
        assert(false); //This case is supposed to throw
    } catch(...) { }
    DRAGONWING::resetFallbacks();
}

void DRAGON_Test::verify_WING_fbcount_cumulative_throw(){
    DRAGONWING::reportFallback(1, 1);
    DRAGONWING::reportFallback(1, 2); //Verify that previous threshold doesn't matter, only the latest one
    try{
        DRAGONWING::reportFallback(1,2);
        assert(false); //This case is supposed to throw
    } catch(...) { }
    DRAGONWING::reportFallback(1, 2112); //But a previous throw doesn't doom it if we raise the threshold
    DRAGONWING::resetFallbacks();
    
}
