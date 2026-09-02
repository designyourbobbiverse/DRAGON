//
//  WING_Test_Main.cpp
//  DRAGON/Testing/Infrastructure/WING
//
//  Created by Bobbie Markwick on 09/06/2026.
//

#include "Testing.hpp"
#include "DragonWing.hpp"

#include <iostream> //For std::cout
using namespace DRAGON_Test;


void DRAGON_Test::verify_WING(bool output){
    if (output) std::cout << "DRAGONWING (Memory Management):\n";
    if (output) std::cout << "- Allocation & Dimensions:\n";
    if (output) std::cout << "\t- Primitive Arrays: ";
    verify_WING_alloc_prim_dimensions();
    if (output) std::cout << "Passed\n";
    if (output) std::cout << "\t- Conservative Arrays: ";
    verify_WING_alloc_flux_dimensions();
    if (output) std::cout << "Passed\n";
    if (output) std::cout << "\t- vec3 Arrays: ";
    verify_WING_alloc_vec3_dimensions();
    if (output) std::cout << "Passed\n";
    
    #ifdef REUSE_AUX_GRIDS
    if (output) std::cout << "- Array Reuse: ";
    verify_WING_array_reuse();
    if (output) std::cout << "Passed\n";
    #endif
    
    if (output) std::cout << "- ArrayGuard:\n";
    if (output) std::cout << "\t- get() == []: ";
    verify_WING_arrayguard_get();
    if (output) std::cout << "Passed\n";
    if (output) std::cout << "\t- Release on destruct: ";
    verify_WING_scope_release();
    if (output) std::cout << "Passed\n";
    if (output) std::cout << "\t- Multiple release() test: ";
    verify_WING_double_release();
    if (output) std::cout << "Passed\n";
    
    if (output) std::cout << "- Active arrays survive purge: ";
    verify_WING_array_purge_survival();
    if (output) std::cout << "Passed\n";
    if (output) std::cout << "All Memory Management Tests Passed.\n\n";

    
    
    if (output) std::cout << "DRAGONWING (Atomics):\n";
    if (output) std::cout << "- Reset: ";
    verify_WING_fbcount_reset();
    if (output) std::cout << "Passed\n";
    if (output) std::cout << "- Zero op: ";
    verify_WING_fbcount_no_op();
    if (output) std::cout << "Passed\n";
    if (output) std::cout << "- Single weight<=threshold doesn't throw: ";
    verify_WING_fbcount_single_no_throw();
    if (output) std::cout << "Passed\n";
    if (output) std::cout << "- Single weight>threshold does throw: ";
    verify_WING_fbcount_single_over_throw();
    if (output) std::cout << "Passed\n";
    if (output) std::cout << "- Cumulative calls throw: ";
    verify_WING_fbcount_cumulative_throw();
    if (output) std::cout << "Passed\n";
    if (output) std::cout << "All Atomics Tests Passed.\n\n";
    
    
    
    if (output) std::cout << "DRAGONWING (Multithreading):\n";
    if (output) std::cout << "- Clear Restart Message: ";
    verify_WING_restart_message();
    if (output) std::cout << "\n"; //This particular test should print out Passed itself.
    if (output) std::cout << "- Launch parallel: ";
    verify_WING_launchParallel_calls_advance();
    if (output) std::cout << "Passed\n";
    if (output) std::cout << "- Error catch: ";
    verify_WING_launchParallel_catches();
    if (output) std::cout << "\n"; //This particular test should print out Passed itself.
    if (output) std::cout << "- Checkpoints: ";
    verify_WING_checkpoints();
    if (output) std::cout << "Passed\n";
    if (output) std::cout << "- Restart blocks phase 2: ";
    verify_WING_checkpoints_restart();
    if (output) std::cout << "Passed\n";
    if (output) std::cout << "All Multithreading Tests Passed.\n\n";

}
