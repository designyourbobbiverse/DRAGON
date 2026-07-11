//
//  DRAGONWING_Config.hpp
//  DRAGONWING
//
//  Created by Bobbie Markwick on 10/07/2026.
//


namespace DRAGONWING::CONFIG{


#define REUSE_AUX_GRIDS //Reuse Auxiliary Grids instead of reallocating them each time


//This controls the number of grids allowed to be in phase 1 of the unsplit update simultaneously
//Phase 1 is computation heavy and uses as many as 19 auxiliary grids per grid
//A smaller number uses less memory, a larger number is faster (to a point)
//A good starting point for this parameter is typically to match the number of cores on your machine
constexpr int phase_1_max_threads = 10000;





//Enable this if for some reason you need to compile DRAGON for a system without multithreading support
//#define MULTITHREAD_UNAVAILABLE

}


