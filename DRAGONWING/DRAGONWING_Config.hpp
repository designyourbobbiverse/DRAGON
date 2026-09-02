//
//  DRAGONWING_Config.hpp
//  DRAGONWING
//
//  Created by Bobbie Markwick on 10/07/2026.
//

#ifndef DRAGONWING_Config_h
#define DRAGONWING_Config_h

namespace DRAGONWING::Config{


#define REUSE_AUX_GRIDS //Reuse Auxiliary Grids instead of reallocating them each time


//phase_1_max_threads controls the number of grids allowed to be in phase 1 of the unsplit update simultaneously
//Phase 1 is computation heavy and may use over a dozen auxiliary grids per grid
//A smaller number uses less memory, a larger number is faster (to a point)
//If you are worried about memory, a good starting point for this parameter is typically to match the number of cores on your machine
#ifndef TESTMODE //The test mode value is set later, and is smaller than you typically want
constexpr int phase_1_max_threads = -1; //Use a nonpositive number means that there is no limit
#endif





//******************************************************************//
//*DO NOT TOUCH anything below this line (but feel free to read it)*//
//******************************************************************//

#ifdef TESTMODE
constexpr int phase_1_max_threads = 5; //Deliberately set low for testing
#endif
}
#endif
