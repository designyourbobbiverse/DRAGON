//
//  DragonWing.hpp
//  DRAGONWING
//
//  Created by Bobbie Markwick on 11/06/2026.
//

#ifndef DRAGON_WING_hpp
#define DRAGON_WING_hpp
#include <string>

class Grid;


namespace DRAGONWING{
    void* launchParallel(Grid* grid, double dt);

    void initialize(int nthreads);
    bool requestRestart(std::string msg = ""); //Returns true if genuinely in multithread mode
    std::string restartMsg();
    void reportCheckpoint1();
    void reportCheckpoint2();
    bool waitForCheckpoint1(); //Returns false iff someone requested a restart
    bool waitForCheckpoint2(); //Returns false iff someone requested a restart
}


#endif
