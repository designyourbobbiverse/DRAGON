//
//  DragonWing.hpp
//  DRAGONWING
//
//  Created by Bobbie Markwick on 11/06/2026.
//

#ifndef DRAGONWING
#define DRAGONWING

class Advanceable;


namespace DRARGONWING{
    void* launchParallel(Advanceable* grid, double dt);

    void initialize(int nthreads);
    bool requestRestart(); //Returns true if genuinely in multithread mode
    void reportCheckpoint1();
    void reportCheckpoint2();
    bool waitForCheckpoint1(); //Returns false iff someone requested a restart
    bool waitForCheckpoint2(); //Returns false iff someone requested a restart
}


#endif
