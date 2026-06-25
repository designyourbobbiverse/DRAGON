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
    void synchronize(void* thread = nullptr);
}


#endif
