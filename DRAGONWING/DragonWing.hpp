//
//  Dragon_Wing.hpp
//  DRAGON
//
//  Created by Bobbie Markwick on 11/06/2026.
//

#ifndef DRAGONWING
#define DRAGONWING

class AMRGrid2D;

namespace DRARGONWING{
    void* launchParallel(AMRGrid2D* grid, double dt);
    void synchronize(void* thread = nullptr);
}


#endif
