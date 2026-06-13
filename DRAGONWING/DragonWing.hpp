//
//  Dragon_Wing.hpp
//  DRAGON
//
//  Created by Bobbie Markwick on 11/06/2026.
//

#ifndef DRAGONWING
#define DRAGONWING

class AMRGrid1D;
class AMRGrid2D;
class AMRGrid3D;


namespace DRARGONWING{
    void* launchParallel(AMRGrid1D* grid, double dt);
    void* launchParallel(AMRGrid2D* grid, double dt);
    void* launchParallel(AMRGrid3D* grid, double dt);
    void synchronize(void* thread = nullptr);
}


#endif
