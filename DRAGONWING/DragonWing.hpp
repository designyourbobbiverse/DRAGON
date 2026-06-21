//
//  Dragon_Wing.hpp
//  DRAGON
//
//  Created by Bobbie Markwick on 11/06/2026.
//

#ifndef DRAGONWING
#define DRAGONWING

class DistGrid1D;
class DistGrid2D;
class DistGrid3D;


namespace DRARGONWING{
    void* launchParallel(DistGrid1D* grid, double dt);
    void* launchParallel(DistGrid2D* grid, double dt);
    void* launchParallel(DistGrid3D* grid, double dt);
    void synchronize(void* thread = nullptr);
}


#endif
