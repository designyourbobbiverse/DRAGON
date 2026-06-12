//
//  Dragon_Wing.hpp
//  DRAGON
//
//  Created by Bobbie Markwick on 11/06/2026.
//

class Grid2D;

namespace DRARGONWING{
    void* launchParallel(Grid2D* grid, double dt);
    void synchronize(void* thread = nullptr);
}
