//
//  Problem.hpp
//  DRAGON/main
//
//  Created by Bobbie Markwick on 24/06/2026.
//

#include "Grid.hpp"

namespace Problem{
    Grid& makeProblem();
    void initializeProblem(Grid& grid);
    void cycleComplete(Grid& problem, int cycle);
}
