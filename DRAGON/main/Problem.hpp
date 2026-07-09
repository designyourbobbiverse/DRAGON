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
    void afterCycle(Grid& problem, int cycle, double t);
    void beforeCycle(Grid& problem, int cycle, double t);
    void problemComplete(Grid& problem, double t);
}
