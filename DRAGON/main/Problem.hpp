//
//  Problem.hpp
//  DRAGON/main
//
//  Created by Bobbie Markwick on 24/06/2026.
//

#include "Grid.hpp"

namespace Problem{
    Grid& makeProblem();

    PrimitiveState initialFluidState(double x, double y, double z);
    vec3 initialMagneticPotential(double x, double y, double z);
    void completeProblemInit(Grid& grid);


    void afterCycle(Grid& problem, int cycle, double t);
    void beforeCycle(Grid& problem, int cycle, double t);
    void problemComplete(Grid& problem, double t);
}
