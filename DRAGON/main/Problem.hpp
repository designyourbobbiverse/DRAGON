//
//  Problem.hpp
//  DRAGON/main
//
//  Created by Bobbie Markwick on 24/06/2026.
//

#ifndef Problem_hpp
#define Problem_hpp

#include "Grid.hpp"

namespace Problem{
    DRAGON::Grid& makeProblem();

    DRAGON::PrimitiveState initialFluidState(double x, double y, double z);
    DRAGON::vec3 initialMagneticPotential(double x, double y, double z);
    void completeProblemInit(DRAGON::Grid& grid);


    void afterCycle(DRAGON::Grid& problem, int cycle, double t);
    void beforeCycle(DRAGON::Grid& problem, int cycle, double t);
    void problemComplete(DRAGON::Grid& problem, double t);
}

#endif
