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
    //Grid object construction
    DRAGON::Grid& makeProblem();

    //User-specified initial state
    DRAGON::PrimitiveState initialFluidState(double x, double y, double z);
    DRAGON::vec3 initialMagneticPotential(double x, double y, double z);
    void completeProblemInit(DRAGON::Grid& grid);

    //Initialization machinery (init.cpp instead of problem.cpp, generally not user-edited)
    void load(DRAGON::Grid& problem, double& time, int& cycle);
    void init(DRAGON::Grid1D& grid);
    void init(DRAGON::Grid2D& grid);
    void init(DRAGON::Grid3D& grid);


    //pre/post-Cycle and post-completion processing
    void beforeCycle(DRAGON::Grid& problem, int cycle, double t);
    void afterCycle(DRAGON::Grid& problem, int cycle, double t);
    void problemComplete(DRAGON::Grid& problem, double t);

}

#endif
