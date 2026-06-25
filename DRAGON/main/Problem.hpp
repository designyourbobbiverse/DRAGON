//
//  Problem.hpp
//  DRAGON/main
//
//  Created by Bobbie Markwick on 24/06/2026.
//

#include "Grid.hpp"

namespace Problem{
    Advanceable& makeProblem();
    void cycleComplete(Advanceable& problem, int cycle, double time);
}
