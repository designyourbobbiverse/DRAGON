//
//  Boundary/BoundarySet.cpp
//  DRAGON
//
//  Created by Bobbie Markwick on 29/05/2026.
//

#include "Boundary.hpp"
#include <algorithm>


using namespace Boundary;
int Boundary::BoundaryType::get_faces() const { return faces; }

void Boundary::BoundarySet::apply(Grid1D &grid) const{
    for (auto& b : boundaries)  b->apply(grid);
}
void Boundary::BoundarySet::apply(Grid2D &grid) const{
    for (auto& b : boundaries)  b->apply(grid);
}
void Boundary::BoundarySet::apply(Grid3D &grid) const{
    for (auto& b : boundaries)  b->apply(grid);
}
