//
//  Boundary/Set.cpp
//  DRAGON
//
//  Created by Bobbie Markwick on 29/05/2026.
//

#include "Boundary.hpp"
#include <algorithm>


using namespace Boundary;

int Boundary::BoundaryType::get_faces(){ return faces; }

//MARK: Initialization

template<typename... Bs> Boundary::BoundarySet::BoundarySet(Bs&&... bs) {
    (addBoundary(std::forward<Bs>(bs)), ...);

    //Check for uncovered faces
    int uncoverdFaces = X | Y | Z;
    for (auto& b : boundaries){
        uncoverdFaces &= ~b->get_faces();
    }
    addBoundary(Outflow(uncoverdFaces));
}

template<typename B> void Boundary::BoundarySet::addBoundary(B&& b) {
      static_assert(std::derived_from<std::decay_t<B>, BoundaryType>);
      boundaries.push_back(std::make_unique<std::decay_t<B>>(std::forward<B>(b)));
}

//MARK: Applying Conditions
void Boundary::BoundarySet::apply(Grid1D &grid) const{
    for (auto& b : boundaries)  b->apply(grid);
}
void Boundary::BoundarySet::apply(Grid2D &grid) const{
    for (auto& b : boundaries)  b->apply(grid);
}
void Boundary::BoundarySet::apply(Grid3D &grid) const{
    for (auto& b : boundaries)  b->apply(grid);
}
