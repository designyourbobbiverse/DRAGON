//
//  Ignore.cpp
//  DRAGON
//
//  Created by Bobbie Markwick on 12/06/2026.
//

#include "Boundary.hpp"

using namespace Boundary;

Boundary::Ignore::Ignore(int faces_, bool corners):  GhostFill(faces_, corners) {}
Boundary::Ignore::Ignore(std::string s, bool corners) : Ignore(face_mask(s),corners) {}

//This class does nothing
void Boundary::Ignore::apply(Grid1D &grid) {}
void Boundary::Ignore::apply(Grid2D &grid) {}
void Boundary::Ignore::apply(Grid3D &grid) {}
