//
//  Boundary/BoundarySet.cpp
//  DRAGON
//
//  Created by Bobbie Markwick on 29/05/2026.
//

#include "Boundary.hpp"
#include <algorithm>


using namespace Boundary;
int ::GhostFill::get_faces() const { return faces; }

void Boundary::Boundaries::apply(Grid1D &grid) const{
    for (auto& b : boundaries)  b->apply(grid);
}
void Boundary::Boundaries::apply(Grid2D &grid) const{
    for (auto& b : boundaries)  b->apply(grid);
}
void Boundary::Boundaries::apply(Grid3D &grid) const{
    for (auto& b : boundaries)  b->apply(grid);
}

int Boundary::face_mask(std::string s){
    int mask = 0;
    s += " ";
    for(auto i = s.begin(); i+1 < s.end(); i++){
        char current = *i, next = *(i + 1);
        if(current == 'X' || current == 'x'){
            switch(next){
            case '+': mask |= X_positive; break;
            case '-': mask |= X_negative; break;
            default: mask |= X;
            }
        } else if(current == 'Y' || current == 'y'){
            switch(next){
            case '+': mask |= Y_positive; break;
            case '-': mask |= Y_negative; break;
            default: mask |= Y;
            }
        } else if(current == 'Z' || current == 'z'){
            switch(next){
            case '+': mask |= Z_positive; break;
            case '-': mask |= Z_negative; break;
            default: mask |= Z;
            }
        }
    }
    return mask;
}
