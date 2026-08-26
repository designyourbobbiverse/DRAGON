//
//  GhostFill.cpp
//  DRAGON/Boundary
//
//  Created by Bobbie Markwick on 29/05/2026.
//

#include "Boundary/Boundary.hpp"

#include "Hydro/Grid.hpp" //For dispatch to 1D/2D/3D
using namespace DRAGON;

int GhostFill::get_faces() const { return faces; }

//Convert a String such as "X-" or "YZ" into a bitmask
int Boundary::face_mask(std::string s){
    int mask = 0;
    s += " ";
    for (auto i = s.begin(); i+1 < s.end(); i++) {
        char current = *i, next = *(i + 1);
        if (current == 'X' || current == 'x') {
            switch(next) {
            case '+': mask |= X_positive; break;
            case '-': mask |= X_negative; break;
            default: mask |= X;
            }
        } else if (current == 'Y' || current == 'y') {
            switch(next) {
            case '+': mask |= Y_positive; break;
            case '-': mask |= Y_negative; break;
            default: mask |= Y;
            }
        } else if (current == 'Z' || current == 'z') {
            switch(next) {
            case '+': mask |= Z_positive; break;
            case '-': mask |= Z_negative; break;
            default: mask |= Z;
            }
        }
    }
    return mask;
}


void GhostFill::apply(Grid& grid){
    Grid3D* grid3D = dynamic_cast<Grid3D*>(&grid);
    if (grid3D) {
        apply(*grid3D);
        return;
    }
    Grid2D* grid2D = dynamic_cast<Grid2D*>(&grid);
    if (grid2D) {
        apply(*grid2D);
        return;
    }
    Grid1D* grid1D = dynamic_cast<Grid1D*>(&grid);
    if (grid1D) {
        apply(*grid1D);
        return;
    }
}
