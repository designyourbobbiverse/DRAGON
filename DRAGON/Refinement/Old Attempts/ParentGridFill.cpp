//
//  ParentGridFill.cpp
//  DRAGON
//
//  Created by Bobbie Markwick on 12/06/2026.
//

#include "ParentGridFill.hpp"

using namespace Boundary;

ParentGridFill::ParentGridFill(AMRGrid2D *p, int x, int y): GhostFill(X|Y, true), parent(p), x_offset(x), y_offset(y) {}
    
void ParentGridFill::apply(Grid1D& grid) {
    
}
void ParentGridFill::apply(Grid2D& grid) {
    int g = grid.getGhosts(), nx = grid.getSizeX(), ny = grid.getSizeY();
    for(int i = -g; i < nx + g; i++){
        for(int j = -g; j < 0; j++) grid[i,j] = (*parent)[i+x_offset,j+y_offset];
        for(int j = ny; j < ny+g; j++) grid[i,j] = (*parent)[i+x_offset,j+y_offset];
    }
    for(int j = 0; j < ny; j++){
        for(int i = -g; i < 0; i++) grid[i,j] = (*parent)[i+x_offset,j+y_offset];
        for(int i = nx; i < nx+g; i++) grid[i,j] = (*parent)[i+x_offset,j+y_offset];
    }
}
void ParentGridFill::apply(Grid3D& grid) {
    
}
