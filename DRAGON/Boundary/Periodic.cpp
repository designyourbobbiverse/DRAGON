//
//  Boundary/Periodic.cpp
//  DRAGON
//
//  Created by Bobbie Markwick on 28/05/2026.
//

#include "Boundary.hpp"

using namespace Boundary;

//MARK: Constructors
Boundary::Periodic::Periodic(int faces_, bool corners) : GhostFill(
    ((faces_ & X) ? X : 0) | ((faces_ & Y) ? Y : 0) | ((faces_ & Z) ? Z : 0)
, corners){} //Ensure that periodic boundaries have matching edges
Boundary::Periodic::Periodic(int faces_) : Periodic(faces_,true) {}
Boundary::Periodic::Periodic(): Periodic(X|Y|Z) {}


//MARK: 1D
void Boundary::Periodic::apply(Grid1D& grid) const {
    if((faces & X) == 0) return;
    int ng = grid.getGhosts(), nx = grid.getSize();
    for(int g = 1; g <= ng; g++){
        grid[-g] = grid[nx-g];
        grid[nx-1+g] = grid[g-1];
    }
}

//MARK: 2D
void Boundary::Periodic::apply(Grid2D& grid) const {
    int ng = grid.getGhosts(), nx = grid.getSizeX(), ny = grid.getSizeY();
    int i0 = (corners ? -ng : 0), in = (corners ? nx + ng : nx);
    int j0 = i0, jn = (corners ? ny + ng : ny);

    if (faces & X){
        for(int j = j0 ; j < jn; j++){
            for(int g = 1; g <= ng; g++){
                grid[-g,j] = grid[nx-g,j];
                grid[nx-1+g,j] = grid[g-1,j];
            }
        }
    }
    if (faces & Y){
        for(int i = i0 ; i < in; i++){
            for(int g = 1; g <= ng; g++){
                grid[i,-g] = grid[i,ny-g];
                grid[i,ny-1+g] = grid[i,g-1];
            }
        }
    }
}
//MARK: 3D
void Boundary::Periodic::apply(Grid3D& grid) const {
    int ng = grid.getGhosts(), nx = grid.getSizeX(), ny = grid.getSizeY(), nz = grid.getSizeZ();
    int i0 = (corners ? -ng : 0), in = (corners ? nx + ng : nx);
    int j0 = i0, jn = (corners ? ny + ng : ny);
    int k0 = i0, kn = (corners ? nz + ng : nz);

    if (faces & X){
        for(int j = j0 ; j < jn; j++){
            for(int k = k0 ; k < kn; k++){
                for(int g = 1; g <= ng; g++){
                    grid[-g,j,k] = grid[nx-g,j,k];
                    grid[nx-1+g,j,k] = grid[g-1,j,k];
                }
            }
        }
    }
    if (faces & Y){
        for(int i = i0 ; i < in; i++){
            for(int k = k0 ; k < kn; k++){
                for(int g = 1; g <= ng; g++){
                    grid[i,-g,k] = grid[i,ny-g,k];
                    grid[i,ny-1+g,k] = grid[i,g-1,k];
                }
            }
        }
    }
    if (faces & Z){
        for(int i = i0 ; i < in; i++){
            for(int j = j0 ; j < jn; j++){
                for(int g = 1; g <= ng; g++){
                    grid[i,j,-g] = grid[i,j,nz-g];
                    grid[i,j,nz-1+g] = grid[i,j,g-1];
                }
            }
        }
    }
}
