//
//  Reflective.cpp
//  Hydro Code
//
//  Created by Bobbie Markwick on 28/05/2026.
//

#include "Boundary.hpp"

using namespace Boundary;

//MARK: Constructors
Boundary::Reflective::Reflective(int faces_, bool corners):  BoundaryType(faces_, corners) {}
Boundary::Reflective::Reflective(int faces_) : Reflective(faces_,true) {}
Boundary::Reflective::Reflective(): Reflective(X|Y|Z) {}

//MARK: 1D
void Boundary::Reflective::apply(Grid1D& grid) const {
    int ng = grid.getGhosts();
    if (faces & X_negative){
        for(int g = 1; g <= ng; g++){
            grid[-g] = grid[g-1];
            grid[-g].vx *= -1;
        }
    }
    if (faces & X_positive){
        int nx = grid.getSize();
        for(int g = 1; g <= ng; g++){
            grid[nx-1+g] = grid[nx-g];
            grid[nx-1+g].vx *= -1;
        }
    }
}

//MARK: 2D
void Boundary::Reflective::apply(Grid2D& grid) const {
    int ng = grid.getGhosts(), nx = grid.getSizeX(), ny = grid.getSizeY();
    int i0 = (corners ? -ng : 0), in = (corners ? nx + ng : nx);
    int j0 = i0, jn = (corners ? ny + ng : ny);

    if (faces & X_negative){
        for(int j = j0 ; j < jn; j++){
            for(int g = 1; g <= ng; g++){
                grid[-g,j] = grid[g-1,j];
                grid[-g,j].vx *= -1;
            }
        }
    }
    if (faces & X_positive){
        for(int j = j0 ; j < jn; j++){
            for(int g = 1; g <= ng; g++){
                grid[nx-1+g,j] = grid[nx-g,j];
                grid[nx-1+g,j].vx *= -1;
            }
        }
    }
    if (faces & Y_negative){
        for(int i = i0 ; i < in; i++){
            for(int g = 1; g <= ng; g++){
                grid[i,-g] = grid[i,g-1];
                grid[i,-g].vy *= -1;
            }
        }
    }
    if (faces & Y_positive){
        for(int i = i0 ; i < in; i++){
            for(int g = 1; g <= ng; g++){
                grid[i,ny-1+g] = grid[i,ny-g];
                grid[i,ny-1+g].vy *= -1;
            }
        }
    }
}
//MARK: 3D
void Boundary::Reflective::apply(Grid3D& grid) const {
    int ng = grid.getGhosts(), nx = grid.getSizeX(), ny = grid.getSizeY(), nz = grid.getSizeZ();
    int i0 = (corners ? -ng : 0), in = (corners ? nx + ng : nx);
    int j0 = i0, jn = (corners ? ny + ng : ny);
    int k0 = i0, kn = (corners ? nz + ng : nz);

    if (faces & X_negative){
        for(int j = j0 ; j < jn; j++){
            for(int k = k0 ; k < kn; k++){
                for(int g = 1; g <= ng; g++){
                    grid[-g,j,k] = grid[g-1,j,k];
                    grid[-g,j,k].vx *= -1;
                }
            }
        }
    }
    if (faces & X_positive){
        for(int j = j0 ; j < jn; j++){
            for(int k = k0 ; k < kn; k++){
                for(int g = 1; g <= ng; g++){
                    grid[nx-1+g,j,k] = grid[nx-g,j,k];
                    grid[nx-1+g,j,k].vx *= -1;
                }
            }
        }
    }
    if (faces & Y_negative){
        for(int i = i0 ; i < in; i++){
            for(int k = k0 ; k < kn; k++){
                for(int g = 1; g <= ng; g++){
                    grid[i,-g,k] = grid[i,g-1,k];
                    grid[i,-g,k].vy *= -1;
                }
            }
        }
    }
    if (faces & Y_positive){
        for(int i = i0 ; i < in; i++){
            for(int k = k0 ; k < kn; k++){
                for(int g = 1; g <= ng; g++){
                    grid[i,ny-1+g,k] = grid[i,ny-g,k];
                    grid[i,ny-1+g,k].vy *= -1;
                }
            }
        }
    }
    if (faces & Z_negative){
        for(int i = i0 ; i < in; i++){
            for(int j = j0 ; j < jn; j++){
                for(int g = 1; g <= ng; g++){
                    grid[i,j,-g] = grid[i,j,g-1];
                    grid[i,j,-g].vz *= -1;
                }
            }
        }
    }
    if (faces & Z_positive){
        for(int i = i0 ; i < in; i++){
            for(int j = j0 ; j < jn; j++){
                for(int g = 1; g <= ng; g++){
                    grid[i,j,nz-1+g] = grid[i,j,nz-g];
                    grid[i,j,nz-1+g].vz *= -1;
                }
            }
        }
    }
}
