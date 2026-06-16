//
//  Boundary/Reflective.cpp
//  DRAGON
//
//  Created by Bobbie Markwick on 28/05/2026.
//

#include "Boundary.hpp"
#include "Grid.hpp"

using namespace Boundary;

Boundary::Reflective::Reflective(int faces_, bool corners):  GhostFill(faces_, corners) {}
Boundary::Reflective::Reflective(std::string s, bool corners) : Reflective(face_mask(s),corners) {}

//MARK: 1D
void Boundary::Reflective::apply(Grid1D& grid) {
    int ng = grid.getGhosts();
    if (faces & X_negative){
        for(int g = 1; g <= ng; g++){
            grid[-g] = grid[g-1];
            grid[-g].v.x *= -1;
            grid[-g].B.x *= -1;
        }
    }
    if (faces & X_positive){
        int nx = grid.getSize();
        for(int g = 1; g <= ng; g++){
            grid[nx-1+g] = grid[nx-g];
            grid[nx-1+g].v.x *= -1;
            grid[nx-1+g].B.x *= -1;
        }
    }
}

//MARK: 2D
void Boundary::Reflective::apply(Grid2D& grid) {
    int ng = grid.getGhosts(), nx = grid.getSizeX(), ny = grid.getSizeY();
    int i0 = (corners ? -ng : 0), in = (corners ? nx + ng : nx);
    int j0 = i0, jn = (corners ? ny + ng : ny);

    if (faces & X_negative){
        for(int j = j0 ; j < jn; j++){
            for(int g = 1; g <= ng; g++){
                grid[-g,j] = grid[g-1,j];
                grid[-g,j].v.x *= -1;
                grid[-g,j].B.x *= -1;
            }
        }
    }
    if (faces & X_positive){
        for(int j = j0 ; j < jn; j++){
            for(int g = 1; g <= ng; g++){
                grid[nx-1+g,j] = grid[nx-g,j];
                grid[nx-1+g,j].v.x *= -1;
                grid[nx-1+g,j].B.x *= -1;
            }
        }
    }
    if (faces & Y_negative){
        for(int i = i0 ; i < in; i++){
            for(int g = 1; g <= ng; g++){
                grid[i,-g] = grid[i,g-1];
                grid[i,-g].v.y *= -1;
                grid[i,-g].B.y *= -1;
            }
        }
    }
    if (faces & Y_positive){
        for(int i = i0 ; i < in; i++){
            for(int g = 1; g <= ng; g++){
                grid[i,ny-1+g] = grid[i,ny-g];
                grid[i,ny-1+g].v.y *= -1;
                grid[i,ny-1+g].B.y *= -1;
            }
        }
    }
}
//MARK: 3D
void Boundary::Reflective::apply(Grid3D& grid) {
    int ng = grid.getGhosts(), nx = grid.getSizeX(), ny = grid.getSizeY(), nz = grid.getSizeZ();
    int i0 = (corners ? -ng : 0), in = (corners ? nx + ng : nx);
    int j0 = i0, jn = (corners ? ny + ng : ny);
    int k0 = i0, kn = (corners ? nz + ng : nz);

    if (faces & X_negative){
        for(int j = j0 ; j < jn; j++){
            for(int k = k0 ; k < kn; k++){
                for(int g = 1; g <= ng; g++){
                    grid[-g,j,k] = grid[g-1,j,k];
                    grid[-g,j,k].v.x *= -1;
                    grid[-g,j,k].B.x *= -1;
                }
            }
        }
    }
    if (faces & X_positive){
        for(int j = j0 ; j < jn; j++){
            for(int k = k0 ; k < kn; k++){
                for(int g = 1; g <= ng; g++){
                    grid[nx-1+g,j,k] = grid[nx-g,j,k];
                    grid[nx-1+g,j,k].v.x *= -1;
                    grid[nx-1+g,j,k].B.x *= -1;
                }
            }
        }
    }
    if (faces & Y_negative){
        for(int i = i0 ; i < in; i++){
            for(int k = k0 ; k < kn; k++){
                for(int g = 1; g <= ng; g++){
                    grid[i,-g,k] = grid[i,g-1,k];
                    grid[i,-g,k].v.y *= -1;
                    grid[i,-g,k].B.y *= -1;
                }
            }
        }
    }
    if (faces & Y_positive){
        for(int i = i0 ; i < in; i++){
            for(int k = k0 ; k < kn; k++){
                for(int g = 1; g <= ng; g++){
                    grid[i,ny-1+g,k] = grid[i,ny-g,k];
                    grid[i,ny-1+g,k].v.y *= -1;
                    grid[i,ny-1+g,k].B.y *= -1;
                }
            }
        }
    }
    if (faces & Z_negative){
        for(int i = i0 ; i < in; i++){
            for(int j = j0 ; j < jn; j++){
                for(int g = 1; g <= ng; g++){
                    grid[i,j,-g] = grid[i,j,g-1];
                    grid[i,j,-g].v.z *= -1;
                    grid[i,j,-g].B.z *= -1;
                }
            }
        }
    }
    if (faces & Z_positive){
        for(int i = i0 ; i < in; i++){
            for(int j = j0 ; j < jn; j++){
                for(int g = 1; g <= ng; g++){
                    grid[i,j,nz-1+g] = grid[i,j,nz-g];
                    grid[i,j,nz-1+g].v.z *= -1;
                    grid[i,j,nz-1+g].B.z *= -1;
                }
            }
        }
    }
}
