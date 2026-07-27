//
//  BoundaryTypes/Periodic.cpp
//  DRAGON/Boundary
//
//  Created by Bobbie Markwick on 28/05/2026.
//

#include "Boundary.hpp"
#include "Grid.hpp"

using namespace Boundary;

//MARK: Constructors
Boundary::Periodic::Periodic(int faces_, bool corners) : //Ensure that periodic boundaries have matching edges
    GhostFill( ((faces_ & X) ? X : 0) | ((faces_ & Y) ? Y : 0) | ((faces_ & Z) ? Z : 0), corners){}
Boundary::Periodic::Periodic(std::string s, bool corners) : Periodic(face_mask(s),corners) {}

//MARK: 1D
//Set all ghosts to be the equivalent cell mod nx
void Boundary::Periodic::apply(Grid1D& grid) {
    if((faces & X) == 0) return;
    int ng = grid.getGhosts(), nx = grid.getSize();
    for(int g = 1; g <= ng; g++){
        grid[-g] = grid[nx-g];
        grid[nx-1+g] = grid[g-1];
    }
}

//MARK: 2D
//Set all ghosts to be the equivalent cell mod nx (or ny)
void Boundary::Periodic::apply(Grid2D& grid) {
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
//MARK: 2D MHD
    #ifdef MHD
    auto& _A = grid._A();
    // A has one more physical point per dimension than w.
    if (faces & X){
        for(int j = j0 ; j <= jn; j++){
            auto dA = _A[nx,j] - _A[0,j]; //B must be periodic, but A has some gauge freedom

            for(int g = 1; g <= ng; g++){
                _A[-g,j] = _A[nx-g,j] - dA;
                _A[nx+g,j] = _A[g,j] + dA;
            }
        }
    }
    if (faces & Y){
        for(int i = i0 ; i <= in; i++){
            auto dA = _A[i,ny] - _A[i,0]; //B must be periodic, but A has some gauge freedom

            for(int g = 1; g <= ng; g++){
                _A[i,-g] = _A[i,ny-g] - dA;
                _A[i,ny+g] = _A[i,g] + dA;
            }
        }
    }
    #endif
}
//MARK: 3D
//Set all ghosts to be the equivalent cell mod nx (or ny or nz)
void Boundary::Periodic::apply(Grid3D& grid) {
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
//MARK: 3D MHD
    #ifdef MHD
    auto& _B = grid._B();
    // Faces have one more physical point per dimension than bodies
    if (faces & X){
        for(int j = j0 ; j <= jn; j++){
            for(int k = k0 ; k <= kn; k++){
                for(int g = 1; g <= ng; g++){
                    _B[-g,j,k] = _B[nx-g,j,k];
                    _B[nx-1+g,j,k] = _B[g-1,j,k];
                }
            }
        }
    }
    if (faces & Y){
        for(int i = i0 ; i <= in; i++){
            for(int k = k0 ; k <= kn; k++){
                for(int g = 1; g <= ng; g++){
                    _B[i,-g,k] = _B[i,ny-g,k];
                    _B[i,ny-1+g,k] = _B[i,g-1,k];
                }
            }
        }
    }
    if (faces & Z){
        for(int i = i0 ; i <= in; i++){
            for(int j = j0 ; j <= jn; j++){
                for(int g = 1; g <= ng; g++){
                    _B[i,j,-g] = _B[i,j,nz-g];
                    _B[i,j,nz-1+g] = _B[i,j,g-1];
                }
            }
        }
    }
    #endif
}
