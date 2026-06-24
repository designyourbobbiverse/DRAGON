//
//  Boundary/Periodic.cpp
//  DRAGON
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
    #ifdef MHD
    auto& _A = grid.getA();
    #endif

    if (faces & X){
        for(int j = j0 ; j < jn; j++){
            #ifdef MHD
            auto dAx = _A[nx,j].x - _A[0,j].x;
            #endif
            for(int g = 1; g <= ng; g++){
                grid[-g,j] = grid[nx-g,j];
                grid[nx-1+g,j] = grid[g-1,j];
                #ifdef MHD
                _A[-g,j] = _A[nx-g,j];
                _A[nx+g-1,j] = _A[g-1,j];
                _A[nx+g-1,j].x += dAx; //B must be periodic, Ax need not be periodic.
                #endif
            }
            #ifdef MHD
            _A[nx+ng,j] = _A[ng,j];
            _A[nx+ng,j].x += dAx;
            #endif
        }
    }
    if (faces & Y){
        for(int i = i0 ; i < in; i++){
            #ifdef MHD
            auto dAy = _A[i,ny].y - _A[i,0].y;
            #endif
            for(int g = 1; g <= ng; g++){
                grid[i,-g] = grid[i,ny-g];
                grid[i,ny-1+g] = grid[i,g-1];
                #ifdef MHD
                _A[i,-g] = _A[i,ny-g];
                _A[i,ny-1+g] = _A[i,g-1];
                _A[i,ny-1+g].y += dAy; //B must be periodic, Ay need not be periodic.
                #endif
            }
            #ifdef MHD
            _A[i,ny+ng] = _A[i,ng];
            _A[i,ny+ng].y += dAy; //B must be periodic, Ax need not be periodic.
            #endif
        }
    }
}
//MARK: 3D
//Set all ghosts to be the equivalent cell mod nx (or ny or nz)
void Boundary::Periodic::apply(Grid3D& grid) {
    int ng = grid.getGhosts(), nx = grid.getSizeX(), ny = grid.getSizeY(), nz = grid.getSizeZ();
    int i0 = (corners ? -ng : 0), in = (corners ? nx + ng : nx);
    int j0 = i0, jn = (corners ? ny + ng : ny);
    int k0 = i0, kn = (corners ? nz + ng : nz);
    #ifdef MHD
    auto& _A = grid.getA();
    #endif
    
    if (faces & X){
        for(int j = j0 ; j < jn; j++){
            for(int k = k0 ; k < kn; k++){
                #ifdef MHD
                auto dAx = _A[nx,j,k].x - _A[0,j,k].x;
                #endif
                for(int g = 1; g <= ng; g++){
                    grid[-g,j,k] = grid[nx-g,j,k];
                    grid[nx-1+g,j,k] = grid[g-1,j,k];
                    #ifdef MHD
                    _A[-g,j,k] = _A[nx-g,j,k];
                    _A[nx-1+g,j,k] = _A[g-1,j,k];
                    _A[nx-1+g,j,k].x += dAx; //B must be periodic, Ax need not be periodic.
                    #endif
                }
                #ifdef MHD
                _A[nx+ng,j,k] = _A[ng,j,k];
                _A[nx+ng,j,k] += dAx;
                #endif
            }
        }
    }
    if (faces & Y){
        for(int i = i0 ; i < in; i++){
            for(int k = k0 ; k < kn; k++){
                #ifdef MHD
                auto dAy = _A[i,ny,k].y - _A[i,0,k].y;
                #endif
                for(int g = 1; g <= ng; g++){
                    grid[i,-g,k] = grid[i,ny-g,k];
                    grid[i,ny-1+g,k] = grid[i,g-1,k];
                    #ifdef MHD
                    _A[i,-g,k] = _A[i,ny-g,k];
                    _A[i,ny-1+g,k] = _A[i,g-1,k];
                    _A[i,ny-1+g,k].y += dAy; //B must be periodic, Ay need not be periodic.
                    #endif
                }
                #ifdef MHD
                _A[i,ny+ng,k] = _A[i,ng,k];
                _A[i,ny+ng,k].y += dAy;
                #endif
            }
        }
    }
    if (faces & Z){
        for(int i = i0 ; i < in; i++){
            for(int j = j0 ; j < jn; j++){
                #ifdef MHD
                auto dAz = _A[i,j,nz].z - _A[i,j,0].z;
                #endif
                for(int g = 1; g <= ng; g++){
                    grid[i,j,-g] = grid[i,j,nz-g];
                    grid[i,j,nz-1+g] = grid[i,j,g-1];
                    #ifdef MHD
                    _A[i,j,-g] = _A[i,j,nz-g];
                    _A[i,j,nz-1+g] = _A[i,j,g-1];
                    _A[i,j,nz-1+g].z += dAz; //B must be periodic, Az need not be periodic.
                    #endif
                }
                #ifdef MHD
                _A[i,j,nz+ng] = _A[i,j,ng];
                #endif
            }
        }
    }
}
