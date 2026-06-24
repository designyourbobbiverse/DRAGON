//
//  Boundary/Fixed.cpp
//  DRAGON
//
//  Created by Bobbie Markwick on 28/05/2026.
//

#include "Boundary.hpp"
#include "Grid.hpp"

using namespace Boundary;

//MARK: Constructors
Boundary::Fixed::Fixed(PrimitiveState w_, int faces_, bool corners):  GhostFill(faces_, corners), state(w_) {}
Boundary::Fixed::Fixed(PrimitiveState w_,std::string s,bool corners) : Fixed(w_,face_mask(s),corners) {}

//MARK: 1D
//All applicable ghost cells are set to equal this->state
void Boundary::Fixed::apply(Grid1D& grid) {
    int ng = grid.getGhosts();
    if (faces & X_negative){
        for(int g = 1; g <= ng; g++){
            grid[-g] = state;
        }
    }
    if (faces & X_positive){
        int nx = grid.getSize();
        for(int g = 1; g <= ng; g++){
            grid[nx-1+g] = state;
        }
    }
}

//MARK: 2D
//All applicable ghost cells are set to equal this->state
void Boundary::Fixed::apply(Grid2D& grid) {
    //Calculate the bounds ahead of time
    int ng = grid.getGhosts(), nx = grid.getSizeX(), ny = grid.getSizeY();
    int i0 = (corners ? -ng : 0), in = (corners ? nx + ng : nx);
    int j0 = i0, jn = (corners ? ny + ng : ny);


    if (faces & X_negative){
        for(int j = j0 ; j < jn; j++){
            for(int g = 1; g <= ng; g++){
                grid[-g,j] = state;
            }
        }
    }
    if (faces & X_positive){
        for(int j = j0 ; j < jn; j++){
            for(int g = 1; g <= ng; g++){
                grid[nx-1+g,j] = state;
            }
        }
    }
    if (faces & Y_negative){
        for(int i = i0 ; i < in; i++){
            for(int g = 1; g <= ng; g++){
                grid[i,-g] = state;
            }
        }
    }
    if (faces & Y_positive){
        for(int i = i0 ; i < in; i++){
            for(int g = 1; g <= ng; g++){
                grid[i,ny-1+g] = state;
            }
        }
    }
//MARK: 2D MHD
    #ifdef MHD
    auto& _A = grid.getA();
    // A has one more physical point per dimension than w.
    if (faces & X_negative){
        for(int j = j0 ; j <= jn; j++){
            for(int g = 1; g <= ng; g++){
                _A[-g,j].z =_A[-g+1,j].z  + state.B.y * grid.dx;
            }
        }
    }
    if (faces & X_positive){
        for(int j = j0 ; j <= jn; j++){
            for(int g = 1; g <= ng; g++){
                _A[nx+g,j].z = _A[nx-1+g,j].z - state.B.y * grid.dx;
            }
        }
    }
    if (faces & Y_negative){
        for(int i = i0 ; i <= in; i++){
            for(int g = 1; g <= ng; g++){
                _A[i,-g].z =_A[i,-g+1].z  - state.B.x * grid.dy;
            }
        }
    }
    if (faces & Y_positive){
        for(int i = i0 ; i <= in; i++){
            for(int g = 1; g <= ng; g++){
                _A[i,ny+g].z = _A[i,ny-1+g].z + state.B.x  * grid.dy;
            }
        }
    }
    #endif
}
//MARK: 3D
//All applicable ghost cells are set to equal this->state
void Boundary::Fixed::apply(Grid3D& grid) {
    //Calculate the bounds ahead of time
    int ng = grid.getGhosts(), nx = grid.getSizeX(), ny = grid.getSizeY(), nz = grid.getSizeZ();
    int i0 = (corners ? -ng : 0), in = (corners ? nx + ng : nx);
    int j0 = i0, jn = (corners ? ny + ng : ny);
    int k0 = i0, kn = (corners ? nz + ng : nz);
    
    if (faces & X_negative){
        for(int j = j0; j < jn; j++){
            for(int k = k0; k < kn; k++){
                for(int g = 1; g <= ng; g++){
                    grid[-g,j,k] = state;
                }
            }
        }
    }
    if (faces & X_positive){
        for(int j = j0; j < jn; j++){
            for(int k = k0; k < kn; k++){
                for(int g = 1; g <= ng; g++){
                    grid[nx-1+g,j,k] = state;
                }
            }
        }
    }
    if (faces & Y_negative){
        for(int i = i0; i < in; i++){
            for(int k = k0; k < kn; k++){
                for(int g = 1; g <= ng; g++){
                    grid[i,-g,k] = state;
                }
            }
        }
    }
    if (faces & Y_positive){
        for(int i = i0; i < in; i++){
            for(int k = k0; k < kn; k++){
                for(int g = 1; g <= ng; g++){
                    grid[i,ny-1+g,k] = state;
                }
            }
        }
    }
    if (faces & Z_negative){
        for(int i = i0; i < in; i++){
            for(int j = j0; j < jn; j++){
                for(int g = 1; g <= ng; g++){
                    grid[i,j,-g] = state;
                }
            }
        }
    }
    if (faces & Z_positive){
        for(int i = i0; i < in; i++){
            for(int j = j0; j < jn; j++){
                for(int g = 1; g <= ng; g++){
                    grid[i,j,nz-1+g] = state;
                }
            }
        }
    }
    
//MARK: 3D MHD
#ifdef MHD
    auto& _A = grid.getA();
    // A has one more physical point per dimension than w.

    double dx = grid.dx, dy = grid.dy, dz = grid.dz;
    if (faces & X_negative){
        for(int j = j0; j <= jn; j++){
            for(int k = k0; k <= kn; k++){
                for(int g = 1; g <= ng; g++){
                    _A[-g,j,k].x = _A[-g+1,j,k].x;
                    double dAx_dz = k==kn ? 0 : (_A[-g+1,j,k+1].x - _A[-g+1,j,k].x) / dz;
                    double dAx_dy = j == jn ? 0 : (_A[-g+1,j+1,k].x - _A[-g+1,j,k].x) / dy;
                    
                    _A[-g,j,k].y = _A[-g+1,j,k].y - (state.B.z+dAx_dy) * dx;
                    _A[-g,j,k].z = _A[-g+1,j,k].z + (state.B.y-dAx_dz) * dx;
                }
            }
        }
    }
    if (faces & X_positive){
        for(int j = j0; j <= jn; j++){
            for(int k = k0; k <= kn; k++){
                for(int g = 1; g <= ng; g++){
                    _A[nx+g-1,j,k].x = _A[nx+g-2,j,k].x;
                    double dAx_dz = k == kn ? 0 : (_A[nx+g-2,j,k+1].x - _A[nx+g-2,j,k].x) / dz;
                    double dAx_dy = j == jn ? 0 : (_A[nx+g-2,j+1,k].x - _A[nx+g-2,j,k].x) / dy;
                    
                    _A[nx+g,j,k].y = _A[nx+g-1,j,k].y + (state.B.z+dAx_dy) * dx;
                    _A[nx+g,j,k].z = _A[nx+g-1,j,k].z - (state.B.y-dAx_dz) * dx;
                }
            }
        }
    }
    if (faces & Y_negative){
        for(int i = i0; i <= in; i++){
            for(int k = k0; k <= kn; k++){
                for(int g = 1; g <= ng; g++){
                    _A[i,-g,k].y = _A[i,-g+1,k].y;
                    double dAy_dz = k == kn ? 0 : (_A[i,-g+1,k+1].y - _A[i,-g+1,k].y) / dz;
                    double dAy_dx = i == in ? 0 : (_A[i+1,-g+1,k].y - _A[i,-g+1,k].y) / dx;

                    _A[i,-g,k].x = _A[i,-g+1,k].x + (state.B.z-dAy_dx) * dy;
                    _A[i,-g,k].z = _A[i,-g+1,k].z - (state.B.x+dAy_dz) * dy;
                }
            }
        }
    }
    if (faces & Y_positive){
        for(int i = i0; i <= in; i++){
            for(int k = k0; k <= kn; k++){
                for(int g = 1; g <= ng; g++){
                    _A[i,ny+g-1,k].y = _A[i,ny+g-2,k].y;
                    double dAy_dz = k == kn ? 0 : (_A[i,ny+g-2,k+1].y - _A[i,ny+g-2,k].y) / dz;
                    double dAy_dx = i == in ? 0 : (_A[i+1,ny+g-2,k].y - _A[i,ny+g-2,k].y) / dx;

                    _A[i,ny+g,k].x = _A[i,ny+g-1,k].x - (state.B.z-dAy_dx) * dy;
                    _A[i,ny+g,k].z = _A[i,ny+g-1,k].z + (state.B.x+dAy_dz) * dy;
                }
            }
        }
    }
    if (faces & Z_negative){
        for(int i = i0; i <= in; i++){
            for(int j = j0; j <= jn; j++){
                for(int g = 1; g <= ng; g++){
                    _A[i,j,-g].z = _A[i,j,-g+1].z;
                    double dAz_dy = j == jn ? 0 : (_A[i,j+1,-g+1].z - _A[i,j,-g+1].z) / dy;
                    double dAz_dx = i == in ? 0 : (_A[i+1,j,-g+1].z - _A[i,j,-g+1].z) / dx;

                    _A[i,j,-g].y = _A[i,j,-g+1].y + (state.B.x-dAz_dy) * dz;
                    _A[i,j,-g].x = _A[i,j,-g+1].x - (state.B.y+dAz_dx) * dz;
                }
            }
        }
    }
    if (faces & Z_positive){
        for(int i = i0; i <= in; i++){
            for(int j = j0; j <= jn; j++){
                for(int g = 1; g <= ng; g++){
                    _A[i,j,nz+g-1].z = _A[i,j,nz+g-2].z;
                    double dAz_dy = j == jn ? 0 : (_A[i,j+1,nz+g-2].z - _A[i,j,nz+g-2].z) / dy;
                    double dAz_dx = i == in ? 0 : (_A[i+1,j,nz+g-2].z - _A[i,j,nz+g-2].z) / dx;

                    _A[i,j,nz+g].y = _A[i,j,nz+g-1].y - (state.B.x-dAz_dy) * dz;
                    _A[i,j,nz+g].x = _A[i,j,nz+g-1].x + (state.B.y+dAz_dx) * dz;
                }
            }
        }
    }
#endif
}
