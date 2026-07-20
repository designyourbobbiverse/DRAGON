//
//  BoundaryTypes/Reflective.cpp
//  DRAGON/Boundary
//
//  Created by Bobbie Markwick on 28/05/2026.
//

#include "Boundary.hpp"
#include "Grid.hpp"

using namespace Boundary;

Boundary::Reflective::Reflective(int faces_, bool corners):  GhostFill(faces_, corners) {}
Boundary::Reflective::Reflective(std::string s, bool corners) : Reflective(face_mask(s),corners) {}

//MARK: 1D
//Set each ghost to be a reflection of the cell which is opposite the boundary from it
void Boundary::Reflective::apply(Grid1D& grid) {
    int ng = grid.getGhosts();
    if (faces & X_negative){
        for(int g = 1; g <= ng; g++){
            grid[-g] = grid[g-1];
            grid[-g].v.x *= -1;//Mirror Normal Velocity
            #ifdef MHD //Mirror Transverse Magnetic Fields
            grid[-g].B *= -1;
            grid[-g].B.x *= -1;
            #endif
        }
    }
    if (faces & X_positive){
        int nx = grid.getSize();
        for(int g = 1; g <= ng; g++){
            grid[nx-1+g] = grid[nx-g];
            grid[nx-1+g].v.x *= -1;//Mirror Normal Velocity
            #ifdef MHD //Mirror Transverse Magnetic Fields
            grid[nx-1+g].B *= -1;
            grid[nx-1+g].B.x *= -1;
            #endif
        }
    }
}

//MARK: 2D
//Set each ghost to be a reflection of the cell which is opposite the boundary from it
void Boundary::Reflective::apply(Grid2D& grid) {
    int ng = grid.getGhosts(), nx = grid.getSizeX(), ny = grid.getSizeY();
    int i0 = (corners ? -ng : 0), in = (corners ? nx + ng : nx);
    int j0 = i0, jn = (corners ? ny + ng : ny);


    if (faces & X_negative){
        for(int j = j0 ; j < jn; j++){
            for(int g = 1; g <= ng; g++){
                grid[-g,j] = grid[g-1,j];
                grid[-g,j].v.x *= -1;//Mirror Normal Velocity
                #ifdef MHD
                grid[-g,j].B *= -1;
                grid[-g,j].B.x *= -1;
                #endif
            }
        }
    }
    if (faces & X_positive){
        for(int j = j0 ; j < jn; j++){
            for(int g = 1; g <= ng; g++){
                grid[nx-1+g,j] = grid[nx-g,j];
                grid[nx-1+g,j].v.x *= -1;//Mirror Normal Velocity
                #ifdef MHD
                grid[nx-1+g,j].B *= -1;
                grid[nx-1+g,j].B.x *= -1;
                #endif
            }
        }
    }
    if (faces & Y_negative){
        for(int i = i0 ; i < in; i++){
            for(int g = 1; g <= ng; g++){
                grid[i,-g] = grid[i,g-1];
                grid[i,-g].v.y *= -1;//Mirror Normal Velocity
                #ifdef MHD
                grid[i,-g].B *= -1;
                grid[i,-g].B.y *= -1;
                #endif
            }
        }
    }
    if (faces & Y_positive){
        for(int i = i0 ; i < in; i++){
            for(int g = 1; g <= ng; g++){
                grid[i,ny-1+g] = grid[i,ny-g];
                grid[i,ny-1+g].v.y *= -1;//Mirror Normal Velocity
                #ifdef MHD
                grid[i,ny-1+g].B *= -1;
                grid[i,ny-1+g].B.y *= -1;
                #endif
            }
        }
    }
//MARK: 2D MHD
    #ifdef MHD //Mirror Transverse Magnetic Fields, preserve normal magnetic fields
    auto& _A = grid._A();
    // A has one more physical point per dimension than w. Transverse components are offset by 1/2 from w
    if (faces & X_negative){
        for(int j = j0 ; j <= jn; j++){
            for(int g = 1; g <= ng; g++){
                _A[-g,j] = _A[g,j]; //Copy the Transverse A field (reflected over 0 instead of -1/2)
                _A[-g,j].x = -_A[g-1,j].x; //Invert the Normal A field (which reflects over -1/2)
            }
        }
    }
    if (faces & X_positive){
        for(int j = j0 ; j <= jn; j++){
            for(int g = 1; g <= ng; g++){
                _A[nx+g,j] = _A[nx-g,j]; //Copy the Transverse A field (reflected over nx instead of nx+1/2)
                _A[nx-1+g,j].x = -_A[nx-g,j].x; //Invert the Normal A field (which reflects over nx+1/2)
            }
        }
    }
    if (faces & Y_negative){
        for(int i = i0 ; i <= in; i++){
            for(int g = 1; g <= ng; g++){
                _A[i,-g] = _A[i,g]; //Copy the Transverse A field (reflected over 0 instead of -1/2)
                _A[i,-g].y = -_A[i,g-1].y; //Invert the Normal A field (which reflects over -1/2)
            }
        }
    }
    if (faces & Y_positive){
        for(int i = i0 ; i <= in; i++){
            for(int g = 1; g <= ng; g++){
                _A[i,ny+g] = _A[i,ny-g]; //Copy the Transverse A field (reflected over ny instead of ny+1/2)
                _A[i,ny-1+g].y = -_A[i,ny-g].y; //Invert the Normal A field (which reflects over ny+1/2)
            }
        }
    }
    #endif

}
//MARK: 3D
//Set each ghost to be a reflection of the cell which is opposite the boundary from it
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
                    grid[-g,j,k].v.x *= -1;//Mirror Normal Velocity
                    #ifdef MHD
                    grid[-g,j,k].B *= -1;
                    grid[-g,j,k].B.x *= -1;
                    #endif
                }
            }
        }
    }
    if (faces & X_positive){
        for(int j = j0 ; j < jn; j++){
            for(int k = k0 ; k < kn; k++){
                for(int g = 1; g <= ng; g++){
                    grid[nx-1+g,j,k] = grid[nx-g,j,k];
                    grid[nx-1+g,j,k].v.x *= -1;//Mirror Normal Velocity
                    #ifdef MHD
                    grid[nx-1+g,j,k].B *= -1;
                    grid[nx-1+g,j,k].B.x *= -1;
                    #endif
                }
            }
        }
    }
    if (faces & Y_negative){
        for(int i = i0 ; i < in; i++){
            for(int k = k0 ; k < kn; k++){
                for(int g = 1; g <= ng; g++){
                    grid[i,-g,k] = grid[i,g-1,k];
                    grid[i,-g,k].v.y *= -1;//Mirror Normal Velocity
                    #ifdef MHD
                    grid[i,-g,k].B *= -1;
                    grid[i,-g,k].B.y *= -1;
                    #endif
                }
            }
        }
    }
    if (faces & Y_positive){
        for(int i = i0 ; i < in; i++){
            for(int k = k0 ; k < kn; k++){
                for(int g = 1; g <= ng; g++){
                    grid[i,ny-1+g,k] = grid[i,ny-g,k];
                    grid[i,ny-1+g,k].v.y *= -1;//Mirror Normal Velocity
                    #ifdef MHD
                    grid[i,ny-1+g,k].B *= -1;
                    grid[i,ny-1+g,k].B.y *= -1;
                    #endif
                }
            }
        }
    }
    if (faces & Z_negative){
        for(int i = i0 ; i < in; i++){
            for(int j = j0 ; j < jn; j++){
                for(int g = 1; g <= ng; g++){
                    grid[i,j,-g] = grid[i,j,g-1];
                    grid[i,j,-g].v.z *= -1;//Mirror Normal Velocity
                    #ifdef MHD
                    grid[i,j,-g].B *= -1;
                    grid[i,j,-g].B.z *= -1;
                    #endif
                }
            }
        }
    }
    if (faces & Z_positive){
        for(int i = i0 ; i < in; i++){
            for(int j = j0 ; j < jn; j++){
                for(int g = 1; g <= ng; g++){
                    grid[i,j,nz-1+g] = grid[i,j,nz-g];
                    grid[i,j,nz-1+g].v.z *= -1;//Mirror Normal Velocity
                    #ifdef MHD
                    grid[i,j,nz-1+g].B *= -1;
                    grid[i,j,nz-1+g].B.z *= -1;
                    #endif
                }
            }
        }
    }
//MARK: 3D MHD
    #ifdef MHD //Mirror Transverse Magnetic Fields, preserve normal magnetic fields
    auto& _A = grid._A();
    // A has one more physical point per dimension than w. Transverse components are offset by 1/2 from w
    if (faces & X_negative){
        for(int j = j0 ; j <= jn; j++){
            for(int k = k0 ; k <= kn; k++){
                for(int g = 1; g <= ng; g++){
                    _A[-g,j,k] = _A[g,j,k]; //Copy the Transverse A field (reflected over 0 instead of -1/2)
                    _A[-g,j,k].x = -_A[g-1,j,k].x; //Invert the Normal A field (which reflects over -1/2)
                }
            }
        }
    }
    if (faces & X_positive){
        for(int j = j0 ; j <= jn; j++){
            for(int k = k0 ; k <= kn; k++){
                for(int g = 1; g <= ng; g++){
                    _A[nx+g,j,k] = _A[nx-g,j,k]; //Copy the Transverse A field (reflected over nx instead of nx+1/2)
                    _A[nx-1+g,j,k].x = -_A[nx-g,j,k].x; //Invert the Normal A field (which reflects over nx+1/2)
                }
            }
        }
    }
    if (faces & Y_negative){
        for(int i = i0 ; i <= in; i++){
            for(int k = k0 ; k <= kn; k++){
                for(int g = 1; g <= ng; g++){
                    _A[i,-g,k] = _A[i,g,k]; //Copy the Transverse A field (reflected over 0 instead of -1/2)
                    _A[i,-g,k].y = -_A[i,g-1,k].y; //Invert the Normal A field (which reflects over -1/2)
                }
            }
        }
    }
    if (faces & Y_positive){
        for(int i = i0 ; i <= in; i++){
            for(int k = k0 ; k <= kn; k++){
                for(int g = 1; g <= ng; g++){
                    _A[i,ny+g,k] = _A[i,ny-g,k]; //Copy the Transverse A field (reflected over ny instead of ny+1/2)
                    _A[i,ny-1+g,k].y = -_A[i,ny-g,k].y; //Invert the Normal A field (which reflects over ny+1/2)
                }
            }
        }
    }
    if (faces & Z_negative){
        for(int i = i0 ; i <= in; i++){
            for(int j = j0 ; j <= jn; j++){
                for(int g = 1; g <= ng; g++){
                    _A[i,j,-g] = _A[i,j,g]; //Copy the Transverse A field (reflected over 0 instead of -1/2)
                    _A[i,j,-g].z = -_A[i,j,g-1].z; //Invert the Normal A field (which reflects over -1/2)
                }
            }
        }
    }
    if (faces & Z_positive){
        for(int i = i0 ; i <= in; i++){
            for(int j = j0 ; j <= jn; j++){
                for(int g = 1; g <= ng; g++){
                    _A[i,j,nz+g] = _A[i,j,nz-g]; //Copy the Transverse A field (reflected over nz instead of nz+1/2)
                    _A[i,j,nz-1+g].z = -_A[i,j,nz-g].z; //Invert the Normal A field (which reflects over nz+1/2)
                }
            }
        }
    }
    #endif
}
