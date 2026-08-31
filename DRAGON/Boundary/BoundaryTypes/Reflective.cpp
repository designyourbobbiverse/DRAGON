//
//  BoundaryTypes/Reflective.cpp
//  DRAGON/Boundary
//
//  Created by Bobbie Markwick on 28/05/2026.
//

#include "Boundary/Boundary.hpp"

#include "Hydro/Grid.hpp" //For filling ghost cells in the grid

using namespace DRAGON;
using namespace Boundary;

Boundary::Reflective::Reflective(int faces_, bool corners):  GhostFill(faces_, corners){}
Boundary::Reflective::Reflective(std::string s, bool corners) : Reflective(face_mask(s),corners){}

//MARK: 1D
//Set each ghost to be a reflection of the cell which is opposite the boundary from it
void Boundary::Reflective::apply(Grid1D& grid){
    int ng = grid.getGhosts();
    auto& q = grid.passives();

    if (faces & X_negative) {
        for (int g = 1; g <= ng; g++) {
            grid[-g] = grid[g-1];
            grid[-g].v.x *= -1;//Mirror Normal Velocity
            #ifdef MHD //Mirror Transverse Magnetic Fields
            grid[-g].B *= -1;
            grid[-g].B.x *= -1;
            #endif
        }
        q[-1] = q[0]; //Copy passives (only first ghost)
    }
    if (faces & X_positive) {
        int nx = grid.getSize();
        for (int g = 1; g <= ng; g++) {
            grid[nx-1+g] = grid[nx-g];
            grid[nx-1+g].v.x *= -1;//Mirror Normal Velocity
            #ifdef MHD //Mirror Transverse Magnetic Fields
            grid[nx-1+g].B *= -1;
            grid[nx-1+g].B.x *= -1;
            #endif
        }
        q[nx] = q[nx-1]; //Copy passives (only first ghost)
    }
}

//MARK: 2D
//Set each ghost to be a reflection of the cell which is opposite the boundary from it
void Boundary::Reflective::apply(Grid2D& grid){
    int ng = grid.getGhosts(), nx = grid.getSizeX(), ny = grid.getSizeY();
    int i0 = (corners ? -ng : 0), in = (corners ? nx + ng : nx);
    int j0 = i0, jn = (corners ? ny + ng : ny);

    auto& q = grid.passives();

    if (faces & X_negative) {
        for (int j = j0 ; j < jn; j++) {
            for (int g = 1; g <= ng; g++) {
                grid[-g,j] = grid[g-1,j];
                grid[-g,j].v.x *= -1;//Mirror Normal Velocity
                #ifdef MHD
                grid[-g,j].B *= -1;
                grid[-g,j].B.x *= -1;
                #endif
            }
            if (j >= -1 && j <= ny) q[-1,j] = q[0,j]; //Copy passives (only first ghost)
        }
    }
    if (faces & X_positive) {
        for (int j = j0 ; j < jn; j++) {
            for (int g = 1; g <= ng; g++) {
                grid[nx-1+g,j] = grid[nx-g,j];
                grid[nx-1+g,j].v.x *= -1;//Mirror Normal Velocity
                #ifdef MHD
                grid[nx-1+g,j].B *= -1;
                grid[nx-1+g,j].B.x *= -1;
                #endif
            }
            if (j >= -1 && j <= ny) q[nx,j] = q[nx-1,j]; //Copy passives (only first ghost)
        }
    }
    if (faces & Y_negative) {
        for (int i = i0 ; i < in; i++) {
            for (int g = 1; g <= ng; g++) {
                grid[i,-g] = grid[i,g-1];
                grid[i,-g].v.y *= -1;//Mirror Normal Velocity
                #ifdef MHD
                grid[i,-g].B *= -1;
                grid[i,-g].B.y *= -1;
                #endif
            }
            if (i >= -1 && i <= nx) q[i,-1] = q[i,0]; //Copy passives (only first ghost)
        }
    }
    if (faces & Y_positive) {
        for (int i = i0 ; i < in; i++) {
            for (int g = 1; g <= ng; g++) {
                grid[i,ny-1+g] = grid[i,ny-g];
                grid[i,ny-1+g].v.y *= -1;//Mirror Normal Velocity
                #ifdef MHD
                grid[i,ny-1+g].B *= -1;
                grid[i,ny-1+g].B.y *= -1;
                #endif
            }
            if (i >= -1 && i <= nx) q[i,ny] = q[i,ny-1]; //Copy passives (only first ghost)
        }
    }
//MARK: 2D MHD
    #ifdef MHD //Mirror Transverse Magnetic Fields, preserve normal magnetic fields
    auto& _B = grid._B(); // B lives on faces, which have one more physical point per dimension than w.
    if (faces & X_negative) {
        for (int j = j0 ; j <= jn; j++) {
            for (int g = 1; g <= ng; g++) {
                _B[-g,j] = _B[g-1,j] * -1; //Invert the Transverse B field
                _B[-g,j].x = _B[g,j].x; //Copy the Normal B field  (reflected over 0 instead of -1/2)
            }
        }
    }
    if (faces & X_positive) {
        for (int j = j0 ; j <= jn; j++) {
            _B[nx,j].y = _B[nx-1,j].y * -1; //Invert the Transverse B field in the first ghost
            _B[nx,j].z = _B[nx-1,j].z * -1; //Invert the Transverse B field in the first ghost
            for (int g = 1; g <= ng; g++) {
                _B[nx+g,j] = _B[nx-g-1,j] * -1; //Invert the Transverse B field
                _B[nx+g,j].x = _B[nx-g,j].x; //Copy the Normal A field
            }
        }
    }
    if (faces & Y_negative) {
        for (int i = i0 ; i <= in; i++) {
            for (int g = 1; g <= ng; g++) {
                _B[i,-g] = _B[i,g-1] * -1; //Invert the Transverse B field
                _B[i,-g].y = _B[i,g].y; //Copy the Normal B field (reflected over 0 instead of -1/2)
            }
        }
    }
    if (faces & Y_positive) {
        for (int i = i0 ; i <= in; i++) {
            _B[i,ny].x = _B[i,ny-1].x * -1; //Invert the Transverse B field in the first ghost
            _B[i,ny].z = _B[i,ny-1].z * -1; //Invert the Transverse B field in the first ghost
            for (int g = 1; g <= ng; g++) {
                _B[i,ny+g] = _B[i,ny-g-1] * -1; //Invert the Transverse B field
                _B[i,ny+g].y = _B[i,ny-g].y; //Copy the Normal B field
            }
        }
    }
    #endif

}
//MARK: 3D
//Set each ghost to be a reflection of the cell which is opposite the boundary from it
void Boundary::Reflective::apply(Grid3D& grid){
    int ng = grid.getGhosts(), nx = grid.getSizeX(), ny = grid.getSizeY(), nz = grid.getSizeZ();
    int i0 = (corners ? -ng : 0), in = (corners ? nx + ng : nx);
    int j0 = i0, jn = (corners ? ny + ng : ny);
    int k0 = i0, kn = (corners ? nz + ng : nz);

    auto& q = grid.passives();

    if (faces & X_negative) {
        for (int j = j0 ; j < jn; j++) {
            for (int k = k0 ; k < kn; k++) {
                for (int g = 1; g <= ng; g++) {
                    grid[-g,j,k] = grid[g-1,j,k];
                    grid[-g,j,k].v.x *= -1;//Mirror Normal Velocity
                    #ifdef MHD
                    grid[-g,j,k].B *= -1;
                    grid[-g,j,k].B.x *= -1;
                    #endif
                }
                if (j >= -1 && j <= ny && k >= -1 && k <= nz)
                    q[-1,j,k] = q[0,j,k]; //Copy passives (only first ghost)
            }
        }
    }
    if (faces & X_positive) {
        for (int j = j0 ; j < jn; j++) {
            for (int k = k0 ; k < kn; k++) {
                for (int g = 1; g <= ng; g++) {
                    grid[nx-1+g,j,k] = grid[nx-g,j,k];
                    grid[nx-1+g,j,k].v.x *= -1;//Mirror Normal Velocity
                    #ifdef MHD
                    grid[nx-1+g,j,k].B *= -1;
                    grid[nx-1+g,j,k].B.x *= -1;
                    #endif
                }
                if (j >= -1 && j <= ny && k >= -1 && k <= nz)
                    q[nx,j,k] = q[nx-1,j,k]; //Copy passives (only first ghost)
            }
        }
    }
    if (faces & Y_negative) {
        for (int i = i0 ; i < in; i++) {
            for (int k = k0 ; k < kn; k++) {
                for (int g = 1; g <= ng; g++) {
                    grid[i,-g,k] = grid[i,g-1,k];
                    grid[i,-g,k].v.y *= -1;//Mirror Normal Velocity
                    #ifdef MHD
                    grid[i,-g,k].B *= -1;
                    grid[i,-g,k].B.y *= -1;
                    #endif
                }
                if (i >= -1 && i <= nx && k >= -1 && k <= nz)
                    q[i,-1,k] = q[i,0,k]; //Copy passives (only first ghost)
            }
        }
    }
    if (faces & Y_positive) {
        for (int i = i0 ; i < in; i++) {
            for (int k = k0 ; k < kn; k++) {
                for (int g = 1; g <= ng; g++) {
                    grid[i,ny-1+g,k] = grid[i,ny-g,k];
                    grid[i,ny-1+g,k].v.y *= -1;//Mirror Normal Velocity
                    #ifdef MHD
                    grid[i,ny-1+g,k].B *= -1;
                    grid[i,ny-1+g,k].B.y *= -1;
                    #endif
                }
                if (i >= -1 && i <= nx && k >= -1 && k <= nz)
                    q[i,ny,k] = q[i,ny-1,k]; //Copy passives (only first ghost)
            }
        }
    }
    if (faces & Z_negative) {
        for (int i = i0 ; i < in; i++) {
            for (int j = j0 ; j < jn; j++) {
                for (int g = 1; g <= ng; g++) {
                    grid[i,j,-g] = grid[i,j,g-1];
                    grid[i,j,-g].v.z *= -1;//Mirror Normal Velocity
                    #ifdef MHD
                    grid[i,j,-g].B *= -1;
                    grid[i,j,-g].B.z *= -1;
                    #endif
                }
                if (i >= -1 && i <= nx && j >= -1 && j <= ny)
                    q[i,j,-1] = q[i,j,0]; //Copy passives (only first ghost)
            }
        }
    }
    if (faces & Z_positive) {
        for (int i = i0 ; i < in; i++) {
            for (int j = j0 ; j < jn; j++) {
                for (int g = 1; g <= ng; g++) {
                    grid[i,j,nz-1+g] = grid[i,j,nz-g];
                    grid[i,j,nz-1+g].v.z *= -1;//Mirror Normal Velocity
                    #ifdef MHD
                    grid[i,j,nz-1+g].B *= -1;
                    grid[i,j,nz-1+g].B.z *= -1;
                    #endif
                }
                if (i >= -1 && i <= nx && j >= -1 && j <= ny)
                    q[i,j,nz] = q[i,j,nz-1]; //Copy passives (only first ghost)
            }
        }
    }
//MARK: 3D MHD
    #ifdef MHD //Mirror Transverse Magnetic Fields, preserve normal magnetic fields
    auto& _B = grid._B();// B lives on faces, which have one more physical point per dimension than w.
    if (faces & X_negative) {
        for (int j = j0 ; j <= jn; j++) {
            for (int k = k0 ; k <= kn; k++) {
                for (int g = 1; g <= ng; g++) {
                    _B[-g,j,k] = _B[g-1,j,k] * -1; //Invert the Transverse B field
                    _B[-g,j,k].x = _B[g,j,k].x; //Copy the Normal B field  (reflected over 0 instead of -1/2)
                }
            }
        }
    }
    if (faces & X_positive) {
        for (int j = j0 ; j <= jn; j++) {
            for (int k = k0 ; k <= kn; k++) {
                _B[nx,j,k].y = _B[nx-1,j,k].y * -1; //Invert the Transverse B field in the first ghost
                _B[nx,j,k].z = _B[nx-1,j,k].z * -1; //Invert the Transverse B field in the first ghost
                for (int g = 1; g <= ng; g++) {
                    _B[nx+g,j,k] = _B[nx-g-1,j,k] * -1; //Invert the Transverse B field
                    _B[nx+g,j,k].x = _B[nx-g,j,k].x; //Copy the Normal A field
                }
            }
        }
    }
    if (faces & Y_negative) {
        for (int i = i0 ; i <= in; i++) {
            for (int k = k0 ; k <= kn; k++) {
                for (int g = 1; g <= ng; g++) {
                    _B[i,-g,k] = _B[i,g-1,k] * -1; //Invert the Transverse B field
                    _B[i,-g,k].y = _B[i,g,k].y; //Copy the Normal B field (reflected over 0 instead of -1/2)
                }
            }
        }
    }
    if (faces & Y_positive) {
        for (int i = i0 ; i <= in; i++) {
            for (int k = k0 ; k <= kn; k++) {
                _B[i,ny,k].x = _B[i,ny-1,k].x * -1; //Invert the Transverse B field in the first ghost
                _B[i,ny,k].z = _B[i,ny-1,k].z * -1; //Invert the Transverse B field in the first ghost
                for (int g = 1; g <= ng; g++) {
                    _B[i,ny+g,k] = _B[i,ny-g-1,k] * -1; //Invert the Transverse B field
                    _B[i,ny+g,k].y = _B[i,ny-g,k].y; //Copy the Normal B field
                }
            }
        }
    }
    if (faces & Z_negative) {
        for (int i = i0 ; i <= in; i++) {
            for (int j = j0 ; j <= jn; j++) {
                for (int g = 1; g <= ng; g++) {
                    _B[i,j,-g] = _B[i,j,g-1] * -1; //Invert the Transverse B field
                    _B[i,j,-g].z = _B[i,j,g].z; //Copy the Normal B field (reflected over 0 instead of -1/2)
                }
            }
        }
    }
    if (faces & Z_positive) {
        for (int i = i0 ; i <= in; i++) {
            for (int j = j0 ; j <= jn; j++) {
                _B[i,j,nz].x = _B[i,j,nz-1].x * -1; //Invert the Transverse B field in the first ghost
                _B[i,j,nz].y = _B[i,j,nz-1].y * -1; //Invert the Transverse B field in the first ghost
                for (int g = 1; g <= ng; g++) {
                    _B[i,j,nz+g] = _B[i,j,nz-g-1] * -1; //Invert the Transverse B field
                    _B[i,j,nz+g].z = _B[i,j,nz-g].z; //Copy the Normal B field
                }
            }
        }
    }
    #endif
}
