//
//  BoundaryTypes/Outflow.cpp
//  DRAGON/Boundary
//
//  Created by Bobbie Markwick on 28/05/2026.
//

#include "Boundary/Boundary.hpp"

#include "Hydro/Grid.hpp" //For filling ghost cells in the grid

using namespace DRAGON;
using namespace Boundary;

//MARK: Constructors
Boundary::Outflow::Outflow(int faces_, bool corners, bool gated_):  GhostFill(faces_, corners), gated(gated_){}
Boundary::Outflow::Outflow(std::string s, bool corners, bool gated): Outflow(face_mask(s),corners,gated){}
Outflow Boundary::Outflow::Gated(int faces, bool corner_ghosts){ return Outflow(faces, corner_ghosts,true); }
Outflow Boundary::Outflow::Gated(std::string s, bool corners){ return Outflow(s,corners,true);}

//MARK: 1D
//Set each relevant ghost to the nearest physical cell, possibly subject to inflow gating
void Boundary::Outflow::apply(Grid1D& grid){
    int ng = grid.getGhosts();
    auto& q = grid.passives();
    if (faces & X_negative) {
        for (int g = 1; g <= ng; g++) {
            grid[-g] = grid[0];
            if (gated && grid[-g].v.x > 0) grid[-g].v.x = 0; //Kill inflows if applicable
        }
        q[-1] = q[0]; //Copy passives (only first ghost)
    }
    if (faces & X_positive) {
        int nx = grid.getSize();
        for (int g = 1; g <= ng; g++) {
            grid[nx-1+g] = grid[nx-1];
            if (gated && grid[nx-1+g].v.x < 0) grid[nx-1+g].v.x = 0; //Kill inflows if applicable
        }
        q[nx] = q[nx-1]; //Copy passives (only first ghost)
    }
}

//MARK: 2D
//Set each relevant ghost to the nearest physical cell, possibly subject to inflow gating
void Boundary::Outflow::apply(Grid2D& grid){
    //Calculate the bounds ahead of time
    int ng = grid.getGhosts(), nx = grid.getSizeX(), ny = grid.getSizeY();
    int i0 = (corners ? -ng : 0), in = (corners ? nx + ng : nx);
    int j0 = i0, jn = (corners ? ny + ng : ny);
    
    auto& q = grid.passives();

    if (faces & X_negative) {
        for (int j = j0 ; j < jn; j++) {
            for (int g = 1; g <= ng; g++) {
                grid[-g,j] = grid[0,j];
                if (gated && grid[-g,j].v.x > 0) grid[-g,j].v.x = 0; //Kill inflows if applicable
            }
            if (j >= -1 && j <= ny) q[-1,j] = q[0,j]; //Copy passives (only first ghost)
        }
    }
    if (faces & X_positive) {
        for (int j = j0 ; j < jn; j++) {
            for (int g = 1; g <= ng; g++) {
                grid[nx-1+g,j] = grid[nx-1,j];
                if (gated && grid[nx-1+g,j].v.x < 0) grid[nx-1+g,j].v.x = 0; //Kill inflows if applicable
            }
            if (j >= -1 && j <= ny) q[nx,j] = q[nx-1,j]; //Copy passives (only first ghost)
        }
    }
    if (faces & Y_negative) {
        for (int i = i0 ; i < in; i++) {
            for (int g = 1; g <= ng; g++) {
                grid[i,-g] = grid[i,0];
                if (gated && grid[i,-g].v.y > 0) grid[i,-g].v.y = 0; //Kill inflows if applicable
            }
            if (i >= -1 && i <= nx) q[i,-1] = q[i,0]; //Copy passives (only first ghost)
        }
    }
    if (faces & Y_positive) {
        for (int i = i0 ; i < in; i++) {
            for (int g = 1; g <= ng; g++) {
                grid[i,ny-1+g] = grid[i,ny-1];
                if (gated && grid[i,ny-1+g].v.y < 0) grid[i,ny-1+g].v.y = 0; //Kill inflows if applicable
            }
            if (i >= -1 && i <= nx) q[i,ny] = q[i,ny-1]; //Copy passives (only first ghost)
        }
    }
//MARK: 2D MHD
    #ifdef MHD //Copy Transverse Fields
    auto& _B = grid._B(); // B lives on faces, which have one more physical point per dimension than w.
    if (faces & X_negative) {
        for (int j = j0; j <= jn; j++) {
            for (int g = 1; g <= ng; g++) {
                _B[-g,j] = _B[0,j];
                _B[-g,j].x = 2*_B[-g+1,j].x - _B[-g+2,j].x;//Linear extrapolation of normal B to keep divB = 0
            }
        }
    }
    if (faces & X_positive) {
        for (int j = j0; j <= jn; j++) {
            _B[nx,j].y = _B[nx-1,j].y;
            _B[nx,j].z = _B[nx-1,j].z;
            for (int g = 1; g <= ng; g++) {
                _B[nx+g,j] = _B[nx,j];
                _B[nx+g,j].x = 2*_B[nx+g-1,j].x - _B[nx+g-2,j].x;//Linear extrapolation of normal B to keep divB = 0
            }
        }
    }
    if (faces & Y_negative) {
        for (int i = i0; i <= in; i++) {
            for (int g = 1; g <= ng; g++) {
                _B[i,-g] = _B[i,0];
                _B[i,-g].y = 2*_B[i,-g+1].y - _B[i,-g+2].y;//Linear extrapolation of normal B to keep divB = 0
            }
        }
    }
    if (faces & Y_positive) {
        for (int i = i0; i <= in; i++) {
            _B[i,ny].x = _B[i,ny-1].x;
            _B[i,ny].z = _B[i,ny-1].z;
            for (int g = 1; g <= ng; g++) {
                _B[i,ny+g] = _B[i,ny];
                _B[i,ny+g].y = 2*_B[i,ny+g-1].y - _B[i,ny+g-2].y;//Linear extrapolation of normal B to keep divB = 0
            }
        }
    }
    #endif
}
//MARK: 3D
//Set each relevant ghost to the nearest physical cell, possibly subject to inflow gating
void Boundary::Outflow::apply(Grid3D& grid){
    int ng = grid.getGhosts(), nx = grid.getSizeX(), ny = grid.getSizeY(), nz = grid.getSizeZ();
    int i0 = (corners ? -ng : 0), in = (corners ? nx + ng : nx);
    int j0 = i0, jn = (corners ? ny + ng : ny);
    int k0 = i0, kn = (corners ? nz + ng : nz);
    
    auto& q = grid.passives();
    
    if (faces & X_negative) {
        for (int j = j0 ; j < jn; j++) {
            for (int k = k0 ; k < kn; k++) {
                for (int g = 1; g <= ng; g++) {
                    grid[-g,j,k] = grid[0,j,k];
                    if (gated && grid[-g,j,k].v.x > 0) grid[-g,j,k].v.x = 0; //Kill inflows if applicable
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
                    grid[nx-1+g,j,k] = grid[nx-1,j,k];
                    if (gated && grid[nx-1+g,j,k].v.x < 0) grid[nx-1+g,j,k].v.x = 0; //Kill inflows if applicable
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
                    grid[i,-g,k] = grid[i,0,k];
                    if (gated && grid[i,-g,k].v.y > 0) grid[i,-g,k].v.y = 0; //Kill inflows if applicable
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
                    grid[i,ny-1+g,k] = grid[i,ny-1,k];
                    if (gated && grid[i,ny-1+g,k].v.y < 0) grid[i,ny-1+g,k].v.y = 0; //Kill inflows if applicable
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
                    grid[i,j,-g] = grid[i,j,0];
                    if (gated && grid[i,j,-g].v.z > 0) grid[i,j,-g].v.z = 0; //Kill inflows if applicable
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
                    grid[i,j,nz-1+g] = grid[i,j,nz-1];
                    if (gated && grid[i,j,nz-1+g].v.z < 0) grid[i,j,nz-1+g].v.z = 0; //Kill inflows if applicable
                }
                if (i >= -1 && i <= nx && j >= -1 && j <= ny)
                    q[i,j,nz] = q[i,j,nz-1]; //Copy passives (only first ghost)
            }
        }
    }
    
//MARK: 3D MHD
    #ifdef MHD //Copy Transverse Fields
    auto& _B = grid._B(); // B lives on faces, which have one more physical point per dimension than w.
    if (faces & X_negative) {
        for (int j = j0; j <= jn; j++) {
            for (int k = k0; k <= kn; k++) {
                for (int g = 1; g <= ng; g++) {
                    _B[-g,j,k] = _B[0,j,k];
                    _B[-g,j,k].x = 2*_B[-g+1,j,k].x - _B[-g+2,j,k].x; //Linear extrapolation of normal B to keep divB = 0
                }
            }
        }
    }
    if (faces & X_positive) {
        for (int j = j0; j <= jn; j++) {
            for (int k = k0; k <= kn; k++) {
                _B[nx,j,k].y = _B[nx-1,j,k].y;
                _B[nx,j,k].z = _B[nx-1,j,k].z;
                for (int g = 1; g <= ng; g++) {
                    _B[nx+g,j,k] = _B[nx,j,k];
                    _B[nx+g,j,k].x = 2*_B[nx+g-1,j,k].x - _B[nx+g-2,j,k].x;//Linear extrapolation of normal B to keep divB = 0
                }
            }
        }
    }
    if (faces & Y_negative) {
        for (int i = i0; i <= in; i++) {
            for (int k = k0; k <= kn; k++) {
                for (int g = 1; g <= ng; g++) {
                    _B[i,-g,k] = _B[i,0,k];
                    _B[i,-g,k].y = 2*_B[i,-g+1,k].y - _B[i,-g+2,k].y;//Linear extrapolation of normal B to keep divB = 0
                }
            }
        }
    }
    if (faces & Y_positive) {
        for (int i = i0; i <= in; i++) {
            for (int k = k0; k <= kn; k++) {
                _B[i,ny,k].x = _B[i,ny-1,k].x;
                _B[i,ny,k].z = _B[i,ny-1,k].z;
                for (int g = 1; g <= ng; g++) {
                    _B[i,ny+g,k] = _B[i,ny,k];
                    _B[i,ny+g,k].y = 2*_B[i,ny+g-1,k].y - _B[i,ny+g-2,k].y;//Linear extrapolation of normal B to keep divB = 0
                }
            }
        }
    }
    if (faces & Z_negative) {
        for (int i = i0; i <= in; i++) {
            for (int j = j0; j <= jn; j++) {
                for (int g = 1; g <= ng; g++) {
                    _B[i,j,-g] = _B[i,j,0];
                    _B[i,j,-g].z = 2*_B[i,j,-g+1].z - _B[i,j,-g+2].z;//Linear extrapolation of normal B to keep divB = 0
                }
            }
        }
    }
    if (faces & Z_positive) {
        for (int i = i0; i <= in; i++) {
            for (int j = j0; j <= jn; j++) {
                _B[i,j,nz].x = _B[i,j,nz-1].x;
                _B[i,j,nz].y = _B[i,j,nz-1].y;
                for (int g = 1; g <= ng; g++) {
                    _B[i,j,nz+g] = _B[i,j,nz];
                    _B[i,j,nz+g].z = 2*_B[i,j,nz+g-1].z - _B[i,j,nz+g-2].z;//Linear extrapolation of normal B to keep divB = 0
                }
            }
        }
    }
    #endif
}
