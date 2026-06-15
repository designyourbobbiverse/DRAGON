//
//  Boundary/Outflow.cpp
//  DRAGON
//
//  Created by Bobbie Markwick on 28/05/2026.
//

#include "Boundary.hpp"
#include "Grid.hpp"

using namespace Boundary;

//MARK: Constructors
Boundary::Outflow::Outflow(int faces_, bool corners, bool gated_):  GhostFill(faces_, corners), gated(gated_) {}
Boundary::Outflow::Outflow(std::string s, bool corners, bool gated): Outflow(face_mask(s),corners,gated) {}
Outflow Boundary::Outflow::Gated(int faces, bool corner_ghosts){ return Outflow(faces, corner_ghosts,true); }
Outflow Boundary::Outflow::Gated(std::string s, bool corners){ return Outflow(s,corners,true);}

//MARK: 1D
void Boundary::Outflow::apply(Grid1D& grid) {
    int ng = grid.getGhosts();
    if (faces & X_negative){
        for(int g = 1; g <= ng; g++){
            grid[-g] = grid[0];
            if(gated && grid[-g].v.x > 0) grid[-g].v.x = 0;
        }
    }
    if (faces & X_positive){
        int nx = grid.getSize();
        for(int g = 1; g <= ng; g++){
            grid[nx-1+g] = grid[nx-1];
            if(gated && grid[nx-1+g].v.x < 0) grid[nx-1+g].v.x = 0;
        }
    }
}

//MARK: 2D
void Boundary::Outflow::apply(Grid2D& grid) {
    int ng = grid.getGhosts(), nx = grid.getSizeX(), ny = grid.getSizeY();
    int i0 = (corners ? -ng : 0), in = (corners ? nx + ng : nx);
    int j0 = i0, jn = (corners ? ny + ng : ny);

    if (faces & X_negative){
        for(int j = j0 ; j < jn; j++){
            for(int g = 1; g <= ng; g++){
                grid[-g,j] = grid[0,j];
                if(gated && grid[-g,j].v.x > 0) grid[-g,j].v.x = 0;
            }
        }
    }
    if (faces & X_positive){
        for(int j = j0 ; j < jn; j++){
            for(int g = 1; g <= ng; g++){
                grid[nx-1+g,j] = grid[nx-1,j];
                if(gated && grid[nx-1+g,j].v.x < 0) grid[nx-1+g,j].v.x = 0;

            }
        }
    }
    if (faces & Y_negative){
        for(int i = i0 ; i < in; i++){
            for(int g = 1; g <= ng; g++){
                grid[i,-g] = grid[i,0];
                if(gated && grid[i,-g].v.y > 0) grid[i,-g].v.y = 0;
            }
        }
    }
    if (faces & Y_positive){
        for(int i = i0 ; i < in; i++){
            for(int g = 1; g <= ng; g++){
                grid[i,ny-1+g] = grid[i,ny-1];
                if(gated && grid[i,ny-1+g].v.y < 0) grid[i,ny-1+g].v.y = 0;
            }
        }
    }
}
//MARK: 3D
void Boundary::Outflow::apply(Grid3D& grid) {
    int ng = grid.getGhosts(), nx = grid.getSizeX(), ny = grid.getSizeY(), nz = grid.getSizeZ();
    int i0 = (corners ? -ng : 0), in = (corners ? nx + ng : nx);
    int j0 = i0, jn = (corners ? ny + ng : ny);
    int k0 = i0, kn = (corners ? nz + ng : nz);
    
    if (faces & X_negative){
        for(int j = j0 ; j < jn; j++){
            for(int k = k0 ; k < kn; k++){
                for(int g = 1; g <= ng; g++){
                    grid[-g,j,k] = grid[0,j,k];
                    if(gated && grid[-g,j,k].v.x > 0) grid[-g,j,k].v.x = 0;
                }
            }
        }
    }
    if (faces & X_positive){
        for(int j = j0 ; j < jn; j++){
            for(int k = k0 ; k < kn; k++){
                for(int g = 1; g <= ng; g++){
                    grid[nx-1+g,j,k] = grid[nx-1,j,k];
                    if(gated && grid[nx-1+g,j,k].v.x < 0) grid[nx-1+g,j,k].v.x = 0;
                }
            }
        }
    }
    if (faces & Y_negative){
        for(int i = i0 ; i < in; i++){
            for(int k = k0 ; k < kn; k++){
                for(int g = 1; g <= ng; g++){
                    grid[i,-g,k] = grid[i,0,k];
                    if(gated && grid[i,-g,k].v.y > 0) grid[i,-g,k].v.y = 0;

                }
            }
        }
    }
    if (faces & Y_positive){
        for(int i = i0 ; i < in; i++){
            for(int k = k0 ; k < kn; k++){
                for(int g = 1; g <= ng; g++){
                    grid[i,ny-1+g,k] = grid[i,ny-1,k];
                    if(gated && grid[i,ny-1+g,k].v.y < 0) grid[i,ny-1+g,k].v.y = 0;
                }
            }
        }
    }
    if (faces & Z_negative){
        for(int i = i0 ; i < in; i++){
            for(int j = j0 ; j < jn; j++){
                for(int g = 1; g <= ng; g++){
                    grid[i,j,-g] = grid[i,j,0];
                    if(gated && grid[i,j,-g].v.z > 0) grid[i,j,-g].v.z = 0;
                }
            }
        }
    }
    if (faces & Z_positive){
        for(int i = i0 ; i < in; i++){
            for(int j = j0 ; j < jn; j++){
                for(int g = 1; g <= ng; g++){
                    grid[i,j,nz-1+g] = grid[i,j,nz-1];
                    if(gated && grid[i,j,nz-1+g].v.z < 0) grid[i,j,nz-1+g].v.z = 0;
                }
            }
        }
    }
}
