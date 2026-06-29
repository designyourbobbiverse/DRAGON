//
//  BoundaryTypes/Jet.cpp
//  DRAGON/Boundary
//
//  Created by Bobbie Markwick on 28/06/2026.
//

#include "Boundary.hpp"
#include "Grid.hpp"
#include <cmath>

using namespace Boundary;
#ifdef MHD
//MARK: Constructors
Boundary::ToroidalJet::ToroidalJet(double rho, double v, double p, double beta, double rm, double rj, int face): Jet(rho,v,p,rj,face), beta(beta), rm(rm) { }
Boundary::ToroidalJet::ToroidalJet(double rho, double v, double p, double beta, double rm, double rj, std::string face): Jet(rho,v,p,rj,face), beta(beta), rm(rm){ }

//MARK: 1D
//All applicable ghost cells are set to equal this->state
void Boundary::ToroidalJet::apply(Grid1D& grid) {
    int ng = grid.getGhosts();
    
    PrimitiveState state;
    state.rho = rho;
    //pm + pt = pa
    //pt/pm = beta
    //pm + pm*beta = pa
    state.p = p;
    state.B.y = p * 8*M_PI/(1.0+ beta);
    
    if (jetface & X_negative){
        state.v = {v,0,0};
        for(int g = 1; g <= ng; g++){
            grid[-g] = state;
            grid[-g].v.x = v;
            grid[-g].p = p;
        }
    }
    if (jetface & X_positive){
        state.v = {-v,0,0};
        int nx = grid.getSize();
        for(int g = 1; g <= ng; g++){
            grid[nx-1+g] = state;
        }
    }
}

//MARK: 2D
//All applicable ghost cells are set to equal this->state
void Boundary::Jet::apply(Grid2D& grid) {
    //Calculate the bounds ahead of time
    double dx = grid.dx, dy = grid.dy;
    int ng = grid.getGhosts(), nx = grid.getSizeX(), ny = grid.getSizeY();
    int i0 = floor(nx*0.5 - rj/dx), in = ceil(nx*0.5 + rj/dx);
    int j0 = floor(ny*0.5 - rj/dy), jn = ceil(ny*0.5 + rj/dy);

    PrimitiveState state;
    state.rho = rho;
    state.p = p;

    if (jetface & X_negative){
        state.v = {v,0,0};
        for(int j = j0 ; j < jn; j++){
            for(int g = 1; g <= ng; g++){
                grid[-g,j] = state;
            }
        }
    }
    if (jetface & X_positive){
        state.v = {-v,0,0};
        for(int j = j0 ; j < jn; j++){
            for(int g = 1; g <= ng; g++){
                grid[nx-1+g,j] = state;
            }
        }
    }
    if (jetface & Y_negative){
        state.v = {0,v,0};
        for(int i = i0 ; i < in; i++){
            for(int g = 1; g <= ng; g++){
                grid[i,-g] = state;
            }
        }
    }
    if (jetface & Y_positive){
        state.v = {0,-v,0};
        for(int i = i0 ; i < in; i++){
            for(int g = 1; g <= ng; g++){
                grid[i,ny-1+g] = state;
            }
        }
    }
//MARK: 2D MHD
    //TODO: Implement
}
//MARK: 3D
//All applicable ghost cells are set to equal this->state
void Boundary::Jet::apply(Grid3D& grid) {
    //Calculate the bounds ahead of time
    double dx = grid.dx, dy = grid.dy, dz = grid.dz;
    int ng = grid.getGhosts(), nx = grid.getSizeX(), ny = grid.getSizeY(), nz = grid.getSizeZ();
    int i0 = floor(nx*0.5 - rj/dx), in = ceil(nx*0.5 + rj/dx);
    int j0 = floor(ny*0.5 - rj/dy), jn = ceil(ny*0.5 + rj/dy);
    int k0 = floor(nz*0.5 - rj/dz), kn = ceil(nz*0.5 + rj/dz);
    
    PrimitiveState state;
    state.rho = rho;
    state.p = p;
    
    if (jetface & X_negative){
        state.v = {v,0,0};
        for(int j = j0; j < jn; j++){
            for(int k = k0; k < kn; k++){
                double y = (j+0.5-ny*0.5)*dy, z = (k+0.5-nz*0.5)*dz;
                if(sqrt(y*y+z*z)>rj) continue;
                
                for(int g = 1; g <= ng; g++){
                    grid[-g,j,k] = state;
                }
            }
        }
    }
    if (jetface & X_positive){
        state.v = {-v,0,0};
        for(int j = j0; j < jn; j++){
            for(int k = k0; k < kn; k++){
                double y = (j+0.5-ny*0.5)*dy, z = (k+0.5-nz*0.5)*dz;
                if(sqrt(y*y+z*z)>rj) continue;
                
                for(int g = 1; g <= ng; g++){
                    grid[nx-1+g,j,k] = state;
                }
            }
        }
    }
    if (jetface & Y_negative){
        state.v = {0,v,0};
        for(int i = i0; i < in; i++){
            for(int k = k0; k < kn; k++){
                double x = (i+0.5-nx*0.5)*dx, z = (k+0.5-nz*0.5)*dz;
                if(sqrt(x*x+z*z)>rj) continue;
                for(int g = 1; g <= ng; g++){
                    grid[i,-g,k] = state;
                }
            }
        }
    }
    if (jetface & Y_positive){
        state.v = {0,-v,0};
        for(int i = i0; i < in; i++){
            for(int k = k0; k < kn; k++){
                double x = (i+0.5-nx*0.5)*dx, z = (k+0.5-nz*0.5)*dz;
                if(sqrt(x*x+z*z)>rj) continue;
                
                for(int g = 1; g <= ng; g++){
                    grid[i,ny-1+g,k] = state;
                }
            }
        }
    }
    if (jetface & Z_negative){
        state.v = {0,0,v};
        for(int i = i0; i < in; i++){
            for(int j = j0; j < jn; j++){
                double x = (i+0.5-nx*0.5)*dx, y = (j+0.5-ny*0.5)*dy;
                if(sqrt(x*x+y*y)>rj) continue;
                
                for(int g = 1; g <= ng; g++){
                    grid[i,j,-g] = state;
                }
            }
        }
    }
    if (jetface & Z_positive){
        state.v = {0,0,-v};
        for(int i = i0; i < in; i++){
            for(int j = j0; j < jn; j++){
                double x = (i+0.5-nx*0.5)*dx, y = (j+0.5-ny*0.5)*dy;
                if(sqrt(x*x+y*y)>rj) continue;
                
                for(int g = 1; g <= ng; g++){
                    grid[i,j,nz-1+g] = state;
                }
            }
        }
    }
    //MARK: 3D MHD
    
    //TODO: Implement
    
}
#endif
