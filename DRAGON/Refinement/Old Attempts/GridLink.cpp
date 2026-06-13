//
//  Boundary/GhostFill.hpp
//  DRAGON
//
//  Created by Bobbie Markwick on 28/05/2026.
//

#include "GridLink.hpp"
#include <cassert>

GridLink::GridLink(Grid1D* grid_, int face): GhostFill(face, true), ref1D(grid_) {
    int count = 0;
    if(faces & Boundary::X_negative) count++;
    if(faces & Boundary::X_positive) count++;
    assert(count == 1);
    reload1D();
}
GridLink::GridLink(Grid2D* grid_, int face): GhostFill(face, true), ref2D(grid_) {
    int count = 0;
    if(faces & Boundary::X_negative) count++;
    if(faces & Boundary::X_positive) count++;
    if(faces & Boundary::Y_negative) count++;
    if(faces & Boundary::Y_positive) count++;
    assert(count == 1);
    reload2D();
}
GridLink::GridLink(Grid3D* grid_, int face): GhostFill(face, true), ref3D(grid_) {
    int count = 0;
    if(faces & Boundary::X_negative) count++;
    if(faces & Boundary::X_positive) count++;
    if(faces & Boundary::Y_negative) count++;
    if(faces & Boundary::Y_positive) count++;
    if(faces & Boundary::Z_negative) count++;
    if(faces & Boundary::Z_positive) count++;
    assert(count == 1);
    reload3D();
}


void GridLink::apply(Grid1D& target){
    int nx = target.getSize();
    int g = target.getGhosts();
    if(g != ghosts1D.size()) {
        reload1D();
        if(g > ghosts1D.size()) g = static_cast<int>(ghosts1D.size());
    }
    
    if(faces & Boundary::X_negative) {
        for(int i=0; i < g; i++) target[-1-i] = ghosts1D[i];
    }
    if(faces & Boundary::X_positive) {
        for(int i=0; i < g; i++) target[nx+i] = ghosts1D[i];
    }
}
void GridLink::apply(Grid2D& target){
    int nx = target.getSizeX(), ny = target.getSizeY();
    int g = target.getGhosts();
    if(g != ghosts2D.size()) {
        reload2D();
        if(g > ghosts2D.size()) g =  static_cast<int>(ghosts2D.size());
    }
    if(faces & Boundary::X_negative) {
        for(int i=0; i< g; i++){
            for(int j=-g; j<ny+g; j++){
                target[-1-i,j] = ghosts2D[i][j+g];
            }
        }
    }
    if(faces & Boundary::X_positive) {
        for(int i=0; i<g; i++){
            for(int j=-g; j<ny+g; j++){
                target[nx+i,j] = ghosts2D[i][j+g];
            }
        }
    }
    if(faces & Boundary::Y_negative) {
        for(int j=0; j< g; j++){
            for(int i=-g; i<nx+g; i++){
                target[i,-1-j] = ghosts2D[j][i+g];
            }
        }
           
    }
    if(faces & Boundary::Y_positive) {
        for(int j=0; j< g; j++){
            for(int i=-g; i<nx+g; i++){
                target[i, ny+j] = ghosts2D[j][i+g];
            }
        }
    }
    
    
}
void GridLink::apply(Grid3D& target) {
    int nx = target.getSizeX(), ny = target.getSizeY(), nz = target.getSizeY();
    int g = target.getGhosts();
    if(g != ghosts3D.size()) {
        reload3D();
        if(g > ghosts3D.size()) g = static_cast<int>(ghosts3D.size());
    }
    if(faces & Boundary::X_negative) {
        for(int i=0; i< g; i++){
            for(int j=-g; j<ny+g; j++){
                for(int k=-g; k<nz+g; k++){
                    target[-1-i,j,k] = ghosts3D[i][j+g][k+g];
                }
            }
        }
           
    }
    if(faces & Boundary::X_positive) {
        for(int i=0; i<g; i++){
            for(int j=-g; j<ny+g; j++){
                for(int k=-g; k<nz+g; k++){
                    target[nx+i,j,k] = ghosts3D[i][j+g][k+g];
                }
            }
        }
    }
    if(faces & Boundary::Y_negative) {
        for(int j=0; j< g; j++){
            for(int i=-g; i<nx+g; i++){
                for(int k=-g; k<nz+g; k++){
                    target[i,-1-j,k] = ghosts3D[j][i+g][k+g];
                }
            }
        }
           
    }
    if(faces & Boundary::Y_positive) {
        for(int j=0; j< g; j++){
            for(int i=-g; i<nx+g; i++){
                for(int k=-g; k<nz+g; k++){
                    target[i, ny+j,k] = ghosts3D[j][i+g][k+g];
                }
            }
        }
    }
    if(faces & Boundary::Z_negative) {
        for(int k=0; k< g; k++){
            for(int i=-g; i<nx+g; i++){
                for(int j=-g; j<ny+g; j++){
                    target[i,j,-1-k] = ghosts3D[k][i+g][j+g];
                }
            }
        }
           
    }
    if(faces & Boundary::Y_positive) {
        for(int k=0; k< g; k++){
            for(int i=-g; i<nx+g; i++){
                for(int j=-g; j<ny+g; j++){
                    target[i,j,nz+k] = ghosts3D[k][i+g][j+g];
                }
            }
        }
    }
}


void GridLink::reload1D(){
    Grid1D& ref = *ref1D;
    int g = ref.getGhosts(), nx = ref.getSize();
    ghosts1D.clear();ghosts2D.reserve(g);
    
    if(faces & Boundary::X_negative) {
        for(int i=1; i<= g; i++) {
            ghosts1D.push_back( ref[nx-i]);
        }
    } else if(faces & Boundary::X_positive) {
        for(int i=0; i< g; i++) {
            ghosts1D.push_back( ref[i]);
        }
    }
}
void GridLink::reload2D(){
    Grid2D& ref = *ref2D;
    int g = ref.getGhosts(), nx = ref.getSizeX(), ny = ref.getSizeY();
    ghosts2D.clear(); ghosts2D.reserve(g);

    
    if(faces & Boundary::X_negative) {
        for(int i=1; i<= g; i++){
            std::vector<PrimitiveState> layer;
            layer.reserve(ny + 2*g);
            
            for(int j=-g; j<ny+g; j++){
                layer.push_back(ref[nx-i,j]);
            }
            ghosts2D.push_back(std::move(layer));
        }
           
    } else if(faces & Boundary::X_positive) {
        for(int i=0; i<g; i++){
            std::vector<PrimitiveState> layer;
            layer.reserve(ny + 2*g);
            
            for(int j=-g; j<ny+g; j++){
                layer.push_back(ref[i,j]);
            }
            ghosts2D.push_back(std::move(layer));
        }
    } else if(faces & Boundary::Y_negative) {
        for(int j=1; j<=g; j++){
            std::vector<PrimitiveState> layer;
            layer.reserve(nx + 2*g);
            
            for(int i=-g; i<nx+g; i++){
                layer.push_back(ref[i,ny-j]);
            }
            ghosts2D.push_back(std::move(layer));
        }
           
    } else if(faces & Boundary::Y_positive) {
        for(int j=0; j<g; j++){
            std::vector<PrimitiveState> layer;
            layer.reserve(nx + 2*g);
            
            for(int i=-g; i<nx+g; i++){
                layer.push_back(ref[i,j]);
            }
            ghosts2D.push_back(std::move(layer));
        }
    }
}
void GridLink::reload3D(){
    Grid3D& ref = *ref3D;
    int g = ref.getGhosts(), nx = ref.getSizeX(), ny = ref.getSizeY(), nz = ref.getSizeZ();
    ghosts3D.clear();
    if(faces & Boundary::X_negative) {
        for(int i=1; i<= g; i++){
            std::vector<std::vector<PrimitiveState>> layer;
            layer.reserve(ny + 2*g);
            for(int j=-g; j<ny+g; j++){
                std::vector<PrimitiveState> row;
                row.reserve(nz + 2*g);
                for(int k=-g; k<nz+g; k++){
                    row.push_back(ref[nx-i,j,k]);
                }
                layer.push_back(std::move(row));
            }
            ghosts3D.push_back(std::move(layer));
        }
    } else if(faces & Boundary::X_positive) {
        for(int i=0; i< g; i++){
            std::vector<std::vector<PrimitiveState>> layer;
            layer.reserve(ny + 2*g);
            for(int j=-g; j<ny+g; j++){
                std::vector<PrimitiveState> row;
                row.reserve(nz + 2*g);
                for(int k=-g; k<nz+g; k++){
                    row.push_back(ref[i,j,k]);
                }
                layer.push_back(std::move(row));
            }
            ghosts3D.push_back(std::move(layer));
        }
    } else if(faces & Boundary::Y_negative) {
        for(int j=1; j<= g; j++){
            std::vector<std::vector<PrimitiveState>> layer;
            layer.reserve(nx + 2*g);
            for(int i=-g; i<nx+g; i++){
                std::vector<PrimitiveState> row;
                row.reserve(nz + 2*g);
                for(int k=-g; k<nz+g; k++){
                    row.push_back(ref[i,ny-j,k]);
                }
                layer.push_back(std::move(row));
            }
            ghosts3D.push_back(std::move(layer));
        }
    } else if(faces & Boundary::Y_positive) {
        for(int j=0; j< g; j++){
            std::vector<std::vector<PrimitiveState>> layer;
            layer.reserve(nx + 2*g);
            for(int i=-g; i<nx+g; i++){
                std::vector<PrimitiveState> row;
                row.reserve(nz + 2*g);
                for(int k=-g; k<nz+g; k++){
                    row.push_back(ref[i,j,k]);
                }
                layer.push_back(std::move(row));
            }
            ghosts3D.push_back(std::move(layer));
        }
    } else if(faces & Boundary::Z_negative) {
        for(int k=1; k <= g; k++){
            std::vector<std::vector<PrimitiveState>> layer;
            layer.reserve(nx + 2*g);
            for(int i=-g; i<nx+g; i++){
                std::vector<PrimitiveState> row;
                row.reserve(ny + 2*g);
                for(int j=-g; j<ny+g; j++){
                    row.push_back(ref[i,j,nz-k]);
                }
                layer.push_back(std::move(row));
            }
            ghosts3D.push_back(std::move(layer));
        }
    } else if(faces & Boundary::Z_positive) {
        for(int k=0; k < g; k++){
            std::vector<std::vector<PrimitiveState>> layer;
            layer.reserve(nx + 2*g);
            for(int i=-g; i<nx+g; i++){
                std::vector<PrimitiveState> row;
                row.reserve(ny + 2*g);
                for(int j=-g; j<ny+g; j++){
                    row.push_back(ref[i,j,k]);
                }
                layer.push_back(std::move(row));
            }
            ghosts3D.push_back(std::move(layer));
        }
    }
    

}
