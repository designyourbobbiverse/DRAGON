//
//  Testing.cpp
//  DRAGON/Testing
//
//  Created by Bobbie Markwick on 12/05/2026.
//
#include "Testing.hpp"
#include <iostream>
#include "GridLink.hpp"

using namespace DRAGON_Test;

/*
void test_gridlink_2D_X_positive(){
    Grid2D A(4, 4, 1.0, 1.0, 2);
    Grid2D B(4, 4, 1.0, 1.0, 2);

    for(int j = 0; j < 4; j++){
        for(int i = 0; i < 4; i++){
            A[i,j] = make_state(100 + 10*i + j, 0, 0, 0, 1);
            B[i,j] = make_state(200 + 10*i + j, 0, 0, 0, 1);
        }
    }

    A.boundary = GridLink(&B, Boundary::X_positive);
    A.boundary.apply(A);

    assert((A[4,0].rho == B[0,0].rho));
    assert((A[5,0].rho == B[1,0].rho));
    assert((A[4,3].rho == B[0,3].rho));
    assert((A[5,3].rho == B[1,3].rho));
}

void AMRGrid2D::debugCheckChildGhosts(){
    if(ncx * ncy == 1) return;

    // Make sure child ghost cells are actually filled.
    for(auto& child : children){
        child->data.boundary.apply(child->data);
    }

    // X seams
    for(int zi = 0; zi < ncx - 1; ++zi){
        for(int zj = 0; zj < ncy; ++zj){
            AMRGrid2D& L = *children[ncy * zi + zj];
            AMRGrid2D& R = *children[ncy * (zi + 1) + zj];

            int g  = L.getGhosts();
            int nx = L.getSizeX();
            int ny = L.getSizeY();

            for(int q = 0; q < g; ++q){
                for(int j = 0; j < ny; ++j){
                    assert((L[nx + q, j] == R[q, j]));
                    assert((R[-1 - q, j] == L[nx - 1 - q, j]));
                }
            }
        }
    }

    // Y seams
    for(int zi = 0; zi < ncx; ++zi){
        for(int zj = 0; zj < ncy - 1; ++zj){
            AMRGrid2D& B = *children[ncy * zi + zj];
            AMRGrid2D& T = *children[ncy * zi + (zj + 1)];

            int g  = B.getGhosts();
            int nx = B.getSizeX();
            int ny = B.getSizeY();

            for(int q = 0; q < g; ++q){
                for(int i = 0; i < nx; ++i){
                    assert((B[i, ny + q] == T[i, q]));
                    assert((T[i, -1 - q] == B[i, ny - 1 - q]));
                }
            }
        }
    }
}
*/

int main() {
    
    verify_fluid_element();
    verify_riemann();
    verify_grid();
    verify_boundary();
    verify_cfl();
    
    verify_godunov_1D();
    verify_godunov_2D_Split();
    verify_godunov_2D_Unsplit();
    verify_godunov_3D_Split();
    verify_godunov_3D_Unsplit();
    std::cout << "All Godunov Scheme tests passed.\n\n";

    
    std::cout << "All tests passed.\n";

    
    
}
