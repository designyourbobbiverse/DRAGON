//
//  Faraday_Tests.cpp
//  DRAGON/Testing/MHD
//
//  This file contains all the steps needed to go from E to an updated B
//  Created by Bobbie Markwick on 19/07/2026.
//

#include "Testing.hpp"
#include "CT.hpp"
#include "Constants.h"

#ifdef MHD
using namespace DRAGON_Test;


//MARK: E = -dA/dt
void DRAGON_Test::verify_ct_E_updates_A_2D(){
    const int nx = 10, ny = 10, ng = 2;
    MagneticArray2D A0(nx,ny,ng);

    for(int i = -ng; i < nx+ng; i++){
        for(int j = -ng; j < ny+ng; j++){
            A0[i,j] = { (rand()%2000001)*1e-3 - 1e3, (rand()%2000001)*1e-3 - 1e3 , (rand()%2000001)*1e-3 - 1e3};
        }
    }
    
    MagneticArray2D A(nx,ny,ng);
    A.clone(A0);
    
    MagneticArray2D E(nx,ny,ng);
    for(int i = -ng; i < nx+ng; i++){
        for(int j = -ng; j < ny+ng; j++){
            E[i,j] = { (rand()%2000001)*1e-3 - 1e3, (rand()%2000001)*1e-3 - 1e3 , (rand()%2000001)*1e-3 - 1e3};
        }
    }
    //Check a single update
    CT::updatePotential(A, E, 1.0,ng);
    for(int i = -ng; i < nx+ng; i++){
        for(int j = -ng; j < ny+ng; j++){
            expect_close(A[i,j], A0[i,j] - E[i,j]);
        }
    }
    //Check another update to catch dt bugs
    CT::updatePotential(A, E, -2.0,ng);
    for(int i = -ng; i < nx+ng; i++){
        for(int j = -ng; j < ny+ng; j++){
            expect_close(A[i,j], A0[i,j] + E[i,j]);
        }
    }
}
void DRAGON_Test::verify_ct_E_updates_A_3D(){
    const int nx = 10, ny = 10, nz = 10, ng = 2;
    MagneticArray3D A0(nx,ny,nz,ng);

    for(int i = -ng; i < nx+ng; i++){
        for(int j = -ng; j < ny+ng; j++){
            for(int k = -ng; k < nz+ng; k++){
                A0[i,j,k] = { (rand()%2000001)*1e-3 - 1e3, (rand()%2000001)*1e-3 - 1e3 , (rand()%2000001)*1e-3 - 1e3};
            }
        }
    }
    
    MagneticArray3D A(nx,ny,nz,ng);
    A.clone(A0);
    
    MagneticArray3D E(nx,ny,nz,ng);
    for(int i = -ng; i < nx+ng; i++){
        for(int j = -ng; j < ny+ng; j++){
            for(int k = -ng; k < nz+ng; k++){
                E[i,j,k] = { (rand()%2000001)*1e-3 - 1e3, (rand()%2000001)*1e-3 - 1e3 , (rand()%2000001)*1e-3 - 1e3};
            }
        }
    }
    //Check a single update
    CT::updatePotential(A, E, 1.0,ng);
    for(int i = -ng; i < nx+ng; i++){
        for(int j = -ng; j < ny+ng; j++){
            for(int k = -ng; k < nz+ng; k++){
                expect_close(A[i,j,k], A0[i,j,k] - E[i,j,k]);
            }
        }
    }
    //Check another update to catch dt bugs
    CT::updatePotential(A, E, -2.0,ng);
    for(int i = -ng; i < nx+ng; i++){
        for(int j = -ng; j < ny+ng; j++){
            for(int k = -ng; k < nz+ng; k++){
                expect_close(A[i,j,k], A0[i,j,k] + E[i,j,k]);
            }
        }
    }
}


//MARK: Uniform E
void DRAGON_Test::verify_ct_uniform_E_2D(){
    const int nx = 10, ny = 10, ng = 2;
    MagneticArray2D A(nx,ny,ng);

    for(int i = -ng; i < nx+ng; i++){
        for(int j = -ng; j < ny+ng; j++){
            A[i,j] = { (rand()%2000001)*1e-3 - 1e3, (rand()%2000001)*1e-3 - 1e3 , (rand()%2000001)*1e-3 - 1e3};
        }
    }
    MagneticArray2D B(nx,ny,ng);
    CT::computeFaceFields(A, B, 1, 1);
    
    
    MagneticArray2D expected(nx,ny,ng);
    expected.clone(B);
    
    MagneticArray2D E(nx,ny,ng);
    for(int i = -ng; i < nx+ng; i++){
        for(int j = -ng; j < ny+ng; j++){
            E[i,j] = {1,2,3};
        }
    }
    CT::updatePotential(A, E, 1.0,ng);
    CT::computeFaceFields(A, B, 1, 1);
    for(int i = -ng; i < nx+ng; i++){
        for(int j = -ng; j < ny+ng; j++){
            expect_close(B[i,j], expected[i,j]);
        }
    }
}
void DRAGON_Test::verify_ct_uniform_E_3D(){
    const int nx = 10, ny = 10, nz = 10, ng = 2;
    MagneticArray3D A(nx,ny,nz,ng);

    for(int i = -ng; i < nx+ng; i++){
        for(int j = -ng; j < ny+ng; j++){
            for(int k = -ng; k < nz+ng; k++){
                A[i,j,k] = { (rand()%2000001)*1e-3 - 1e3, (rand()%2000001)*1e-3 - 1e3 , (rand()%2000001)*1e-3 - 1e3};
            }
        }
    }
    MagneticArray3D B(nx,ny,nz,ng);
    CT::computeFaceFields(A, B, 1, 1,1);
    
    
    MagneticArray3D expected(nx,ny,nz,ng);
    expected.clone(B);
    
    MagneticArray3D E(nx,ny,nz,ng);
    for(int i = -ng; i < nx+ng; i++){
        for(int j = -ng; j < ny+ng; j++){
            for(int k = -ng; k < nz+ng; k++){
                E[i,j,k] = {1,2,3};
            }
        }
    }
    CT::updatePotential(A, E, 1.0,ng);
    CT::computeFaceFields(A, B, 1, 1,1);
    for(int i = -ng; i < nx+ng; i++){
        for(int j = -ng; j < ny+ng; j++){
            for(int k = -ng; k < nz+ng; k++){
                expect_close(B[i,j,k], expected[i,j,k]);
            }
        }
    }
    

}

//Tests above this line, I wrote myself from scratch
//Tests below this line, AI originally drafted.  I have since revised the tests to better match my coding style.

//MARK: B = curl(A)
//Theses tests check that Stokes Theorem Holds for B = curl(A)
void DRAGON_Test::verify_ct_compute_faces_2D(){
    const int nx = 3, ny = 4, ng = 1;
    const double dx = 0.5, dy = 2.0;
    MagneticArray2D A(nx+1, ny+1, ng);
    for(int i = -ng; i < nx+1+ng; i++){
        for(int j = -ng; j < ny+1+ng; j++){
            // Ax = 2i+j, Ay = i-3j, Az = i²+3ij+2j²
            A[i,j] = vec3{2.0*i + j, 1.0*i - 3.0*j, 1.0*i*i + 3.0*i*j + 2.0*j*j};
        }
    }

    MagneticArray2D B(nx+1, ny+1, ng);
    CT::computeFaceFields(A, B, dx, dy);

    for(int i = -ng; i <= nx+ng; i++){
        for(int j = -ng; j < ny+ng; j++){
            double line_integral = A[i,j+1].z - A[i,j].z;
            double flux = B[i,j].x * dy;
            assert(approx(flux, line_integral));
        }
    }
    for(int i = -ng; i < nx+ng; i++){
        for(int j = -ng; j <= ny+ng; j++){
            double line_integral = A[i,j].z - A[i+1,j].z;
            double flux = B[i,j].y * dx;
            assert(approx(flux, line_integral));
        }
    }
    for(int i = -ng; i < nx+ng; i++){
        for(int j = -ng; j < ny+ng; j++){
            double line_integral = A[i,j].x*dx + A[i+1,j].y*dy - A[i,j+1].x*dx - A[i,j].y*dy;
            double flux = B[i,j].z * dx * dy;
            assert(approx(flux, line_integral));
        }
    }
}

void DRAGON_Test::verify_ct_compute_faces_3D(){
    const int nx = 2, ny = 3, nz = 4, ng = 1;
    const double dx = 0.5, dy = 2.0, dz = 0.25;
    MagneticArray3D A(nx+1, ny+1, nz+1, ng);
    for(int i = -ng; i < nx+1+ng; i++){
        for(int j = -ng; j < ny+1+ng; j++){
            for(int k = -ng; k < nz+1+ng; k++){
                // Ax = ij+2k², Ay = jk+3i², Az = ki+4j²
                A[i,j,k] = vec3{1.0*i*j + 2.0*k*k, 1.0*j*k + 3.0*i*i, 1.0*k*i + 4.0*j*j};
            }
        }
    }

    MagneticArray3D B(nx+1, ny+1, nz+1, ng);
    CT::computeFaceFields(A, B, dx, dy, dz);

    // _B[i,j,k].x = (_A[i,j+1,k].z - _A[i,j,k].z)*_dy - (_A[i,j,k+1].y - _A[i,j,k].y) * _dz;
    for(int i = -ng; i <= nx+ng; i++){
        for(int j = -ng; j < ny+ng; j++){
            for(int k = -ng; k < nz+ng; k++){
                double line_integral = A[i,j,k].y*dy + A[i,j+1,k].z*dz - A[i,j,k+1].y*dy - A[i,j,k].z*dz;
                double flux = B[i,j,k].x * dy * dz;
                assert(approx(flux, line_integral));
            }
        }
    }
    for(int i = -ng; i < nx+ng; i++){
        for(int j = -ng; j <= ny+ng; j++){
            for(int k = -ng; k < nz+ng; k++){
                double line_integral = A[i,j,k].z*dz + A[i,j,k+1].x*dx - A[i+1,j,k].z*dz - A[i,j,k].x*dx;
                double flux = B[i,j,k].y * dx * dz;
                assert(approx(flux, line_integral));
            }
        }
    }
    for(int i = -ng; i < nx+ng; i++){
        for(int j = -ng; j < ny+ng; j++){
            for(int k = -ng; k <= nz+ng; k++){
                double line_integral = A[i,j,k].x*dx + A[i+1,j,k].y*dy - A[i,j+1,k].x*dx - A[i,j,k].y*dy;
                double flux = B[i,j,k].z * dx * dy;
                assert(approx(flux, line_integral));
            }
        }
    }
}


//MARK: Body-Centred Fields
void DRAGON_Test::verify_ct_body_fields_2D(){
    const int nx = 3, ny = 4, ng = 1;
    const double dx = 0.5, dy = 2.0;
    Grid2D grid(nx,ny,dx,dy,ng);
    for(int i = -ng; i < nx+1+ng; i++){
        for(int j = -ng; j < ny+1+ng; j++){
            double Az = i*i + 3*i*j + 2*j*j;
            grid._A()[i,j] = vec3{7, -4, Az};
        }
    }
    for(int i = 0; i < nx; i++){
        for(int j = 0; j < ny; j++){
            grid[i,j].B = vec3{-7,-8,-9};
        }
    }

    MagneticArray2D B_face(nx+1,ny+1,ng);
    CT::computeFaceFields(grid._A(), B_face, dx, dy);
    grid.initialize_B_fields();

    for(int i = 0; i < nx; i++){
        for(int j = 0; j < ny; j++){
            vec3 expected  = vec3{-7,-8,-9};
            expected.x =  (B_face[i,j].x + B_face[i+1,j].x)/2;
            expected.y = (B_face[i,j].y + B_face[i,j+1].y)/2;
            expected.z = B_face[i,j].z;
            expect_close(grid[i,j].B, expected);
        }
    }
}

void DRAGON_Test::verify_ct_body_fields_3D(){
    const int nx = 2, ny = 3, nz = 4, ng = 1;
    const double dx = 0.5, dy = 2.0, dz = 4.0;
    Grid3D grid(nx,ny,nz,dx,dy,dz,ng);

    for(int i = -ng; i < nx+1+ng; i++){
        for(int j = -ng; j < ny+1+ng; j++){
            for(int k = -ng; k < nz+1+ng; k++){
                grid._A()[i,j,k] = vec3{ i*j + 2.0*k*k, j*k + 3.0*i*i, k*i + 4.0*j*j };
            }
        }
    }

    MagneticArray3D B_face(nx+1,ny+1,nz+1,ng);
    CT::computeFaceFields(grid._A(), B_face, dx, dy, dz);
    grid.initialize_B_fields();

    for(int i = 0; i < nx; i++){
        for(int j = 0; j < ny; j++){
            for(int k = 0; k < nz; k++){
                vec3 expected;
                expected.x = (B_face[i,j,k].x + B_face[i+1,j,k].x)/2;
                expected.y = (B_face[i,j,k].y + B_face[i,j+1,k].y)/2;
                expected.z = (B_face[i,j,k].z + B_face[i,j,k+1].z)/2;
                expect_close(grid[i,j,k].B, expected);
            }
        }
    }
}


//MARK: Face-Centred Field Copy
void DRAGON_Test::verify_ct_copy_face_fields_2D(){
    const int nx = 2, ny = 3, ng = 1;
    MagneticArray2D B(nx+1,ny+1,ng);
    for(int i = -ng; i < nx+1+ng; i++){
        for(int j = -ng; j < ny+1+ng; j++){
            B[i,j] = vec3{100*i + 10*j+1.0, 100*i + 10*j+2.0, 100*i + 10*j+3.0};
        }
    }

    //X test
    FluidArray2D L(nx,ny,ng), R(nx,ny,ng);
    for(int i = -ng; i < nx+ng; i++){
        for(int j = -ng; j < ny+ng; j++){
            L[i,j].B = vec3{-1,-2,-3};
            R[i,j].B = vec3{-4,-5,-6};
        }
    }
    
    CT::copyFaceFields_X(L, B, R);
    
    for(int i = -ng; i < nx+ng; i++){
        for(int j = -ng; j < ny+ng; j++){
            assert((L[i,j].B.x == B[i,j].x));
            assert((R[i,j].B.x == B[i+1,j].x));
            assert((L[i,j].B.y == -2 && L[i,j].B.z == -3));
            assert((R[i,j].B.y == -5 && R[i,j].B.z == -6));
        }
    }

    //Y test
    for(int i=-ng; i<nx+ng; i++){
        for(int j=-ng; j<ny+ng; j++){
            L[i,j].B = vec3{-1,-2,-3};
            R[i,j].B = vec3{-4,-5,-6};
        }
    }
    
    CT::copyFaceFields_Y(L, B, R);
    
    for(int i = -ng; i < nx+ng; i++){
        for(int j = -ng; j < ny+ng; j++){
            assert((L[i,j].B.y == B[i,j].y));
            assert((R[i,j].B.y == B[i,j+1].y));
            assert((L[i,j].B.x == -1 && L[i,j].B.z == -3));
            assert((R[i,j].B.x == -4 && R[i,j].B.z == -6));
        }
    }
}

void DRAGON_Test::verify_ct_copy_face_fields_3D(){
    const int nx = 2, ny = 3, nz = 4, ng = 1;
    MagneticArray3D B(nx+1,ny+1,nz+1,ng);
    for(int i=-ng; i<B.getSizeX()+ng; i++){
        for(int j=-ng; j<B.getSizeY()+ng; j++){
            for(int k=-ng; k<B.getSizeZ()+ng; k++){
                double tag = 10000*i + 100*j + 10*k;
                B[i,j,k] = vec3{tag+1, tag+2, tag+3};
            }
        }
    }

    //X test
    FluidArray3D L(nx,ny,nz,ng), R(nx,ny,nz,ng);
    for(int i = -ng; i < nx+ng; i++){
        for(int j = -ng; j < ny+ng; j++){
            for(int k = -ng; k < nz+ng; k++){
                L[i,j,k].B = vec3{-1,-2,-3};
                R[i,j,k].B = vec3{-4,-5,-6};
            }
        }
    }
    
    CT::copyFaceFields_X(L, B, R);
    
    for(int i = -ng; i < nx+ng; i++){
        for(int j = -ng; j < ny+ng; j++){
            for(int k = -ng; k < nz+ng; k++){
                assert((L[i,j,k].B.x == B[i,j,k].x));
                assert((R[i,j,k].B.x == B[i+1,j,k].x));
                assert((L[i,j,k].B.y == -2 && L[i,j,k].B.z == -3));
                assert((R[i,j,k].B.y == -5 && R[i,j,k].B.z == -6));
            }
        }
    }

    //Y test
    for(int i = -ng; i < nx+ng; i++){
        for(int j = -ng; j < ny+ng; j++){
            for(int k = -ng; k < nz+ng; k++){
                L[i,j,k].B = vec3{-1,-2,-3};
                R[i,j,k].B = vec3{-4,-5,-6};
            }
        }
    }
    
    CT::copyFaceFields_Y(L, B, R);
    
    for(int i = -ng; i < nx+ng; i++){
        for(int j = -ng; j < ny+ng; j++){
            for(int k = -ng; k < nz+ng; k++){
                assert((L[i,j,k].B.y == B[i,j,k].y));
                assert((R[i,j,k].B.y == B[i,j+1,k].y));
                assert((L[i,j,k].B.x == -1 && L[i,j,k].B.z == -3));
                assert((R[i,j,k].B.x == -4 && R[i,j,k].B.z == -6));
            }
        }
    }
    //Z test
    for(int i = -ng; i < nx+ng; i++){
        for(int j = -ng; j < ny+ng; j++){
            for(int k = -ng; k < nz+ng; k++){
                L[i,j,k].B = vec3{-1,-2,-3};
                R[i,j,k].B = vec3{-4,-5,-6};
            }
        }
    }
    
    CT::copyFaceFields_Z(L, B, R);
    
    for(int i = -ng; i < nx+ng; i++){
        for(int j = -ng; j < ny+ng; j++){
            for(int k = -ng; k < nz+ng; k++){
                assert((L[i,j,k].B.z == B[i,j,k].z));
                assert((R[i,j,k].B.z == B[i,j,k+1].z));
                assert((L[i,j,k].B.x == -1 && L[i,j,k].B.y == -2));
                assert((R[i,j,k].B.x == -4 && R[i,j,k].B.y == -5));
            }
        }
    }
}
#endif
