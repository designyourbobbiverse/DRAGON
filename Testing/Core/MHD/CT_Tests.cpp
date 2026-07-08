//
//  CT_Tests.cpp
//  DRAGON/Testing/MHD
//
//  Created by Bobbie Markwick on 25/06/2026.
//

#include "Testing.hpp"
#include "CT.hpp"
#include "DistGrid.hpp"
#include "DragonWing.hpp"
#include "Constants.h"
#include <iostream>

#ifdef MHD
using namespace DRAGON_Test;

void DRAGON_Test::verify_ct_2D(bool output){
    if(output) std::cout << "Constrained Transport (2D):\n";
    if(output) std::cout << "- Zero Divergence: ";
    verify_ct_divergence_2D();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "- Internal Function Checks: ";
    verify_ct_copy_face_fields_2D();
    verify_ct_body_fields_2D();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "- Stationary Field: ";
    verify_ct_stationary_2D();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "- Loop Advection: ";
    verify_ct_loop_advection_2D();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "- Alfven Wave: ";
    verify_ct_alfven_wave_2D();
    if(output) std::cout << "Passed\n";
}
void DRAGON_Test::verify_ct_3D(bool output){
    if(output) std::cout << "Constrained Transport (3D): \n";
    if(output) std::cout << "- Zero Divergence: ";
    verify_ct_divergence_3D();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "- Internal Function Checks: ";
    verify_ct_copy_face_fields_3D();
    verify_ct_body_fields_3D();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "- Stationary Field: ";
    verify_ct_stationary_3D();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "- Loop Advection: ";
    verify_ct_loop_advection_3D();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "- Alfven Wave: ";
    verify_ct_alfven_wave_3D();
    if(output) std::cout << "Passed\n";

}




//MARK: Div B = 0
void assert_divergenceless(const MagneticArray2D& A, double dx, double dy){
    const int nx = A.getSizeX(), ny = A.getSizeY(), ng = A.getGhosts();
    MagneticArray2D B(nx,ny,ng);
    CT::computeFaceFields(A, B, dx,dy);
    for(int i = -ng; i < nx+ng - 1 ; i++){
        for(int j = -ng; j < ny+ng - 1; j++){
            double dBx = (B[i+1,j].x - B[i,j].x)/dx;
            double dBy = (B[i,j+1].y - B[i,j].y)/dy;
            assert(fabs(dBx + dBy) < 1e-12);
        }
    }
}
void assert_divergenceless(const MagneticArray3D& A, double dx, double dy, double dz){
    const int nx = A.getSizeX(), ny = A.getSizeY(),nz = A.getSizeZ(), ng = A.getGhosts();
    MagneticArray3D B(nx,ny,nz,ng);
    CT::computeFaceFields(A, B, dx,dy,dz);
    for(int i = -ng; i < nx+ng - 1 ; i++){
        for(int j = -ng; j < ny+ng - 1; j++){
            for(int k = -ng; k < nz+ng - 1; k++){
                double dBx = (B[i+1,j,k].x - B[i,j,k].x)/dx;
                double dBy = (B[i,j+1,k].y - B[i,j,k].y)/dy;
                double dBz = (B[i,j,k+1].z - B[i,j,k].z)/dz;
                assert(fabs(dBx + dBy + dBz) < 1e-12);
            }
        }
    }
}

void DRAGON_Test::verify_ct_divergence_2D(){
    const int nx = 10, ny = 10, ng = 2;
    MagneticArray2D A(nx,ny,ng);
    for(int i = -ng; i < nx+ng; i++){
        for(int j = -ng; j < ny+ng; j++){
            A[i,j] = { (rand()%1000000)*1e-3 - 1e3, (rand()%1000000)*1e-3 - 1e3 , (rand()%1000000)*1e-3 - 1e3};
        }
    }
    assert_divergenceless(A, 1, 2);
}


void DRAGON_Test::verify_ct_divergence_3D(){
    const int nx = 10, ny = 10, nz = 10, ng = 2;
    MagneticArray3D A(nx,ny,nz,ng);
    for(int i = -ng; i < nx+ng; i++){
        for(int j = -ng; j < ny+ng; j++){
            for(int k = -ng; k < nz+ng; k++){
                A[i,j,k] = { (rand()%1000000)*1e-3 - 1e3, (rand()%1000000)*1e-3 - 1e3 , (rand()%1000000)*1e-3 - 1e3};
            }
        }
    }
    
    assert_divergenceless(A, 1, 2,3);
}

//MARK: Stationary Field
void DRAGON_Test::verify_ct_stationary_2D(){
    double dx = M_PI/5;
    Grid2D grid(10,10,dx,dx, 2), expected(10,10,dx,dx, 2);
    double p0 = 5.0;
    PrimitiveState W = make_state(1.0, 0.0, 0.0, 0.0, p0);

    for (int i = 0; i <= grid.getSizeX(); i++){
        for (int j = 0; j <= grid.getSizeY(); j++){
            grid[i,j] = W;
            grid._A()[i,j] = vec3{0,0,-cos(i*dx)-cos(j*dx)};
        }
    }
    grid.boundary = Boundary::Periodic();
    grid.boundary.apply(grid);
    grid.initialize_B_fields();

    for (int i = 0; i <= grid.getSizeX(); i++){
        for (int j = 0; j <= grid.getSizeY(); j++){
            grid[i,j].p = p0 - grid[i,j].B * grid[i,j].B * inv_8pi;
            expected[i,j] = grid[i,j];
        }
    }
    assert_divergenceless(grid._A(),dx,dx);

    
    grid.advance(1.0);
    
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            expect_close(grid[i,j].B, expected[i,j].B);
        }
    }
    assert_divergenceless(grid._A(),dx,dx);
}

void DRAGON_Test::verify_ct_stationary_3D(){
    double dx = M_PI/5;
    Grid3D grid(10,10,10,dx,dx,dx, 2), expected(10,10,10,dx,dx,dx, 2);
    double p0 = 5.0;
    PrimitiveState W = make_state(1.0, 0.0, 0.0, 0.0, p0);

    for (int i = 0; i <= grid.getSizeX(); i++){
        for (int j = 0; j <= grid.getSizeY(); j++){
            for (int k = 0; k <= grid.getSizeZ(); k++){
                grid[i,j,k] = W;
                grid._A()[i,j,k] = vec3{-cos(k*dx)-cos(j*dx),-cos(i*dx)-cos(k*dx),-cos(i*dx)-cos(j*dx)};
            }
        }
    }
    grid.boundary = Boundary::Periodic();
    grid.boundary.apply(grid);
    grid.initialize_B_fields();

    for (int i = 0; i <= grid.getSizeX(); i++){
        for (int j = 0; j <= grid.getSizeY(); j++){
            for (int k = 0; k <= grid.getSizeZ(); k++){
                grid[i,j,k].p = p0 - grid[i,j,k].B * grid[i,j,k].B * inv_8pi;
                expected[i,j,k] = grid[i,j,k];
            }
        }
    }
    assert_divergenceless(grid._A(),dx,dx,dx);

    
    grid.advance(0.2);
    
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            for (int k = 0; k < grid.getSizeZ(); k++){
                expect_close(grid[i,j,k].B, expected[i,j,k].B);
            }
        }
    }
    assert_divergenceless(grid._A(),dx,dx,dx);
}
//MARK: Loop Advection
void DRAGON_Test::verify_ct_loop_advection_2D(){
    double dx = 1.0;
    Grid2D grid(10,10,dx,dx, 2), expected(10,10,dx,dx, 2);
    double p0 = 50.0;
    PrimitiveState W = make_state(1.0, 5.0, 5.0, 0.0, p0);

    for (int i = 0; i <= grid.getSizeX(); i++){
        for (int j = 0; j <= grid.getSizeY(); j++){
            grid[i,j] = W;
            double x = i*dx-5, y = j*dx - 5, r2 = x*x + y*y;
            double Az = fmax(4.0-r2, 0);
            grid._A()[i,j] = vec3{0,0,Az};
        }
    }
    grid.boundary = Boundary::Periodic();
    grid.boundary.apply(grid);
    grid.initialize_B_fields();

    for (int i = 0; i <= grid.getSizeX(); i++){
        for (int j = 0; j <= grid.getSizeY(); j++){
            grid[i,j].p = p0 - grid[i,j].B * grid[i,j].B  * inv_8pi;
            expected[i,j] = grid[i,j];
        }
    }
    assert_divergenceless(grid._A(),dx,dx);

    
    grid.advance(2.0);
    
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            expect_close(grid[i,j].B, expected[i,j].B);
        }
    }
    assert_divergenceless(grid._A(),dx,dx);
}

void DRAGON_Test::verify_ct_loop_advection_3D(){
    double dx = 1.0;
    Grid3D grid(10,10,10,dx,dx,dx, 2), expected(10,10,10,dx,dx,dx, 2);
    double p0 = 50.0;
    PrimitiveState W = make_state(1.0, 5.0, 5.0,5.0, p0);

    for (int i = 0; i <= grid.getSizeX(); i++){
        for (int j = 0; j <= grid.getSizeY(); j++){
            for (int k = 0; k <= grid.getSizeZ(); k++){
                grid[i,j,k] = W;
                double x = i*dx-5, y = j*dx - 5, z = k*dx - 5, r2 = x*x + y*y * z*z;
                double Az = fmax(4.0-r2, 0);
                grid._A()[i,j,k] = vec3{0,0,Az};
            }
        }
    }
    grid.boundary = Boundary::Periodic();
    grid.boundary.apply(grid);
    grid.initialize_B_fields();

    for (int i = 0; i <= grid.getSizeX(); i++){
        for (int j = 0; j <= grid.getSizeY(); j++){
            for (int k = 0; k <= grid.getSizeZ(); k++){
                grid[i,j,k].p = p0 - grid[i,j,k].B * grid[i,j,k].B * inv_8pi;
                expected[i,j,k] = grid[i,j,k];
            }
        }
    }
    assert_divergenceless(grid._A(),dx,dx,dx);

    
    grid.advance(2.0);
    
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            for (int k = 0; k < grid.getSizeZ(); k++){
                expect_close(grid[i,j,k].B, expected[i,j,k].B);
            }
        }
    }
    assert_divergenceless(grid._A(),dx,dx,dx);
}
//MARK: Alfven Wave
void DRAGON_Test::verify_ct_alfven_wave_2D(){
    double dx = 0.1, amp = 0.01;
    Grid2D grid(10,10,dx,dx, 2), expected(10,10,dx,dx, 2);
    PrimitiveState W = make_state(1.0, 0.0,0.0,0.0, 1.0);
    

    for (int i = 0; i <= grid.getSizeX(); i++){
        for (int j = 0; j <= grid.getSizeY(); j++){
            grid[i,j] = W;
            double x = i*dx, y = j*dx;
            grid._A()[i,j] = vec3{0,0,y};
            //Perturb
            grid[i,j].v.y += amp*sin(2*M_PI*(x+0.5));
            grid._A()[i,j].z -= amp/sqrt(M_PI) * cos(2*M_PI*x);
        }
    }
    grid.boundary = Boundary::Periodic();
    grid.initialize_B_fields();

    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            expected[i,j] = grid[i,j];
        }
    }
    
    assert_divergenceless(grid._A(),dx,dx);
    
    grid.advance(sq4pi);
    
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            expect_close(grid[i,j].B, expected[i,j].B);
        }
    }
    assert_divergenceless(grid._A(),dx,dx);
}

void DRAGON_Test::verify_ct_alfven_wave_3D(){
    double dx = 0.1, amp = 0.01;
    DistGrid3D grid(10,10,10,dx,dx,dx, 2), expected(10,10,10,dx,dx,dx, 2);
    PrimitiveState W = make_state(1.0, 0.0,0.0,0.0, 1.0);
    

    for (int i = 0; i <= grid.getSizeX(); i++){
        for (int j = 0; j <= grid.getSizeY(); j++){
            for (int k = 0; k <= grid.getSizeZ(); k++){
                grid[i,j,k] = W;
                double x = i*dx, y = j*dx, z = k*dx;
                grid._A()[i,j,k] = vec3{0,-z,y};
                //Perturb
                grid[i,j,k].v.y += amp*sin(2*M_PI*(x+0.5));
                grid[i,j,k].v.z += amp*cos(2*M_PI*(x+0.5));
                grid._A()[i,j,k].y -= amp/sqrt(M_PI) * sin(2*M_PI*x);
                grid._A()[i,j,k].z -= amp/sqrt(M_PI) * cos(2*M_PI*x);
            }
        }
    }
    grid.boundary = Boundary::Periodic(Boundary::X);
    grid.initialize_B_fields();

    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            for (int k = 0; k < grid.getSizeZ(); k++){
                expected[i,j,k] = grid[i,j,k];
            }
        }
    }
    
    assert_divergenceless(grid._A(),dx,dx,dx);
    
    grid.advance(10*sq4pi);
    
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            for (int k = 0; k < grid.getSizeZ(); k++){
                expect_close(grid[i,j,k].B, expected[i,j,k].B);
            }
        }
    }
    assert_divergenceless(grid._A(),dx,dx,dx);
    
    DRAGONWING::initialize(0);
}

//Tests above this line, I wrote myself from scratch
//Tests below this line, AI originally drafted.  I have since revised the tests to better match my coding style.

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
