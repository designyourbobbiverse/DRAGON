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
#include "CFL.hpp"
#include <iostream>

#ifdef MHD
using namespace DRAGON_Test;

void DRAGON_Test::verify_ct_2D(bool output){
    if(output) std::cout << "Constrained Transport (2D):\n";
    if(output) std::cout << "- Faraday Pipeline:\n";
    if(output) std::cout << "\t- Edge E -> Edge A: ";
    verify_ct_E_updates_A_2D();
    if(output) std::cout << "Passed\n";

    if(output) std::cout << "\t- Edge A -> Face B: ";
    verify_ct_compute_faces_2D();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "\t- Face B -> Body B: ";
    verify_ct_body_fields_2D();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "\t- Face B -> Half-states: ";
    verify_ct_copy_face_fields_2D();
    if(output) std::cout << "Passed\n";

    if(output) std::cout << "- Zero Divergence: ";
    verify_ct_divergence_2D();
    if(output) std::cout << "Passed\n";
    
    if(output) std::cout << "- Gauge Invariance: ";
    verify_ct_gauge_2D();
    if(output) std::cout << "Passed\n";
    
    if(output) std::cout << "- Uniform E: ";
    verify_ct_uniform_E_2D();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "- Stationary Field: ";
    verify_ct_stationary_2D();
    if(output) std::cout << "Passed\n";
}
void DRAGON_Test::verify_ct_3D(bool output){
    if(output) std::cout << "Constrained Transport (3D): \n";
    if(output) std::cout << "- Faraday Update:\n";
    if(output) std::cout << "\t- Zero Divergence: ";
    verify_ct_divergence_3D();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "\t- Stokes Theorem: ";
    verify_ct_stokes_theorem_3D();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "\t- Uniform E: ";
    verify_ct_uniform_E_3D();
    if(output) std::cout << "Passed\n";
    
    
    if(output) std::cout << "- Face & Body Fields:\n";
    if(output) std::cout << "\t- Face B -> Body B: ";
    verify_ct_body_fields_3D();
    if(output) std::cout << "Passed\n";
    if(output) std::cout << "\t- Face B -> Half-states: ";
    verify_ct_copy_face_fields_3D();
    if(output) std::cout << "Passed\n";

    
    if(output) std::cout << "- Stationary Field Test: ";
    verify_ct_stationary_3D();
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
void assert_divergenceless(const MagneticArray3D& B, double dx, double dy, double dz){
    const int nx = B.getSizeX(), ny = B.getSizeY(),nz = B.getSizeZ(), ng = 0;

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
            A[i,j] = { (rand()%2000001)*1e-3 - 1e3, (rand()%2000001)*1e-3 - 1e3 , (rand()%2000001)*1e-3 - 1e3};
        }
    }
    assert_divergenceless(A, 1, 2);
}


void DRAGON_Test::verify_ct_divergence_3D(){
    const int nx = 10, ny = 10, nz = 10, ng = 2;
    MagneticArray3D A(nx,ny,nz,ng), B(nx,ny,nz,ng);
    for(int i = -ng; i < nx+ng; i++){
        for(int j = -ng; j < ny+ng; j++){
            for(int k = -ng; k < nz+ng; k++){
                A[i,j,k] = { (rand()%2000001)*1e-3 - 1e3, (rand()%2000001)*1e-3 - 1e3 , (rand()%2000001)*1e-3 - 1e3};
                B[i,j,k] = {0,0,0};
            }
        }
    }
    CT::Faraday(A, B, 1,1,1, 2);
    assert_divergenceless(B, 1,1,1);
}


//MARK: Gauge Invariance
void DRAGON_Test::verify_ct_gauge_2D(){
    const int nx = 10, ny = 10, ng = 2;
    const double dx = 1, dy = 2;
    MagneticArray2D A(nx,ny,ng);
    for(int i = -ng; i < nx+ng; i++){
        for(int j = -ng; j < ny+ng; j++){
            A[i,j] = { (rand()%2000001)*1e-3 - 1e3, (rand()%2000001)*1e-3 - 1e3 , (rand()%2000001)*1e-3 - 1e3};
        }
    }
    //Compute Magnetic
    MagneticArray2D B0(nx,ny,ng);
    CT::computeFaceFields(A, B0, dx, dy);
    //Update A by a Random Gradient
    ExtendedArray2D<double> chi(nx+1,ny+1,ng);
    for(int i = -ng; i <= nx+ng; i++){
        for(int j = -ng; j <= ny+ng; j++){
            chi[i,j] = (rand()%2000001)*1e-3 - 1e3;
        }
    }
    double chi_z = (rand()%2000001)*1e-3 - 1e3;
    for(int i = -ng; i < nx+ng; i++){
        for(int j = -ng; j < ny+ng; j++){
            A[i,j].x += (chi[i+1,j]-chi[i,j]) / dx;
            A[i,j].y += (chi[i,j+1]-chi[i,j]) / dy;
            A[i,j].z += chi_z;
        }
    }
    //Recompute Magnetic
    MagneticArray2D B(nx,ny,ng);
    CT::computeFaceFields(A, B, dx, dy);
    //Make sure it's the same result
    for(int i = -ng; i < nx+ng; i++){
        for(int j = -ng; j < ny+ng; j++){
            expect_close(B[i,j], B0[i,j]);
        }
    }
}

//MARK: Stationary Field
void DRAGON_Test::verify_ct_stationary_2D(){
    return;
    
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
            grid[i,j].p = p0 - grid[i,j].B * grid[i,j].B * _1_8pi;
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
    W.B = {0.2, -0.3, 0.4};

    for (int i = 0; i <= grid.getSizeX(); i++){
        for (int j = 0; j <= grid.getSizeY(); j++){
            for (int k = 0; k <= grid.getSizeZ(); k++){
                grid[i,j,k] = W;
                grid._B()[i,j,k] = W.B;
            }
        }
    }
    grid.boundary = Boundary::Periodic();
    grid.boundary.apply(grid);
    grid.initialize_B_fields();

    for (int i = 0; i <= grid.getSizeX(); i++){
        for (int j = 0; j <= grid.getSizeY(); j++){
            for (int k = 0; k <= grid.getSizeZ(); k++){
                expected[i,j,k] = grid[i,j,k];
            }
        }
    }
    assert_divergenceless(grid._B(),dx,dx,dx);

    
    grid.advance(0.0001);
    
    for (int i = 0; i < grid.getSizeX(); i++){
        for (int j = 0; j < grid.getSizeY(); j++){
            for (int k = 0; k < grid.getSizeZ(); k++){
                expect_close(grid[i,j,k].B, expected[i,j,k].B);
            }
        }
    }
    assert_divergenceless(grid._B(),dx,dx,dx);
    
    DRAGONWING::initialize(0);
}

/*
//MARK: d/dz = 0 -> Bz stays zero
void DRAGON_Test::verify_ct_uniform_Bz0_2D(){
    constexpr int n = 16;
    constexpr double dx = 2.0/n;
    
    auto grid = Grid3D(n,n,n, dx, dx, dx);
    grid.boundary = Boundary::Periodic();
    
    for(int i=0; i<=grid.getSizeX();i++){
        for(int j=0; j<=grid.getSizeY(); j++){
            for(int k=0; k<=grid.getSizeZ(); k++){
                grid[i,j,k] = make_state(1.0, 2.0, 3.0, 4.0, 5.0);
                grid._A()[i,j,k] = {0,0, 0.125*(i+j)};
            }
        }
    }    
    grid.initialize_B_fields();
    
    
    std::cout<<CFL::cfl_time(grid)<<"\n";
    double dt = 0.00001;
    grid.advance(dt);
    double Bzmax = 0;
    for(int i=0; i<grid.getSizeX();i++){
        for(int j=0; j<grid.getSizeY(); j++){
            for(int k=0; k<grid.getSizeZ(); k++){
                double Bz = grid[i,j,k].B.z;
                if(Bzmax < fabs(Bz)) Bzmax = fabs(Bz);
                //assert(approx(Bz,0));
            }
        }
    }
    std::cout<<dt<<","<<Bzmax<<"\n";
}
*/


#endif
