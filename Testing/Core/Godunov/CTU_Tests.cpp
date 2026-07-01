//
//  CTU_Tests.cpp
//  DRAGON/Testing/Core/Godunov
//
//  Created by Bobbie Markwick on 23/06/2026.
//

#include "Testing.hpp"
#include "DistGrid.hpp"
#include "CFL.hpp"
#include "DragonWing.hpp"
#include <iostream>

using namespace DRAGON_Test;
using namespace Boundary;

//MARK: Diagonal Contact
void DRAGON_Test::verify_ctu_diagonal_contact_2D() {
    const int nx = 40, ny = 40, ghosts = 2;
    const double dx = 1.0 / nx;
    const double dy = 1.0 / ny;

    Grid2D grid(nx, ny, dx, dy, ghosts);
    Grid2D initial(nx, ny, dx, dy, ghosts);

    grid.boundary = Periodic("XY");
    initial.boundary = Periodic("XY");

    const double rhoL = 2.0;
    const double rhoR = 1.0;
    const double vx = 1.0;
    const double vy = 1.0;
    const double p  = 1.0;

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            double x = (i + 0.5) * dx;
            double y = (j + 0.5) * dy;
            double rho = ((x - y) > 0.0) ? rhoL : rhoR;
            grid[i,j] = make_state(rho, vx, vy, 0.0, p);
            initial[i,j] = grid[i,j];
        }
    }

    double dt = 0.4 * CFL::cfl_time(grid);
    for (int n = 0; n < 20; ++n) {
        grid.advance_unsplit(dt);
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) { //Physicality
                assert((grid[i,j].isPhysical()));
            }
        }
    }

    // For v.x = v.y and interface x-y = const,
    // the discontinuity should advect parallel to itself.
    // The exact solution is stationary in this frame because:
    // d/dt(x-y) = v.x - v.y = 0
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            assert(approx(grid[i,j].v.x, vx, 1e-10, 1e-10));
            assert(approx(grid[i,j].v.y, vy, 1e-10, 1e-10));
            assert(approx(grid[i,j].p,  p,  1e-10, 1e-10));
        }
    }
}

void DRAGON_Test::verify_ctu_diagonal_contact_3D() {
    const int nx = 24, ny = 24, nz = 24, ghosts = 2;
    const double dx = 1.0 / nx;
    const double dy = 1.0 / ny;
    const double dz = 1.0 / nz;

    Grid3D grid(nx, ny, nz, dx, dy, dz, ghosts);
    Grid3D initial(nx, ny, nz, dx, dy, dz, ghosts);

    grid.boundary = Periodic(X|Y|Z);
    initial.boundary = Periodic(X|Y|Z);

    const double rhoL = 2.0;
    const double rhoR = 1.0;
    const double vx = 1.0;
    const double vy = 1.0;
    const double vz = 1.0;
    const double p  = 1.0;

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                double x = (i + 0.5) * dx;
                double y = (j + 0.5) * dy;
                double z = (k + 0.5) * dz;

                double rho = ((x - y + 0.5*z) > 0.25) ? rhoL : rhoR;

                grid[i,j,k] = make_state(rho, vx, vy, vz, p);
                initial[i,j,k] = grid[i,j,k];
            }
        }
    }

    double dt = 0.15 * CFL::cfl_time(grid);

    for (int n = 0; n < 20; ++n) {
        grid.advance_unsplit(dt);
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                for (int k = 0; k < nz; ++k) {
                    assert((grid[i,j,k].isPhysical()));
                }
            }
        }
    }
}
//MARK: Blast Helpers
constexpr double blast_sym_tol = 1e-6;
static void expect_blast_mirror_x(const PrimitiveState& L, const PrimitiveState& R, double rel = blast_sym_tol, double abs = blast_sym_tol) {
    assert(approx(L.rho, R.rho, rel, abs));
    assert(approx(L.p, R.p, rel, abs));
    assert(approx(L.v.x, -R.v.x, rel, abs));
    assert(approx(L.v.y, R.v.y, rel, abs));
    assert(approx(L.v.z, R.v.z, rel, abs));
#ifdef MHD
    assert(approx(L.B.x, -R.B.x, rel, abs));
    assert(approx(L.B.y, R.B.y, rel, abs));
    assert(approx(L.B.z, R.B.z, rel, abs));
#endif
}

static void expect_blast_mirror_y(const PrimitiveState& L, const PrimitiveState& R, double rel = blast_sym_tol, double abs = blast_sym_tol) {
    assert(approx(L.rho, R.rho, rel, abs));
    assert(approx(L.p, R.p, rel, abs));
    assert(approx(L.v.x, R.v.x, rel, abs));
    assert(approx(L.v.y, -R.v.y, rel, abs));
    assert(approx(L.v.z, R.v.z, rel, abs));
#ifdef MHD
    assert(approx(L.B.x, R.B.x, rel, abs));
    assert(approx(L.B.y, -R.B.y, rel, abs));
    assert(approx(L.B.z, R.B.z, rel, abs));
#endif
}

static void expect_blast_mirror_z(const PrimitiveState& L, const PrimitiveState& R, double rel = blast_sym_tol, double abs = blast_sym_tol) {
    assert(approx(L.rho, R.rho, rel, abs));
    assert(approx(L.p, R.p, rel, abs));
    assert(approx(L.v.x, R.v.x, rel, abs));
    assert(approx(L.v.y, R.v.y, rel, abs));
    assert(approx(L.v.z, -R.v.z, rel, abs));
#ifdef MHD
    assert(approx(L.B.x, R.B.x, rel, abs));
    assert(approx(L.B.y, R.B.y, rel, abs));
    assert(approx(L.B.z, -R.B.z, rel, abs));
#endif
}

//MARK: Blast Tests
void DRAGON_Test::verify_ctu_blast_2D() {
    const int nx = 64, ny = 64;
    const double dx = 1.0 / nx, dy = 1.0 / ny;

    DistGrid2D grid(nx, ny, dx, dy);
    DistGrid2D initial(nx, ny, dx, dy);

    grid.boundary = Outflow();
    
    const double rho0 = 1.0;
    const double p_ambient = 0.1;
    const double p_blast   = 10.0;
    const double r0 = 0.1;

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            double x = (i + 0.5) * dx - 0.5;
            double y = (j + 0.5) * dy - 0.5;
            double r = sqrt(x*x + y*y);

            double p = (r < r0) ? p_blast : p_ambient;

            grid[i,j] = make_state(rho0, 0.0, 0.0, 0.0, p);
            initial[i,j] = grid[i,j];
        }
    }

    grid.advance(0.05);

    //Verify Physicality
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            assert((grid[i,j].isPhysical()));

        }
    }
    // Mass Conservation
    double mass0 = 0, mass1 = 0;
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            mass0 += initial[i,j].rho;
            mass1 += grid[i,j].rho;
        }
    }
    assert(approx(mass1, mass0, 1e-10, 1e-10));
    // The center should depressurize and nearby gas should start moving outward.
    assert((grid[nx/2, ny/2].p < p_blast));
    // Symmetry across the blast center.
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            expect_blast_mirror_x(grid[i,j], grid[nx-1-i,j]);
            expect_blast_mirror_y(grid[i,j], grid[i,ny-1-j]);
        }
    }
    
    DRAGONWING::initialize(0);

}

void DRAGON_Test::verify_ctu_blast_3D() {
    const int nx = 64, ny = 64, nz = 64;
    const double dx = 1.0 / nx, dy = 1.0 / ny, dz = 1.0 / nz;

    DistGrid3D grid(nx, ny, nz, dx, dy, dz);
    DistGrid3D initial(nx, ny,nz, dx, dy, dz);

    grid.boundary = Outflow();


    const double rho0 = 1.0;
    const double p_ambient = 0.1;
    const double p_blast   = 10.0;
    const double r0 = 0.1;

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                double x = (i + 0.5) * dx - 0.5;
                double y = (j + 0.5) * dy - 0.5;
                double z = (k + 0.5) * dz - 0.5;
                double r = sqrt(x*x + y*y + z*z);
                
                double p = (r < r0) ? p_blast : p_ambient;
                
                grid[i,j,k] = make_state(rho0, 0.0, 0.0, 0.0, p);
                initial[i,j,k] = grid[i,j,k];
            }
        }
    }

    grid.advance(0.02);

    //Verify Physicality
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                assert((grid[i,j,k].isPhysical()));
            }
        }
    }
    // Mass Conservation
    double mass0 = 0, mass1 = 0;
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                mass0 += initial[i,j,k].rho;
                mass1 += grid[i,j,k].rho;
            }
        }
    }
    assert(approx(mass1, mass0, 1e-10, 1e-10));
    // The center should depressurize and nearby gas should start moving outward.
    assert((grid[nx/2, ny/2, nz/2].p < p_blast));
    // Symmetry across the blast center.
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                expect_blast_mirror_x(grid[i,j,k], grid[nx-1-i,j,k]);
                expect_blast_mirror_y(grid[i,j,k], grid[i,ny-1-j,k]);
                expect_blast_mirror_z(grid[i,j,k], grid[i,j,nz-1-k]);
            }
        }
    }
    
    DRAGONWING::initialize(0);
    
}
