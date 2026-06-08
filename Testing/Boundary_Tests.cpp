//
//  Riemann_Tests.cpp
//  DRAGON/Testing
//
//  Created by Bobbie Markwick on 12/05/2026.
//

#include "Testing.hpp"
#include "Boundary.hpp"
#include "Grid.hpp"
#include <iostream>

using namespace DRAGON_Test;
using namespace Boundary;


void DRAGON_Test::verify_boundary(){
    verify_boundary_constructors();
    verify_boundary_reflective();
    verify_boundary_periodic();
    verify_boundary_outflow();
    verify_boundary_fixed();
}

//MARK: Helpers
PrimitiveState make_tagged_state(double tag){
    PrimitiveState W;
    W.rho = tag;
    W.vx = tag + 10.0;
    W.vy = tag + 20.0;
    W.vz = tag + 30.0;
    W.p = tag + 40.0;
    return W;
}
const PrimitiveState G = make_tagged_state(-666);


void fill_1D(Grid1D& grid) {
    double s = 1.0 / grid.getSize();
    for (int i = 0; i < grid.getSize(); i++) {
        grid[i] = make_tagged_state(i*s);
    }
    for (int g = 1; g <= grid.getGhosts(); g++){
        grid[-g] = G;
        grid[grid.getSize()+g-1] = G;
    }
}
void fill_2D(Grid2D& grid) {
    double s = 1.0 / fmax(grid.getSizeX(), grid.getSizeY());
    for (int i = -grid.getGhosts(); i < grid.getSizeX()+grid.getGhosts(); i++) {
        for (int j = -grid.getGhosts(); j < grid.getSizeY()+grid.getGhosts(); j++) {
            grid[i, j] = G;
        }
    }
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int j = 0; j < grid.getSizeY(); j++) {
            grid[i, j] = make_tagged_state(i*s + j*s*s);
        }
    }
    
}
void fill_3D(Grid3D& grid) {
    double s = 1.0 / fmax(fmax(grid.getSizeX(), grid.getSizeY()), grid.getSizeZ());
    for (int i = -grid.getGhosts(); i < grid.getSizeX()+grid.getGhosts(); i++) {
        for (int j = -grid.getGhosts(); j < grid.getSizeY()+grid.getGhosts(); j++) {
            for (int k = -grid.getGhosts(); k < grid.getSizeZ()+grid.getGhosts(); k++) {
                grid[i, j, k] = G;
            }
        }
    }
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int j = 0; j < grid.getSizeY(); j++) {
            for (int k = 0; k < grid.getSizeZ(); k++) {
                grid[i, j, k] = make_tagged_state(i*s + j*s*s + k*s*s*s + 1);
            }
        }
    }
}

//MARK: Boundary Consturction
void DRAGON_Test::verify_boundary_constructors(){
    PrimitiveState W = make_tagged_state(999);

    //Reflective Constructors
    Reflective r0;
    Reflective r1(X);
    Reflective r2(X_negative, false);
    assert(&r0);
    assert(&r1);
    assert(&r2);
    
    //Outflow Constructors
    Outflow o0;
    Outflow o1(X);
    Outflow o2("X");
    Outflow o3(X, false);
    Outflow o4(X_positive, true, true);
    Outflow o5("X+", true, true);
    assert(&o0);
    assert(&o1);
    assert(&o2);
    assert(&o3);
    assert(&o4);
    assert(&o5);

    //Gated Outflow shorthand
    Outflow og0 = Outflow::Gated();
    Outflow og1 = Outflow::Gated(X_positive);
    Outflow og2 = Outflow::Gated(X_positive, false);
    Outflow og3 = Outflow::Gated("X+", false);
    assert(&og0);
    assert(&og1);
    assert(&og2);
    assert(&og3);

    //Periodic Constructors
    Periodic p0;
    Periodic p1(X);
    Periodic p2(Y, false);
    Periodic p3("Z");
    assert(&p0);
    assert(&p1);
    assert(&p2);
    assert(&p3);


    //Fixed State Constructors
    Fixed f0(W);
    Fixed f1(W, X);
    Fixed f2(W, Y, false);
    Fixed f3(W, "YZ");
    assert(&f0);
    assert(&f1);
    assert(&f2);
    assert(&f3);

    //Composite Boundary Constructors
    BoundaryList b0(r1, p2, Outflow(Z));
    BoundaryList b1(Fixed(W, "X-"), Outflow::Gated(X_positive));
    assert(&b0);
    assert(&b1);

    //Test that faces are intiallized properly and also test get_faces()
    assert(Reflective(X).get_faces() == X);
    assert(Reflective("X-").get_faces() == X_negative);
    assert(Reflective("X+").get_faces() == X_positive);
    assert(Reflective(X | Y).get_faces() == (X | Y));
    assert(Periodic(Y).get_faces() == Y);
    assert(Outflow(Z_negative | Z_positive).get_faces() == Z);
}

//MARK: Composition
void DRAGON_Test::verify_boundary_composition(){
    verify_boundary_set_missing_faces_outflow_2D();
    verify_boundary_composition_3D();
    verify_boundary_composition_order();
}
void DRAGON_Test::verify_boundary_set_missing_faces_outflow_2D() {

    Grid2D grid(3, 4, 1.0, 1.0, 1);
    fill_2D(grid);
    auto boundary = BoundaryList( Reflective(X) );
    boundary.apply(grid);
    
    // X should be reflective
    PrimitiveState wxL = grid[0,0]; wxL.vx *= -1;
    expect_close(grid[-1,0], wxL);
    PrimitiveState wxR = grid[2,0]; wxR.vx *= -1;
    expect_close(grid[3,0], wxR);
    // Y should be default outflow
    expect_close(grid[0,-1], grid[0,0]);
    expect_close(grid[0, 4], grid[0,3]);
}

void DRAGON_Test::verify_boundary_composition_3D() {

    Grid3D grid(3, 4, 5, 1.0, 1.0,1.0, 1);
    fill_3D(grid);
    auto boundary = Reflective(Z_negative) +  Periodic(Y) + Reflective(Z_positive);
    boundary.apply(grid);
    

    // X should be default outflow
    expect_close(grid[-1,0,0], grid[0,0,0]);
    expect_close(grid[3,0,0], grid[2,0,0]);
    // Y should be Periodic
    expect_close(grid[0,-1,0], grid[0,3,0]);
    expect_close(grid[0, 4,0], grid[0,0,0]);
    // Z should be reflective
    PrimitiveState wxL = grid[0,0,0]; wxL.vz *= -1;
    expect_close(grid[0,0,-1], wxL);
    PrimitiveState wxR = grid[0,0,4]; wxR.vz *= -1;
    expect_close(grid[0,0,5], wxR);
    //Corner
    PrimitiveState w = grid[0,3,0]; w.vz *= -1;
    expect_close(grid[-1,-1,-1], w);
}
void DRAGON_Test::verify_boundary_composition_order() {
    Grid3D grid(3, 4, 5, 1.0, 1.0,1.0, 1);
    fill_3D(grid);
    BoundaryList boundary = Reflective();
    boundary += Periodic("Y");
    boundary.apply(grid);
    
    // Y should be Periodic
    expect_close(grid[0,-1,0], grid[0,3,0]);
    expect_close(grid[0, 4,0], grid[0,0,0]);
    // Z should be reflective
    PrimitiveState wxL = grid[0,0,0]; wxL.vz *= -1;
    expect_close(grid[0,0,-1], wxL);
    PrimitiveState wxR = grid[0,0,4]; wxR.vz *= -1;
    expect_close(grid[0,0,5], wxR);
    //Corner
    PrimitiveState w = grid[0,3,0]; w.vx *= -1;w.vz *= -1;
    expect_close(grid[-1,-1,-1], w);
}

//MARK: Fixed
void DRAGON_Test::verify_boundary_fixed(){
    verify_boundary_fixed_1D();
    verify_boundary_fixed_2D();
    verify_boundary_fixed_3D();
}
void DRAGON_Test::verify_boundary_fixed_1D() {
    Grid1D grid(4, 1.0, 2);
    PrimitiveState W = make_tagged_state(2112);
    
    fill_1D(grid);

    Boundary::Fixed(W,X).apply(grid);
    //Left
    expect_close(grid[-1], W);
    expect_close(grid[-2], W);
    //Right
    expect_close(grid[4], W);
    expect_close(grid[5], W);
    //Single face filling: Left
    fill_1D(grid);
    Boundary::Fixed(W,"X-").apply(grid);
    expect_close(grid[-1], W);
    expect_close(grid[-2], W);
    expect_close(grid[4], G);
    expect_close(grid[5], G);
    //Single face filling: Right
    fill_1D(grid);
    Boundary::Fixed(W,"X+").apply(grid);
    expect_close(grid[-1], G);
    expect_close(grid[-2], G);
    expect_close(grid[4], W);
    expect_close(grid[5], W);
}
//MARK: Fixed - 2D
void DRAGON_Test::verify_boundary_fixed_2D() {
    Grid2D grid(3, 4, 1.0, 1.0, 1);
    PrimitiveState W = make_tagged_state(2112);

    //X
    fill_2D(grid);
    Boundary::Fixed(W,X,false).apply(grid);
    for (int j = 0; j < grid.getSizeY(); j++) {
        expect_close(grid[-1, j], W);
        expect_close(grid[3, j],  W);
    }
    //No Corners = No Corners
    expect_close(grid[-1,-1],G);
    //Y
    fill_2D(grid);
    Boundary::Fixed(W,"Y").apply(grid);
    for (int i = 0; i < grid.getSizeX(); i++) {
        expect_close(grid[i,-1], W);
        expect_close(grid[i,4],  W);
    }
    //Corner
    fill_2D(grid);
    Boundary::Fixed(W,X | Y).apply(grid);
    expect_close(grid[-1,-1], W);
}
//MARK: Fixed - 3D
void DRAGON_Test::verify_boundary_fixed_3D() {
    Grid3D grid(3, 4, 5, 1.0, 1.0, 1.0, 1);
    PrimitiveState W = make_tagged_state(2112);

    //X
    fill_3D(grid);
    Boundary::Fixed(W,"X").apply(grid);
    for (int j = 0; j < grid.getSizeY(); j++) {
        for (int k = 0; k < grid.getSizeZ(); k++) {
            expect_close(grid[-1, j,k], W);
            expect_close(grid[3, j,k],  W);
        }
    }
    //Y
    fill_3D(grid);
    Boundary::Fixed(W,"Y",false).apply(grid);
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int k = 0; k < grid.getSizeZ(); k++) {
            expect_close(grid[i,-1, k], W);
            expect_close(grid[i,4, k],  W);
        }
    }
    //No corners = no corners
    expect_close(grid[-1,-1,1],G);
    //Z
    fill_3D(grid);
    Boundary::Fixed(W,"Z").apply(grid);
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int j = 0; j < grid.getSizeY(); j++) {
            expect_close(grid[i,j,-1], W);
            expect_close(grid[i,j,5],  W);
        }
    }
    //Corner
    fill_3D(grid);
    Boundary::Fixed(W,X | Y).apply(grid);
    expect_close(grid[-1,-1,-1], W);
}

//MARK: Outflow
void DRAGON_Test::verify_boundary_outflow(){
    verify_boundary_outflow_1D();
    verify_boundary_outflow_2D();
    verify_boundary_outflow_3D();
    verify_boundary_outflow_1D_gated();
    verify_boundary_outflow_2D_gated();
    verify_boundary_outflow_3D_gated_X();
    verify_boundary_outflow_3D_gated_Y();
    verify_boundary_outflow_3D_gated_Z();
}
void DRAGON_Test::verify_boundary_outflow_1D() {
    Grid1D grid(4, 1.0, 2);
    fill_1D(grid);

    Outflow(X).apply(grid);
    //Left
    expect_close(grid[-1], grid[0]);
    expect_close(grid[-2], grid[0]);
    //Right
    expect_close(grid[4], grid[3]);
    expect_close(grid[5], grid[3]);
    //Single face filling: Left
    fill_1D(grid);
    Outflow("X-").apply(grid);
    expect_close(grid[-1], grid[0]);
    expect_close(grid[-2], grid[0]);
    expect_close(grid[4], G);
    expect_close(grid[5], G);
    //Single face filling: Right
    fill_1D(grid);
    Outflow("X+").apply(grid);
    expect_close(grid[4], grid[3]);
    expect_close(grid[5], grid[3]);
    expect_close(grid[-1], G);
    expect_close(grid[-2], G);
}
//MARK: Outflow - 2D
void DRAGON_Test::verify_boundary_outflow_2D() {
    Grid2D grid(3, 4, 1.0, 1.0, 1);
    
    //X
    fill_2D(grid);
    Outflow(X).apply(grid);
    for (int j = 0; j < grid.getSizeY(); j++) {
        expect_close(grid[-1, j], grid[0, j]);
        expect_close(grid[3, j],  grid[2, j]);
    }
    //Y
    fill_2D(grid);
    Outflow(Y,false).apply(grid);
    for (int i = 0; i < grid.getSizeX(); i++) {
        expect_close(grid[i,-1], grid[i,0]);
        expect_close(grid[i,4],  grid[i,3]);
    }
    //No corners = no corners
    expect_close(grid[-1,-1],G);
    //Corner
    fill_2D(grid);
    Outflow("XY").apply(grid);
    expect_close(grid[-1,-1], grid[0,0]);
}
//MARK: Outflow - 3D
void DRAGON_Test::verify_boundary_outflow_3D() {
    Grid3D grid(3, 4, 5, 1.0, 1.0, 1.0, 1);
    
    //X
    fill_3D(grid);
    Outflow(X).apply(grid);
    for (int j = 0; j < grid.getSizeY(); j++) {
        for (int k = 0; k < grid.getSizeZ(); k++) {
            expect_close(grid[-1, j,k], grid[0, j,k]);
            expect_close(grid[3, j,k],  grid[2, j,k]);
        }
    }
    //Y
    fill_3D(grid);
    Outflow(Y_negative).apply(grid);
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int k = 0; k < grid.getSizeZ(); k++) {
            expect_close(grid[i,-1, k], grid[i,0, k]);
            expect_close(grid[i,4, k],  G);
        }
    }
    fill_3D(grid);
    Outflow(Y_positive).apply(grid);
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int k = 0; k < grid.getSizeZ(); k++) {
            expect_close(grid[i,-1, k], G);
            expect_close(grid[i,4, k],  grid[i,3,k]);
        }
    }
    //Z
    fill_3D(grid);
    Outflow("Z",false).apply(grid);
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int j = 0; j < grid.getSizeY(); j++) {
            expect_close(grid[i,j,-1], grid[i,j,0]);
            expect_close(grid[i,j,5],  grid[i,j,4]);
        }
    }
    //No corners = no corners
    expect_close(grid[-1,1,-1],G);
    //Corner
    fill_3D(grid);
    Outflow(X | Y | Z).apply(grid);
    expect_close(grid[-1,-1,-1], grid[0,0,0]);
}


//MARK: Periodic
void DRAGON_Test::verify_boundary_periodic(){
    verify_boundary_periodic_1D();
    verify_boundary_periodic_2D();
    verify_boundary_periodic_3D();
}
void DRAGON_Test::verify_boundary_periodic_1D() {
    Grid1D grid(4, 1.0, 2);
    fill_1D(grid);

    Periodic(X).apply(grid);
    //Left
    expect_close(grid[-1], grid[3]);
    expect_close(grid[-2], grid[2]);
    //Right
    expect_close(grid[4], grid[0]);
    expect_close(grid[5], grid[1]);
}
//MARK: Periodic - 2D
void DRAGON_Test::verify_boundary_periodic_2D() {
    Grid2D grid(3, 4, 1.0, 1.0, 1);
    
    //X
    fill_2D(grid);
    Periodic(X).apply(grid);
    for (int j = 0; j < grid.getSizeY(); j++) {
        expect_close(grid[-1, j], grid[2, j]);
        expect_close(grid[3, j],  grid[0, j]);
    }
    //Y
    fill_2D(grid);
    Periodic(Y,false).apply(grid);
    for (int i = 0; i < grid.getSizeX(); i++) {
        expect_close(grid[i,-1], grid[i,3]);
        expect_close(grid[i,4],  grid[i,0]);
    }
    //No corners = no corners
    expect_close(grid[-1,-1],G);
    //Corner
    fill_2D(grid);
    Periodic(X | Y).apply(grid);
    expect_close(grid[-1,-1], grid[2,3]);
}
//MARK: Periodic - 3D
void DRAGON_Test::verify_boundary_periodic_3D() {
    Grid3D grid(3, 4, 5, 1.0, 1.0, 1.0, 1);
    
    //X
    fill_3D(grid);
    Periodic("X",false).apply(grid);
    for (int j = 0; j < grid.getSizeY(); j++) {
        for (int k = 0; k < grid.getSizeZ(); k++) {
            expect_close(grid[-1, j,k], grid[2, j,k]);
            expect_close(grid[3, j,k],  grid[0, j,k]);
        }
    }
    //No corners = no corners
    expect_close(grid[-1,1,-1],G);
    //Y
    fill_3D(grid);
    Periodic(Y).apply(grid);
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int k = 0; k < grid.getSizeZ(); k++) {
            expect_close(grid[i,-1, k], grid[i,3, k]);
            expect_close(grid[i,4, k],  grid[i,0,k]);
        }
    }
    //Z
    fill_3D(grid);
    Periodic("Z").apply(grid);
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int j = 0; j < grid.getSizeY(); j++) {
            expect_close(grid[i,j,-1], grid[i,j,4]);
            expect_close(grid[i,j,5],  grid[i,j,0]);
        }
    }
    //Corner
    fill_3D(grid);
    Periodic(X | Y | Z).apply(grid);
    expect_close(grid[-1,-1, -1], grid[2,3,4]);
}

//MARK: Reflecting
void DRAGON_Test::verify_boundary_reflective(){
    verify_boundary_reflective_1D();
    verify_boundary_reflective_2D();
    verify_boundary_reflective_3D();
}
void DRAGON_Test::verify_boundary_reflective_1D() {

    Grid1D grid(4, 1.0, 2);
    fill_1D(grid);

    Reflective(X).apply(grid);
    //Left
    PrimitiveState w = grid[0]; w.vx *= -1;
    expect_close(grid[-1], w);
    w = grid[1]; w.vx *= -1;
    expect_close(grid[-2], w);
    //Right
    w = grid[3]; w.vx *= -1;
    expect_close(grid[4], w);
    w = grid[2]; w.vx *= -1;
    expect_close(grid[5], w);
    //Single face filling: Left
    fill_1D(grid);
    Reflective("X-").apply(grid);
    w = grid[0]; w.vx *= -1;
    expect_close(grid[-1], w);
    w = grid[1]; w.vx *= -1;
    expect_close(grid[-2], w);
    //Single face filling: Right
    fill_1D(grid);
    Reflective("X+").apply(grid);
    w = grid[3]; w.vx *= -1;
    expect_close(grid[4], w);
    w = grid[2]; w.vx *= -1;
    expect_close(grid[5], w);
}
//MARK: Reflecting - 2D
void DRAGON_Test::verify_boundary_reflective_2D(){
    Grid2D grid(3, 4, 1.0, 1.0, 1);
    //X
    fill_2D(grid);
    Reflective(X,false).apply(grid);
    for (int j = 0; j < grid.getSizeY(); j++) {
        PrimitiveState w = grid[0,j]; w.vx *= -1;
        expect_close(grid[-1, j], w);
       
        w = grid[2,j]; w.vx *= -1;
        expect_close(grid[3, j], w);
    }
    //No corners = no corners
    expect_close(grid[-1,-1],G);
    //Y
    fill_2D(grid);
    Reflective(Y).apply(grid);
    for (int i = 0; i < grid.getSizeX(); i++) {
        PrimitiveState w = grid[i,0]; w.vy *= -1;
        expect_close(grid[i,-1], w);
       
        w = grid[i,3]; w.vy *= -1;
        expect_close(grid[i,4], w);
    }
    //Corner
    fill_2D(grid);
    Reflective(X | Y).apply(grid);
    auto w = grid[0,0]; w.vx *= -1; w.vy *= -1;
    expect_close(grid[-1,-1], w);
}
//MARK: Reflecting - 3D
void DRAGON_Test::verify_boundary_reflective_3D(){
    Grid3D grid(3, 4,5,1.0, 1.0, 1.0, 1);
    //X
    fill_3D(grid);
    Reflective(X).apply(grid);
    for (int j = 0; j < grid.getSizeY(); j++) {
        for (int k = 0; k < grid.getSizeZ(); k++) {
            PrimitiveState w = grid[0,j,k]; w.vx *= -1;
            expect_close(grid[-1, j,k], w);
            
            w = grid[2,j,k]; w.vx *= -1;
            expect_close(grid[3, j,k], w);
        }
    }
    //Y
    fill_3D(grid);
    Reflective(Y).apply(grid);
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int k = 0; k < grid.getSizeZ(); k++) {
            PrimitiveState w = grid[i,0,k]; w.vy *= -1;
            expect_close(grid[i,-1,k], w);
            
            w = grid[i,3,k]; w.vy *= -1;
            expect_close(grid[i,4,k], w);
        }
    }
    //Z
    fill_3D(grid);
    Reflective("Z",false).apply(grid);
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int j = 0; j < grid.getSizeY(); j++) {
            PrimitiveState w = grid[i,j,0]; w.vz *= -1;
            expect_close(grid[i,j,-1], w);
            
            w = grid[i,j,4]; w.vz *= -1;
            expect_close(grid[i,j,5], w);
        }
    }
    //No corners = no corners
    expect_close(grid[1,-1,-1],G);
    //Corner
    fill_3D(grid);
    Reflective().apply(grid);
    auto w = grid[0,0,0]; w.vx *= -1; w.vy *= -1; w.vz *= -1;
    expect_close(grid[-1,-1,-1], w);
}



//MARK: Gated Outflow
void DRAGON_Test::verify_boundary_outflow_1D_gated(){
    Grid1D grid(4, 1.0, 2);
   
    //Left blocks inflow
    fill_1D(grid);
    grid[0].vx = +5.0;
    Outflow::Gated(X_negative).apply(grid);
    assert(approx(grid[-1].vx, 0.0));
    assert(approx(grid[-2].vx, 0.0));
    //Left allows outflow
    fill_1D(grid);
    grid[0].vx = -5.0;
    Outflow::Gated(X_negative).apply(grid);
    assert(approx(grid[-1].vx, -5));
    assert(approx(grid[-2].vx, -5));
    
    //Right blocks inflow
    fill_1D(grid);
    grid[3].vx = -5.0;
    Outflow::Gated(X_positive).apply(grid);
    assert(approx(grid[4].vx, 0.0));
    assert(approx(grid[5].vx, 0.0));
    //Right allows outflow
    fill_1D(grid);
    grid[3].vx = +5.0;
    Outflow::Gated(X_positive).apply(grid);
    assert(approx(grid[4].vx, +5));
    assert(approx(grid[5].vx, +5));
    
    //Normal components copy
    fill_1D(grid);
    grid[0].vx = +5.0;
    Outflow::Gated(X_negative).apply(grid);
    PrimitiveState w=grid[-1]; w.vx += grid[0].vx;
    expect_close(w, grid[0]);

}

//MARK: Gated Outflow - 2D
void DRAGON_Test::verify_boundary_outflow_2D_gated(){
    Grid2D grid(3,4, 1.0,1.0, 2);
   
    
    //X- blocks inflow
    fill_2D(grid);
    for(int j=0;j<grid.getSizeY();j++) grid[0,j].vx = +5.0;
    Outflow::Gated(X_negative).apply(grid);
    for(int j=0;j<grid.getSizeY();j++){
        assert(approx(grid[-1,j].vx, 0.0));
        assert(approx(grid[-2,j].vx, 0.0));
    }
    //X- allows outflow
    fill_2D(grid);
    for(int j=0;j<grid.getSizeY();j++) grid[0,j].vx = -5.0;
    Outflow::Gated(X_negative).apply(grid);
    for(int j=0;j<grid.getSizeY();j++){
        assert(approx(grid[-1,j].vx, -5));
        assert(approx(grid[-2,j].vx, -5));
    }
    
    //X+ blocks inflow
    fill_2D(grid);
    for(int j=0;j<grid.getSizeY();j++) grid[2,j].vx = -5.0;
    Outflow::Gated(X_positive).apply(grid);
    for(int j=0;j<grid.getSizeY();j++){
        assert(approx(grid[3,j].vx, 0.0));
        assert(approx(grid[4,j].vx, 0.0));
    }
    //X+ allows outflow
    fill_2D(grid);
    for(int j=0;j<grid.getSizeY();j++) grid[2,j].vx = +5.0;
    Outflow::Gated(X_positive).apply(grid);
    for(int j=0;j<grid.getSizeY();j++){
        assert(approx(grid[3,j].vx, +5));
        assert(approx(grid[4,j].vx, +5));
    }
    
    //X - Normal components copy
    fill_2D(grid);
    for(int j=0;j<grid.getSizeY();j++) grid[0,j].vx = +5.0;
    Outflow::Gated(X_negative).apply(grid);
    for(int j=0;j<grid.getSizeY();j++){
        PrimitiveState w=grid[-1,j]; w.vx += grid[0,j].vx;
        expect_close(w, grid[0,j]);
    }
    
    
    //Y- blocks inflow
    fill_2D(grid);
    for(int i=0;i<grid.getSizeX();i++) grid[i,0].vy = +5.0;
    Outflow::Gated("Y-").apply(grid);
    for(int i=0;i<grid.getSizeX();i++){
        assert(approx(grid[i,-1].vy, 0.0));
        assert(approx(grid[i,-2].vy, 0.0));
    }
    //Y- allows outflow
    fill_2D(grid);
    for(int i=0;i<grid.getSizeX();i++) grid[i,0].vy = -5.0;
    Outflow::Gated(Y_negative).apply(grid);
    for(int i=0;i<grid.getSizeX();i++){
        assert(approx(grid[i,-1].vy, -5.0));
        assert(approx(grid[i,-2].vy, -5.0));
    }
    
    //Y+ blocks inflow
    fill_2D(grid);
    for(int i=0;i<grid.getSizeX();i++) grid[i,3].vy = -5.0;
    Outflow::Gated("Y+").apply(grid);
    for(int i=0;i<grid.getSizeX();i++){
        assert(approx(grid[i,4].vy, 0.0));
        assert(approx(grid[i,5].vy, 0.0));
    }
    //Y+ allows outflow
    fill_2D(grid);
    for(int i=0;i<grid.getSizeX();i++) grid[i,3].vy = +5.0;
    Outflow::Gated(Y_positive).apply(grid);
    for(int i=0;i<grid.getSizeX();i++){
        assert(approx(grid[i,4].vy, +5));
        assert(approx(grid[i,5].vy, +5));
    }
    
    //Normal components copy
    fill_2D(grid);
    for(int i=0;i<grid.getSizeX();i++) grid[i,0].vy = +5.0;
    Outflow::Gated(Y_negative).apply(grid);
    for(int i=0;i<grid.getSizeX();i++){
        PrimitiveState w=grid[i,-1]; w.vy += grid[i,0].vy;
        expect_close(w, grid[i,0]);
    }
}
//MARK: Gated Outflow - 3D
void DRAGON_Test::verify_boundary_outflow_3D_gated_X(){
    Grid3D grid(3,4, 5,1.0,1.0,1.0, 2);
    
    
    //X- blocks inflow
    fill_3D(grid);
    for(int j=0;j<grid.getSizeY();j++){
        for(int k=0;k<grid.getSizeZ();k++) grid[0,j,k].vx = +5.0;
    }
    Outflow::Gated(X_negative).apply(grid);
    for(int j=0;j<grid.getSizeY();j++){
        for(int k=0;k<grid.getSizeZ();k++){
            assert(approx(grid[-1,j,k].vx, 0.0));
            assert(approx(grid[-2,j,k].vx, 0.0));
        }
    }
    //X- allows outflow
    fill_3D(grid);
    for(int j=0;j<grid.getSizeY();j++){
        for(int k=0;k<grid.getSizeZ();k++) grid[0,j,k].vx = -5.0;
    }
    Outflow::Gated("X-").apply(grid);
    for(int j=0;j<grid.getSizeY();j++){
        for(int k=0;k<grid.getSizeZ();k++){
            assert(approx(grid[-1,j,k].vx, -5));
            assert(approx(grid[-2,j,k].vx, -5));
        }
    }
    
    //X+ blocks inflow
    fill_3D(grid);
    for(int j=0;j<grid.getSizeY();j++){
        for(int k=0;k<grid.getSizeZ();k++) grid[2,j,k].vx = -5.0;
    }
    Outflow::Gated(X_positive).apply(grid);
    for(int j=0;j<grid.getSizeY();j++){
        for(int k=0;k<grid.getSizeZ();k++){
            assert(approx(grid[3,j,k].vx, 0.0));
            assert(approx(grid[4,j,k].vx, 0.0));
        }
    }
    //X+ allows outflow
    fill_3D(grid);
    for(int j=0;j<grid.getSizeY();j++){
        for(int k=0;k<grid.getSizeZ();k++) grid[2,j,k].vx = +5.0;
    }
    Outflow::Gated("X+").apply(grid);
    for(int j=0;j<grid.getSizeY();j++){
        for(int k=0;k<grid.getSizeZ();k++){
            assert(approx(grid[3,j,k].vx, +5));
            assert(approx(grid[4,j,k].vx, +5));
        }
    }
    
    //X - Normal components copy
    fill_3D(grid);
    for(int j=0;j<grid.getSizeY();j++) {
        for(int k=0;k<grid.getSizeZ();k++) grid[0,j,k].vx = +5.0;
    }
    Outflow::Gated(X_negative).apply(grid);
    for(int j=0;j<grid.getSizeY();j++){
        for(int k=0;k<grid.getSizeZ();k++){
            PrimitiveState w=grid[-1,j,k]; w.vx += grid[0,j,k].vx;
            expect_close(w, grid[0,j,k]);
        }
    }
    
}

void DRAGON_Test::verify_boundary_outflow_3D_gated_Y(){
    Grid3D grid(3,4, 5,1.0,1.0,1.0, 2);
    
    //Y- blocks inflow
    fill_3D(grid);
    for(int i=0;i<grid.getSizeX();i++){
        for(int k=0;k<grid.getSizeZ();k++) grid[i,0,k].vy = +5.0;
    }
    Outflow::Gated(Y_negative).apply(grid);
    for(int i=0;i<grid.getSizeX();i++){
        for(int k=0;k<grid.getSizeZ();k++){
            assert(approx(grid[i,-1,k].vy, 0.0));
            assert(approx(grid[i,-2,k].vy, 0.0));
        }
    }
    //Y- allows outflow
    fill_3D(grid);
    for(int i=0;i<grid.getSizeX();i++){
        for(int k=0;k<grid.getSizeZ();k++) grid[i,0,k].vy = -5.0;
    }
    Outflow::Gated(Y_negative).apply(grid);
    for(int i=0;i<grid.getSizeX();i++){
        for(int k=0;k<grid.getSizeZ();k++){
            assert(approx(grid[i,-1,k].vy, -5.0));
            assert(approx(grid[i,-2,k].vy, -5.0));
        }
    }
    
    //Y+ blocks inflow
    fill_3D(grid);
    for(int i=0;i<grid.getSizeX();i++) {
        for(int k=0;k<grid.getSizeZ();k++) grid[i,3,k].vy = -5.0;
    }
    Outflow::Gated(Y_positive).apply(grid);
    for(int i=0;i<grid.getSizeX();i++){
        for(int k=0;k<grid.getSizeZ();k++){
            assert(approx(grid[i,4,k].vy, 0.0));
            assert(approx(grid[i,5,k].vy, 0.0));
        }
    }
    //Y+ allows outflow
    fill_3D(grid);
    for(int i=0;i<grid.getSizeX();i++) {
        for(int k=0;k<grid.getSizeZ();k++) grid[i,3,k].vy = +5.0;
    }
    Outflow::Gated(Y_positive).apply(grid);
    for(int i=0;i<grid.getSizeX();i++){
        for(int k=0;k<grid.getSizeZ();k++){
            assert(approx(grid[i,4,k].vy, +5));
            assert(approx(grid[i,5,k].vy, +5));
        }
    }
    
    //Normal components copy
    fill_3D(grid);
    for(int i=0;i<grid.getSizeX();i++){
        for(int k=0;k<grid.getSizeZ();k++) grid[i,0,k].vy = +5.0;
    }
    Outflow::Gated(Y_negative).apply(grid);
    for(int i=0;i<grid.getSizeX();i++){
        for(int k=0;k<grid.getSizeZ();k++){
            PrimitiveState w=grid[i,-1,k]; w.vy += grid[i,0,k].vy;
            expect_close(w, grid[i,0,k]);
        }
    }
    
}

void DRAGON_Test::verify_boundary_outflow_3D_gated_Z(){
    Grid3D grid(3,4, 5,1.0,1.0,1.0, 2);
    //Z- blocks inflow
    fill_3D(grid);
    for(int i=0;i<grid.getSizeX();i++){
        for(int j=0;j<grid.getSizeY();j++) grid[i,j,0].vz = +5.0;
    }
    Outflow::Gated("Z-").apply(grid);
    for(int i=0;i<grid.getSizeX();i++){
        for(int j=0;j<grid.getSizeY();j++){
            assert(approx(grid[i,j,-1].vz, 0.0));
            assert(approx(grid[i,j,-2].vz, 0.0));
        }
    }
    //Z- allows outflow
    fill_3D(grid);
    for(int i=0;i<grid.getSizeX();i++){
        for(int j=0;j<grid.getSizeY();j++) grid[i,j,0].vz = -5.0;
    }
    Outflow::Gated(Z_negative).apply(grid);
    for(int i=0;i<grid.getSizeX();i++){
        for(int j=0;j<grid.getSizeY();j++){
            assert(approx(grid[i,j,-1].vz, -5.0));
            assert(approx(grid[i,j,-2].vz, -5.0));
        }
    }
    //Z+ blocks inflow
    fill_3D(grid);
    for(int i=0;i<grid.getSizeX();i++) {
        for(int j=0;j<grid.getSizeY();j++) grid[i,j,4].vz = -5.0;
    }
    Outflow::Gated("Z+").apply(grid);
    for(int i=0;i<grid.getSizeX();i++){
        for(int j=0;j<grid.getSizeY();j++){
            assert(approx(grid[i,j,5].vz, 0.0));
            assert(approx(grid[i,j,6].vz, 0.0));
        }
    }
    //Z+ allows outflow
    fill_3D(grid);
    for(int i=0;i<grid.getSizeX();i++) {
        for(int j=0;j<grid.getSizeY();j++) grid[i,j,4].vz = +5.0;
    }
    Outflow::Gated(Z_positive).apply(grid);
    for(int i=0;i<grid.getSizeX();i++){
        for(int j=0;j<grid.getSizeY();j++){
            assert(approx(grid[i,j,5].vz, +5));
            assert(approx(grid[i,j,6].vz, +5));
        }
    }
    
    //Normal components copy
    fill_3D(grid);
    for(int i=0;i<grid.getSizeX();i++){
        for(int j=0;j<grid.getSizeY();j++) grid[i,j,0].vz = +5.0;
    }
    Outflow::Gated(Z_negative).apply(grid);
    for(int i=0;i<grid.getSizeX();i++){
        for(int j=0;j<grid.getSizeY();j++){
            PrimitiveState w=grid[i,j,-1]; w.vz += grid[i,j,0].vz;
            expect_close(w, grid[i,j,0]);
        }
    }
}


