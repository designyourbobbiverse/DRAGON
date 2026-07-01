//
//  Boundary_Tests.cpp
//  DRAGON/Testing/Core/Boundary
//
//  Created by Bobbie Markwick on 29/05/2026.
//

#include "Testing.hpp"
#include "Boundary.hpp"
#include "Grid.hpp"
#include <iostream>

using namespace DRAGON_Test;
using namespace Boundary;


void DRAGON_Test::verify_boundary(bool output){
    if(output) std::cout << "Boundary Constrction:\n";
    verify_boundary_constructors();
    //Outflow
    if(output) std::cout << "Outflow Boundary:\n";
    verify_boundary_outflow();
    //Reflective
    if(output) std::cout << "Reflective Boundary:\n";
    verify_boundary_reflective();
    //Periodic
    if(output) std::cout << "Periodic Boundary:\n";
    verify_boundary_periodic();
    //Fixed
    if(output) std::cout << "Fixed State Boundary:\n";
    verify_boundary_fixed();
    //Ignore
    if(output) std::cout << "No-op Boundary:\n";
    verify_boundary_ignore();
    //Composition
    if(output) std::cout << "Boundary Composition:\n";
    verify_boundary_composition();

    if(output) std::cout << "All Boundary Setup tests passed\n\n";
}

//MARK: Helpers
vec3 make_tagged_vec(double tag){
    return vec3(tag + 100, tag + 200, tag + 300);
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
#ifdef MHD
    for (int i = -grid.getGhosts(); i <= grid.getSizeX()+grid.getGhosts(); i++) {
        for (int j = -grid.getGhosts(); j <= grid.getSizeY()+grid.getGhosts(); j++) {
            grid._A()[i, j] = make_tagged_vec(10*i + j);
        }
    }
#endif
    
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
#ifdef MHD
    for (int i = -grid.getGhosts(); i <= grid.getSizeX()+grid.getGhosts(); i++) {
        for (int j = -grid.getGhosts(); j <= grid.getSizeY()+grid.getGhosts(); j++) {
            for (int k = -grid.getGhosts(); k <= grid.getSizeZ()+grid.getGhosts(); k++) {
                grid._A()[i, j, k] = make_tagged_vec(100*i + 10*j + k);
            }
        }
    }
#endif
}

//MARK: Boundary Consturction
void DRAGON_Test::verify_boundary_constructors(bool output){
    PrimitiveState W = make_tagged_state(999);

    //Reflective Constructors
    if(output) std::cout<<"- Reflective: ";
    Reflective r0;
    Reflective r1(X);
    Reflective r2(X_negative, false);
    assert(&r0);
    assert(&r1);
    assert(&r2);
    if(output) std::cout<<"Passed\n";

    //Outflow Constructors
    if(output) std::cout<<"- Outflow: ";
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
    if(output) std::cout<<"Passed\n";


    //Gated Outflow shorthand
    if(output) std::cout<<"- Gated Outflow: ";
    Outflow og0 = Outflow::Gated();
    Outflow og1 = Outflow::Gated(X_positive);
    Outflow og2 = Outflow::Gated(X_positive, false);
    Outflow og3 = Outflow::Gated("X+", false);
    assert(&og0);
    assert(&og1);
    assert(&og2);
    assert(&og3);
    if(output) std::cout<<"Passed\n";

    //Periodic Constructors
    if(output) std::cout<<"- Periodic: ";
    Periodic p0;
    Periodic p1(X);
    Periodic p2(Y, false);
    Periodic p3("Z");
    assert(&p0);
    assert(&p1);
    assert(&p2);
    assert(&p3);
    if(output) std::cout<<"Passed\n";

    //Fixed State Constructors
    if(output) std::cout<<"- Fixed State: ";
    Fixed f0(W);
    Fixed f1(W, X);
    Fixed f2(W, Y, false);
    Fixed f3(W, "YZ");
    assert(&f0);
    assert(&f1);
    assert(&f2);
    assert(&f3);
    if(output) std::cout<<"Passed\n";

    //Composite Boundary Constructors
    if(output) std::cout<<"- Composite: ";
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
    if(output) std::cout<<"Passed\n";

}

//MARK: Composition
void DRAGON_Test::verify_boundary_composition(bool output){
    if(output) std::cout<<"- 2D: ";
    verify_boundary_set_missing_faces_outflow_2D();
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- 3D: ";
    verify_boundary_composition_3D();
    verify_boundary_composition_order();
    verify_boundary_composition_overlap_order();
    if(output) std::cout<<"Passed\n";
}
void reflectX2D(PrimitiveState& w, bool conductive=true);
void DRAGON_Test::verify_boundary_set_missing_faces_outflow_2D() {

    Grid2D grid(3, 4, 1.0, 1.0, 1);
    fill_2D(grid);
    auto boundary = BoundaryList( Reflective(X) );
    boundary.apply(grid);
    
    // X should be reflective
    PrimitiveState wxL = grid[0,0]; reflectX2D(wxL);
    expect_close(grid[-1,0], wxL);
    PrimitiveState wxR = grid[2,0]; reflectX2D(wxR);
    expect_close(grid[3,0], wxR);
    // Y should be default outflow
    expect_close(grid[0,-1], grid[0,0]);
    expect_close(grid[0, 4], grid[0,3]);
}

void reflectZ3D(PrimitiveState& w);
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
    PrimitiveState wxL = grid[0,0,0]; reflectZ3D(wxL);
    expect_close(grid[0,0,-1], wxL);
    PrimitiveState wxR = grid[0,0,4]; reflectZ3D(wxR);
    expect_close(grid[0,0,5], wxR);
    //Corner
    PrimitiveState w = grid[0,3,0]; reflectZ3D(w);
    expect_close(grid[-1,-1,-1], w);
}
void reflectX3D(PrimitiveState& w);
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
    PrimitiveState wxL = grid[0,0,0]; reflectZ3D(wxL);
    expect_close(grid[0,0,-1], wxL);
    PrimitiveState wxR = grid[0,0,4]; reflectZ3D(wxR);
    expect_close(grid[0,0,5], wxR);
    //Corner
    PrimitiveState w = grid[0,3,0]; reflectX3D(w);reflectZ3D(w);
    expect_close(grid[-1,-1,-1], w);
}
void DRAGON_Test::verify_boundary_composition_overlap_order() {
    Grid2D grid(3, 4, 1.0, 1.0, 1);
    fill_2D(grid);
    PrimitiveState W = make_tagged_state(2112);
    BoundaryList boundary = Fixed(W, X) + Reflective(X_negative);
    boundary.apply(grid);

    PrimitiveState reflected = grid[0,0]; reflectX2D(reflected);
    expect_close(grid[-1,0], reflected);
    expect_close(grid[3,0], W);
}


//MARK: Ignore (No-op) Boundary
void DRAGON_Test::verify_boundary_ignore(bool output){
    if(output) std::cout<<"- 1D: ";
    verify_boundary_ignore_1D();
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- 2D: ";
    verify_boundary_ignore_2D();
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- 3D: ";
    verify_boundary_ignore_3D();
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- Implicit Outflow: ";
    verify_boundary_ignore_blocks_implicit_outflow();
    if(output) std::cout<<"Passed\n";
}

void DRAGON_Test::verify_boundary_ignore_1D(){
    Grid1D grid(4, 1.0, 2);
    fill_1D(grid);

    Ignore(X).apply(grid);

    for(int i = 0; i < grid.getSize(); i++)
        expect_close(grid[i], make_tagged_state(i * (1.0 / grid.getSize())));
    expect_close(grid[-1], G);
    expect_close(grid[-2], G);
    expect_close(grid[4], G);
    expect_close(grid[5], G);
}

void DRAGON_Test::verify_boundary_ignore_2D(){
    Grid2D grid(3, 4, 1.0, 1.0, 1);
    fill_2D(grid);

    Ignore(X | Y).apply(grid);

    for(int i = -grid.getGhosts(); i < grid.getSizeX() + grid.getGhosts(); i++) {
        for(int j = -grid.getGhosts(); j < grid.getSizeY() + grid.getGhosts(); j++) {
            if(i < 0 || i >= grid.getSizeX() || j < 0 || j >= grid.getSizeY()) expect_close(grid[i,j], G);
        }
    }
#ifdef MHD
    for(int i = -grid.getGhosts(); i <= grid.getSizeX() + grid.getGhosts(); i++) {
        for(int j = -grid.getGhosts(); j <= grid.getSizeY() + grid.getGhosts(); j++) {
            expect_close(grid._A()[i,j], make_tagged_vec(10*i + j));
        }
    }
#endif
}

void DRAGON_Test::verify_boundary_ignore_3D(){
    Grid3D grid(3, 4, 5, 1.0, 1.0, 1.0, 1);
    fill_3D(grid);

    Ignore(X | Y | Z).apply(grid);

    for(int i = -grid.getGhosts(); i < grid.getSizeX() + grid.getGhosts(); i++) {
        for(int j = -grid.getGhosts(); j < grid.getSizeY() + grid.getGhosts(); j++) {
            for(int k = -grid.getGhosts(); k < grid.getSizeZ() + grid.getGhosts(); k++) {
                if(i < 0 || i >= grid.getSizeX() || j < 0 || j >= grid.getSizeY() || k < 0 || k >= grid.getSizeZ()) expect_close(grid[i,j,k], G);
            }
        }
    }
#ifdef MHD
    for(int i = -grid.getGhosts(); i <= grid.getSizeX() + grid.getGhosts(); i++) {
        for(int j = -grid.getGhosts(); j <= grid.getSizeY() + grid.getGhosts(); j++) {
            for(int k = -grid.getGhosts(); k <= grid.getSizeZ() + grid.getGhosts(); k++) {
                expect_close(grid._A()[i,j,k], make_tagged_vec(100*i + 10*j + k));
            }
        }
    }
#endif
}

void DRAGON_Test::verify_boundary_ignore_blocks_implicit_outflow(){
    Grid2D grid(3, 4, 1.0, 1.0, 1);
    fill_2D(grid);

    BoundaryList(Ignore(X)).apply(grid);

    //Make sure it blocks the X outflow
    for(int j = 0; j < grid.getSizeY(); j++) {
        expect_close(grid[-1,j], G);
        expect_close(grid[3,j], G);
    }
    //Make sure allows the Y outflow
    for(int i = 0; i < grid.getSizeX(); i++) {
        expect_close(grid[i,-1], grid[i,0]);
        expect_close(grid[i,4], grid[i,3]);
    }
}



