
//
//  Outflow_Tests.cpp
//  DRAGON/Testing/Core/Boundary
//
//  Created by Bobbie Markwick on 20/06/2026.
//

#include "Testing.hpp"
#include "Boundary.hpp"
#include "Grid.hpp"
#include <iostream>

using namespace DRAGON_Test;
using namespace Boundary;



void DRAGON_Test::verify_boundary_outflow(bool output){
    if(output) std::cout<<"- 1D: ";
    verify_boundary_outflow_1D();
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- 2D: ";
    verify_boundary_outflow_2D();
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- 3D: ";
    verify_boundary_outflow_3D();
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- Gated 1D: ";
    verify_boundary_outflow_1D_gated();
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- Gated 2D: ";
    verify_boundary_outflow_2D_gated();
    if(output) std::cout<<"Passed\n";
    if(output) std::cout<<"- Gated 3D: ";
    verify_boundary_outflow_3D_gated_X();
    verify_boundary_outflow_3D_gated_Y();
    verify_boundary_outflow_3D_gated_Z();
    if(output) std::cout<<"Passed\n";
}



//MARK: Helpers
static vec3 make_tagged_vec(double tag){ return vec3(tag + 100, tag + 200, tag + 300); }

static PrimitiveState G = make_tagged_state(-666);

void fill_1D(Grid1D& grid);
void fill_2D(Grid2D& grid);
void fill_3D(Grid3D& grid);

//MARK: Outflow - 1D
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
#ifdef MHD
      
#endif
    }
    //Y
    fill_2D(grid);
    Outflow(Y,false).apply(grid);
    for (int i = 0; i < grid.getSizeX(); i++) {
        expect_close(grid[i,-1], grid[i,0]);
        expect_close(grid[i,4],  grid[i,3]);
#ifdef MHD
       
#endif
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
#ifdef MHD
            
#endif
        }
    }
    //Y
    fill_3D(grid);
    Outflow(Y_negative).apply(grid);
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int k = 0; k < grid.getSizeZ(); k++) {
            expect_close(grid[i,-1, k], grid[i,0, k]);
            expect_close(grid[i,4, k],  G);
#ifdef MHD
           
#endif
        }
    }
    fill_3D(grid);
    Outflow(Y_positive).apply(grid);
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int k = 0; k < grid.getSizeZ(); k++) {
            expect_close(grid[i,-1, k], G);
            expect_close(grid[i,4, k],  grid[i,3,k]);
#ifdef MHD
            
#endif
        }
    }
    //Z
    fill_3D(grid);
    Outflow("Z",false).apply(grid);
    for (int i = 0; i < grid.getSizeX(); i++) {
        for (int j = 0; j < grid.getSizeY(); j++) {
            expect_close(grid[i,j,-1], grid[i,j,0]);
            expect_close(grid[i,j,5],  grid[i,j,4]);
#ifdef MHD
           
#endif
        }
    }
    //No corners = no corners
    expect_close(grid[-1,1,-1],G);
    //Corner
    fill_3D(grid);
    Outflow(X | Y | Z).apply(grid);
    expect_close(grid[-1,-1,-1], grid[0,0,0]);
}



//MARK: Gated Outflow - 1D
void DRAGON_Test::verify_boundary_outflow_1D_gated(){
    Grid1D grid(4, 1.0, 2);
   
    //Left blocks inflow
    fill_1D(grid);
    grid[0].v.x = +5.0;
    Outflow::Gated(X_negative).apply(grid);
    assert(approx(grid[-1].v.x, 0.0));
    assert(approx(grid[-2].v.x, 0.0));
    PrimitiveState blocked = grid[-1]; blocked.v.x = grid[0].v.x;
    expect_close(blocked, grid[0]);
    //Left allows outflow
    fill_1D(grid);
    grid[0].v.x = -5.0;
    Outflow::Gated(X_negative).apply(grid);
    assert(approx(grid[-1].v.x, -5));
    assert(approx(grid[-2].v.x, -5));
    
    //Right blocks inflow
    fill_1D(grid);
    grid[3].v.x = -5.0;
    Outflow::Gated(X_positive).apply(grid);
    assert(approx(grid[4].v.x, 0.0));
    assert(approx(grid[5].v.x, 0.0));
    //Right allows outflow
    fill_1D(grid);
    grid[3].v.x = +5.0;
    Outflow::Gated(X_positive).apply(grid);
    assert(approx(grid[4].v.x, +5));
    assert(approx(grid[5].v.x, +5));
    
    //Normal components copy
    fill_1D(grid);
    grid[0].v.x = +5.0;
    Outflow::Gated(X_negative).apply(grid);
    PrimitiveState w=grid[-1]; w.v.x += grid[0].v.x;
    expect_close(w, grid[0]);

}

//MARK: Gated Outflow - 2D
void DRAGON_Test::verify_boundary_outflow_2D_gated(){
    Grid2D grid(3,4, 1.0,1.0, 2);
   
    
    //X- blocks inflow
    fill_2D(grid);
    for(int j=0;j<grid.getSizeY();j++) grid[0,j].v.x = +5.0;
    Outflow::Gated(X_negative).apply(grid);
    for(int j=0;j<grid.getSizeY();j++){
        assert(approx(grid[-1,j].v.x, 0.0));
        assert(approx(grid[-2,j].v.x, 0.0));
        PrimitiveState blocked = grid[-1,j]; blocked.v.x = grid[0,j].v.x;
        expect_close(blocked, grid[0,j]);
    }
    //X- allows outflow
    fill_2D(grid);
    for(int j=0;j<grid.getSizeY();j++) grid[0,j].v.x = -5.0;
    Outflow::Gated(X_negative).apply(grid);
    for(int j=0;j<grid.getSizeY();j++){
        assert(approx(grid[-1,j].v.x, -5));
        assert(approx(grid[-2,j].v.x, -5));
    }
    
    //X+ blocks inflow
    fill_2D(grid);
    for(int j=0;j<grid.getSizeY();j++) grid[2,j].v.x = -5.0;
    Outflow::Gated(X_positive).apply(grid);
    for(int j=0;j<grid.getSizeY();j++){
        assert(approx(grid[3,j].v.x, 0.0));
        assert(approx(grid[4,j].v.x, 0.0));
    }
    //X+ allows outflow
    fill_2D(grid);
    for(int j=0;j<grid.getSizeY();j++) grid[2,j].v.x = +5.0;
    Outflow::Gated(X_positive).apply(grid);
    for(int j=0;j<grid.getSizeY();j++){
        assert(approx(grid[3,j].v.x, +5));
        assert(approx(grid[4,j].v.x, +5));
    }
    
    //X - Normal components copy
    fill_2D(grid);
    for(int j=0;j<grid.getSizeY();j++) grid[0,j].v.x = +5.0;
    Outflow::Gated(X_negative).apply(grid);
    for(int j=0;j<grid.getSizeY();j++){
        PrimitiveState w=grid[-1,j]; w.v.x += grid[0,j].v.x;
        expect_close(w, grid[0,j]);
    }
    
    
    //Y- blocks inflow
    fill_2D(grid);
    for(int i=0;i<grid.getSizeX();i++) grid[i,0].v.y = +5.0;
    Outflow::Gated("Y-").apply(grid);
    for(int i=0;i<grid.getSizeX();i++){
        assert(approx(grid[i,-1].v.y, 0.0));
        assert(approx(grid[i,-2].v.y, 0.0));
        PrimitiveState blocked = grid[i,-1]; blocked.v.y = grid[i,0].v.y;
        expect_close(blocked, grid[i,0]);
    }
    //Y- allows outflow
    fill_2D(grid);
    for(int i=0;i<grid.getSizeX();i++) grid[i,0].v.y = -5.0;
    Outflow::Gated(Y_negative).apply(grid);
    for(int i=0;i<grid.getSizeX();i++){
        assert(approx(grid[i,-1].v.y, -5.0));
        assert(approx(grid[i,-2].v.y, -5.0));
    }
    
    //Y+ blocks inflow
    fill_2D(grid);
    for(int i=0;i<grid.getSizeX();i++) grid[i,3].v.y = -5.0;
    Outflow::Gated("Y+").apply(grid);
    for(int i=0;i<grid.getSizeX();i++){
        assert(approx(grid[i,4].v.y, 0.0));
        assert(approx(grid[i,5].v.y, 0.0));
    }
    //Y+ allows outflow
    fill_2D(grid);
    for(int i=0;i<grid.getSizeX();i++) grid[i,3].v.y = +5.0;
    Outflow::Gated(Y_positive).apply(grid);
    for(int i=0;i<grid.getSizeX();i++){
        assert(approx(grid[i,4].v.y, +5));
        assert(approx(grid[i,5].v.y, +5));
    }
    
    //Normal components copy
    fill_2D(grid);
    for(int i=0;i<grid.getSizeX();i++) grid[i,0].v.y = +5.0;
    Outflow::Gated(Y_negative).apply(grid);
    for(int i=0;i<grid.getSizeX();i++){
        PrimitiveState w=grid[i,-1]; w.v.y += grid[i,0].v.y;
        expect_close(w, grid[i,0]);
    }
}
//MARK: Gated Outflow - 3D
void DRAGON_Test::verify_boundary_outflow_3D_gated_X(){
    Grid3D grid(3,4, 5,1.0,1.0,1.0, 2);
    
    
    //X- blocks inflow
    fill_3D(grid);
    for(int j=0;j<grid.getSizeY();j++){
        for(int k=0;k<grid.getSizeZ();k++) grid[0,j,k].v.x = +5.0;
    }
    Outflow::Gated(X_negative).apply(grid);
    for(int j=0;j<grid.getSizeY();j++){
        for(int k=0;k<grid.getSizeZ();k++){
            assert(approx(grid[-1,j,k].v.x, 0.0));
            assert(approx(grid[-2,j,k].v.x, 0.0));
            PrimitiveState blocked = grid[-1,j,k]; blocked.v.x = grid[0,j,k].v.x;
            expect_close(blocked, grid[0,j,k]);
        }
    }
    //X- allows outflow
    fill_3D(grid);
    for(int j=0;j<grid.getSizeY();j++){
        for(int k=0;k<grid.getSizeZ();k++) grid[0,j,k].v.x = -5.0;
    }
    Outflow::Gated("X-").apply(grid);
    for(int j=0;j<grid.getSizeY();j++){
        for(int k=0;k<grid.getSizeZ();k++){
            assert(approx(grid[-1,j,k].v.x, -5));
            assert(approx(grid[-2,j,k].v.x, -5));
        }
    }
    
    //X+ blocks inflow
    fill_3D(grid);
    for(int j=0;j<grid.getSizeY();j++){
        for(int k=0;k<grid.getSizeZ();k++) grid[2,j,k].v.x = -5.0;
    }
    Outflow::Gated(X_positive).apply(grid);
    for(int j=0;j<grid.getSizeY();j++){
        for(int k=0;k<grid.getSizeZ();k++){
            assert(approx(grid[3,j,k].v.x, 0.0));
            assert(approx(grid[4,j,k].v.x, 0.0));
        }
    }
    //X+ allows outflow
    fill_3D(grid);
    for(int j=0;j<grid.getSizeY();j++){
        for(int k=0;k<grid.getSizeZ();k++) grid[2,j,k].v.x = +5.0;
    }
    Outflow::Gated("X+").apply(grid);
    for(int j=0;j<grid.getSizeY();j++){
        for(int k=0;k<grid.getSizeZ();k++){
            assert(approx(grid[3,j,k].v.x, +5));
            assert(approx(grid[4,j,k].v.x, +5));
        }
    }
    
    //X - Normal components copy
    fill_3D(grid);
    for(int j=0;j<grid.getSizeY();j++) {
        for(int k=0;k<grid.getSizeZ();k++) grid[0,j,k].v.x = +5.0;
    }
    Outflow::Gated(X_negative).apply(grid);
    for(int j=0;j<grid.getSizeY();j++){
        for(int k=0;k<grid.getSizeZ();k++){
            PrimitiveState w=grid[-1,j,k]; w.v.x += grid[0,j,k].v.x;
            expect_close(w, grid[0,j,k]);
        }
    }
    
}

void DRAGON_Test::verify_boundary_outflow_3D_gated_Y(){
    Grid3D grid(3,4, 5,1.0,1.0,1.0, 2);
    
    //Y- blocks inflow
    fill_3D(grid);
    for(int i=0;i<grid.getSizeX();i++){
        for(int k=0;k<grid.getSizeZ();k++) grid[i,0,k].v.y = +5.0;
    }
    Outflow::Gated(Y_negative).apply(grid);
    for(int i=0;i<grid.getSizeX();i++){
        for(int k=0;k<grid.getSizeZ();k++){
            assert(approx(grid[i,-1,k].v.y, 0.0));
            assert(approx(grid[i,-2,k].v.y, 0.0));
            PrimitiveState blocked = grid[i,-1,k]; blocked.v.y = grid[i,0,k].v.y;
            expect_close(blocked, grid[i,0,k]);
        }
    }
    //Y- allows outflow
    fill_3D(grid);
    for(int i=0;i<grid.getSizeX();i++){
        for(int k=0;k<grid.getSizeZ();k++) grid[i,0,k].v.y = -5.0;
    }
    Outflow::Gated(Y_negative).apply(grid);
    for(int i=0;i<grid.getSizeX();i++){
        for(int k=0;k<grid.getSizeZ();k++){
            assert(approx(grid[i,-1,k].v.y, -5.0));
            assert(approx(grid[i,-2,k].v.y, -5.0));
        }
    }
    
    //Y+ blocks inflow
    fill_3D(grid);
    for(int i=0;i<grid.getSizeX();i++) {
        for(int k=0;k<grid.getSizeZ();k++) grid[i,3,k].v.y = -5.0;
    }
    Outflow::Gated(Y_positive).apply(grid);
    for(int i=0;i<grid.getSizeX();i++){
        for(int k=0;k<grid.getSizeZ();k++){
            assert(approx(grid[i,4,k].v.y, 0.0));
            assert(approx(grid[i,5,k].v.y, 0.0));
        }
    }
    //Y+ allows outflow
    fill_3D(grid);
    for(int i=0;i<grid.getSizeX();i++) {
        for(int k=0;k<grid.getSizeZ();k++) grid[i,3,k].v.y = +5.0;
    }
    Outflow::Gated(Y_positive).apply(grid);
    for(int i=0;i<grid.getSizeX();i++){
        for(int k=0;k<grid.getSizeZ();k++){
            assert(approx(grid[i,4,k].v.y, +5));
            assert(approx(grid[i,5,k].v.y, +5));
        }
    }
    
    //Normal components copy
    fill_3D(grid);
    for(int i=0;i<grid.getSizeX();i++){
        for(int k=0;k<grid.getSizeZ();k++) grid[i,0,k].v.y = +5.0;
    }
    Outflow::Gated(Y_negative).apply(grid);
    for(int i=0;i<grid.getSizeX();i++){
        for(int k=0;k<grid.getSizeZ();k++){
            PrimitiveState w=grid[i,-1,k]; w.v.y += grid[i,0,k].v.y;
            expect_close(w, grid[i,0,k]);
        }
    }
    
}

void DRAGON_Test::verify_boundary_outflow_3D_gated_Z(){
    Grid3D grid(3,4, 5,1.0,1.0,1.0, 2);
    //Z- blocks inflow
    fill_3D(grid);
    for(int i=0;i<grid.getSizeX();i++){
        for(int j=0;j<grid.getSizeY();j++) grid[i,j,0].v.z = +5.0;
    }
    Outflow::Gated("Z-").apply(grid);
    for(int i=0;i<grid.getSizeX();i++){
        for(int j=0;j<grid.getSizeY();j++){
            assert(approx(grid[i,j,-1].v.z, 0.0));
            assert(approx(grid[i,j,-2].v.z, 0.0));
            PrimitiveState blocked = grid[i,j,-1]; blocked.v.z = grid[i,j,0].v.z;
            expect_close(blocked, grid[i,j,0]);
        }
    }
    //Z- allows outflow
    fill_3D(grid);
    for(int i=0;i<grid.getSizeX();i++){
        for(int j=0;j<grid.getSizeY();j++) grid[i,j,0].v.z = -5.0;
    }
    Outflow::Gated(Z_negative).apply(grid);
    for(int i=0;i<grid.getSizeX();i++){
        for(int j=0;j<grid.getSizeY();j++){
            assert(approx(grid[i,j,-1].v.z, -5.0));
            assert(approx(grid[i,j,-2].v.z, -5.0));
        }
    }
    //Z+ blocks inflow
    fill_3D(grid);
    for(int i=0;i<grid.getSizeX();i++) {
        for(int j=0;j<grid.getSizeY();j++) grid[i,j,4].v.z = -5.0;
    }
    Outflow::Gated("Z+").apply(grid);
    for(int i=0;i<grid.getSizeX();i++){
        for(int j=0;j<grid.getSizeY();j++){
            assert(approx(grid[i,j,5].v.z, 0.0));
            assert(approx(grid[i,j,6].v.z, 0.0));
        }
    }
    //Z+ allows outflow
    fill_3D(grid);
    for(int i=0;i<grid.getSizeX();i++) {
        for(int j=0;j<grid.getSizeY();j++) grid[i,j,4].v.z = +5.0;
    }
    Outflow::Gated(Z_positive).apply(grid);
    for(int i=0;i<grid.getSizeX();i++){
        for(int j=0;j<grid.getSizeY();j++){
            assert(approx(grid[i,j,5].v.z, +5));
            assert(approx(grid[i,j,6].v.z, +5));
        }
    }
    
    //Normal components copy
    fill_3D(grid);
    for(int i=0;i<grid.getSizeX();i++){
        for(int j=0;j<grid.getSizeY();j++) grid[i,j,0].v.z = +5.0;
    }
    Outflow::Gated(Z_negative).apply(grid);
    for(int i=0;i<grid.getSizeX();i++){
        for(int j=0;j<grid.getSizeY();j++){
            PrimitiveState w=grid[i,j,-1]; w.v.z += grid[i,j,0].v.z;
            expect_close(w, grid[i,j,0]);
        }
    }
}
