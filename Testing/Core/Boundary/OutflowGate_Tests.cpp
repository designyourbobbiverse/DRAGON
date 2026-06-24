
//
//  OutflowGate_Tests.cpp
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

static PrimitiveState G = make_tagged_state(-666);

void fill_1D(Grid1D& grid);
void fill_2D(Grid2D& grid);
void fill_3D(Grid3D& grid);

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
