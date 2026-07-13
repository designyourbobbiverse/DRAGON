//
//  Problem.cpp
//  DRAGON
//
//  Created by Bobbie Markwick on 24/06/2026.
//

#include "Problem.hpp"
#include "DistGrid.hpp"
#include "Jets.hpp"
#include <cmath>
#include <iostream>
#include <fstream>

typedef DistGrid3D MyGrid;
constexpr double rho_jet = 3;
constexpr double rho_amb = rho_jet/6.0;
constexpr double p_amb = rho_amb * 1.5e-5 ;

Grid& Problem::makeProblem(){
    int res = 2;
    double dx = 0.1/res, dy = 0.1/res, dz = 0.1/res;
    auto _grid = new MyGrid(64*res, 64*res,64*res, dx, dy, dz);
    MyGrid& grid = *_grid;
    
    grid.boundary = MHDJet(rho_jet, 7e-2, p_amb, 5.0, 1.0, 2.4, "Z");
    
    for(int i=0; i<=grid.getSizeX(); i++){
        for(int j=0; j<=grid.getSizeY(); j++){
            for(int k=0; k<=grid.getSizeZ(); k++){
                grid[i,j,k].rho = rho_amb;
                grid[i,j,k].p = p_amb;
                grid[i,j,k].v = {0,0,0};
                grid._A()[i,j,k] = {0,0,0};
            }
        }
    }
    grid.initialize_B_fields();
        
    return grid;
}


void Problem::beforeCycle(Grid &problem, int cycle, double t){
    MyGrid& grid = *dynamic_cast<MyGrid*>(&problem);
    
    //Do any processing you need to do between cycles
    //This is called before the Godunov scheme. It isn't called for the initial frame
    
}


void Problem::afterCycle(Grid &problem, int cycle, double t){
    MyGrid& grid = *dynamic_cast<MyGrid*>(&problem);
    
    //Do any processing you need to do between cycles
    //This is called after the Godunov scheme but before the file write. It isn't called for the initial frame
    
}

void Problem::problemComplete(Grid& problem, double t){
    MyGrid& grid = *dynamic_cast<MyGrid*>(&problem);

    //This is called only after the final time is reached.
}
