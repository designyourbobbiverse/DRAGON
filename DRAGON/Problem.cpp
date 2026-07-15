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
    auto _grid = new MyGrid(100, 100,100, 0.1, 0.1, 0.1);
    MyGrid& grid = *_grid;
    
    grid.boundary = Boundary::Periodic();
    
    return grid;
}

void Problem::initializeProblem(Grid& problem){
    MyGrid& grid = *dynamic_cast<MyGrid*>(&problem);
    //This is where you should initialize the initial data of the grid

    PrimitiveState W{};
    W.rho = 1.0;
    W.v = {5,0,5};
    W.p = 50;

    for (int i = 0; i <= grid.getSizeX(); i++){
        for (int j = 0; j <= grid.getSizeY(); j++){
            for (int k = 0; k <= grid.getSizeZ(); k++){
                grid[i,j,k] = W;
                double x = i*0.1-5, y = j*0.1 - 5, z = k*0.1 - 5;
                double r2 = x*x + y*y + z*z;
                double Az = fmax(4.0-r2, 0);
                grid._A()[i,j,k] = vec3{0,0,Az};
            }
        }
    }
    grid.initialize_B_fields();
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
