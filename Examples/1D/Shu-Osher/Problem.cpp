//
//  Problem_Shu_Osher.cpp
//  DRAGON/Examples/1D/Shu-Osher
//
//  Created by Bobbie Markwick on 7/07/2026.
//

#include "Problem.hpp"
#include "DistGrid.hpp"
#include <cmath>

typedef DistGrid1D MyGrid;//Choose the dimension of your grid here

constexpr double rho_L = 3.857143;
constexpr double rho_R = 1.0;
constexpr double p_L = 10.333333;
constexpr double p_R = 1.0;
constexpr double vxL = 2.629369;
constexpr int nx = 1280*8;
constexpr double dx = 10.0/nx;

Grid& Problem::makeProblem(){
    auto grid = new MyGrid(nx, dx);
    grid->boundary = Boundary::Outflow();
    return *grid;
}

void Problem::initializeProblem(Grid& problem){
    MyGrid& grid = *dynamic_cast<MyGrid*>(&problem);
    
    for(int i = 0; i<nx; i++){
        grid[i].rho = i < nx/10 ? rho_L : rho_R;
        grid[i].p = i < nx/10 ? p_L : p_R;
        grid[i].v = i < nx/10 ? vec3{vxL,0,0} : vec3{0,0,0};
        
        if(i > nx/10) grid[i].rho += 0.2 * sin(5.0 * (i-nx/2)*dx);
    }
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
