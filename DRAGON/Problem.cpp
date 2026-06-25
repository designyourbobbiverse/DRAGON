//
//  Problem.cpp
//  DRAGON
//
//  Created by Bobbie Markwick on 24/06/2026.
//

#include "main/Problem.hpp"
#include "Refinement/DistGrid.hpp"
#include <cmath>
#include <iostream>
#include <fstream>
#include <chrono>


const double rho0 = 1.0;
const double p_ambient = 0.1;
const double p_blast   = 10.0;
const double r0 = 0.3;

Advanceable& Problem::makeProblem(){
    double dx = 0.01, dy = 0.01;
    auto _grid = new DistGrid2D(1000, 1000, dx, dy);
    DistGrid2D& grid = *_grid;
    
    grid.boundary = Boundary::Reflective();
    
    for(int i=0; i<grid.getSizeX(); i++){
        for(int j=0; j<grid.getSizeY(); j++){
            grid[i,j].rho = rho0;
            grid[i,j].v.x = 0; grid[i,j].v.y = 0; grid[i,j].v.z = 0;
            
            double x = (i + 0.5 - grid.getSizeX()/2) * dx;
            double y = (j + 0.5 - grid.getSizeY()/2) * dy;
            double r = sqrt(x*x + y*y);
            grid[i,j].p = (r < r0) ? p_blast : p_ambient;
        }
    }
    
    return grid;
}


void Problem::cycleComplete(Advanceable& problem, int cycle, double time){
    DistGrid2D& grid = *dynamic_cast<DistGrid2D*>(&problem);
    
    //Frame Number
    auto numStr = std::string(cycle<1000 ? "0" : "") + std::string(cycle<100 ? "0" : "") + std::string(cycle<10 ? "0" : "") +  std::to_string(cycle);
    std::cout<<"Frame "<<numStr<<" computed, ";
    //Time
    int h = floor(time/3600.0), m = floor((time-h*3600.0)/60.0);
    double s = round((time - h*3600 - m*60)*100)/100.0;
    std::cout << "Time: "<< h << "h "<< m <<"m " << s <<"s, ";
    //Output
    std::string filename = "/Users/bobbiemarkwick/DRAGON_OUT/frame-" + numStr  +".csv";
    std::ofstream out;
    out.open (filename);
    
    double rho_min = rho0,rho_max = 0;
    for(int i=0; i<grid.getSizeX(); i++){
        for(int j=0; j<grid.getSizeY(); j++){
            out<< grid[i,j].rho << (j+1 == grid.getSizeY() ? "\n" : ",");
            rho_min = std::min(rho_min, grid[i,j].rho);
            rho_max = std::max(rho_max, grid[i,j].rho);
        }
    }
    std::cout<< "rho range: ("<<rho_min<<","<<rho_max<<")\n";
    
    out.close();
}
