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

typedef DistGrid3D MyGrid;
constexpr double rho_jet = 3;
constexpr double rho_amb = rho_jet/6.0;
constexpr double p_amb = rho_amb * 1.5e-5 ;

Grid& Problem::makeProblem(){
    int res = 1;
    double dx = 0.1/res, dy = 0.1/res, dz = 0.1/res;
    auto _grid = new MyGrid(64*res, 64*res,64*res, dx, dy, dz);
    MyGrid& grid = *_grid;
    
    grid.boundary = Boundary::ToroidalJet(rho_jet, 7e-2, p_amb, 10.0, 1.0, 2.4, "Z");
    
    for(int i=0; i<=grid.getSizeX(); i++){
        for(int j=0; j<=grid.getSizeY(); j++){
            for(int k=0; k<=grid.getSizeZ(); k++){
                grid[i,j,k].rho = rho_amb;
                grid[i,j,k].p = p_amb;
                grid[i,j,k].v = {0,0,0};
                grid.A()[i,j,k] = {0,0,0};
            }
        }
    }
    grid.initialize_B_fields();
        
    return grid;
}


void Problem::cycleComplete(Grid& problem, int cycle, double time){
    MyGrid& grid = *dynamic_cast<MyGrid*>(&problem);
    
    //Frame Number
    auto numStr = std::string(cycle<1000 ? "0" : "") + std::string(cycle<100 ? "0" : "") + std::string(cycle<10 ? "0" : "") +  std::to_string(cycle);
    std::cout<<"Frame "<<numStr<<" computed, ";
    //Time
    int h = floor(time/3600.0), m = floor((time-h*3600.0)/60.0);
    double s = round((time - h*3600 - m*60)*100)/100.0;
    std::cout << "Time: "<< h << "h "<< m <<"m " << s <<"s\n";
    //Output
    std::string filename = "/Users/bobbiemarkwick/DRAGON_OUT/frame-" + numStr  +".csv";
    std::ofstream out;
    out.open (filename);
    
    for(int k=0; k<grid.getSizeZ(); k++){
        for(int i=0; i<grid.getSizeX(); i++){
            out<< grid[i,grid.getSizeY()/2,k].rho  * (3/rho_jet) << (i+1 == grid.getSizeX() ? "\n" : ",");
        }
    }
    out.close();
}
