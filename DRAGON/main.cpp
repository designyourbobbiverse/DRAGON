//  DRAGON
//  Created by Bobbie Markwick
//

#include <iostream>
#include <math.h>
#include <fstream>
#include <chrono>

#include "AMRGrid.hpp"

int main(int argc, const char * argv[]) {    
    double dx = 0.01, dy = 0.01;
    AMRGrid2D grid(1000, 1000, dx, dy);
    grid.boundary = Boundary::Reflective();
    
    const double rho0 = 1.0;
    const double p_ambient = 0.1;
    const double p_blast   = 10.0;
    const double r0 = 0.3;
    
    auto start = std::chrono::system_clock::now();
    
    for(int i=0; i<grid.getSizeX(); i++){
        for(int j=0; j<grid.getSizeY(); j++){
            grid[i,j].rho = rho0;
            grid[i,j].vx = 0; grid[i,j].vy = 0; grid[i,j].vz = 0;
            
            double x = (i + 0.5 - grid.getSizeX()/2) * dx;
            double y = (j + 0.5 - grid.getSizeY()/2) * dy;
            double r = sqrt(x*x + y*y);
            grid[i,j].p = (r < r0) ? p_blast : p_ambient;
        }
    }
   
    for(int n = 0; n<10000; n++){
        if(n>0) grid.advance(0.01);
        
        auto numStr = std::string(n<1000 ? "0" : "") + std::string(n<100 ? "0" : "") + std::string(n<10 ? "0" : "") +  std::to_string(n);
        std::cout<<"Frame "<<numStr<<" computed, ";
        //Time
        auto now = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = now - start;
        int h = floor(elapsed_seconds.count()/3600.0), m = floor((elapsed_seconds.count()-h*3600.0)/60.0);
        double s = round((elapsed_seconds.count() - h*3600 - m*60)*100)/100.0;
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
    
    
    
    return 0;
}
