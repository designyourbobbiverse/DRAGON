//  DRAGON
//  Created by Bobbie Markwick
//

#include <iostream>
#include <math.h>
#include <fstream>

#include "Grid.hpp"

int main(int argc, const char * argv[]) {
    
    Grid2D grid(100, 100, 0.1, 0.1);
    grid.boundary = Boundary::Reflective();
    
    const double rho0 = 1.0;
    const double p_ambient = 0.1;
    const double p_blast   = 10.0;
    const double r0 = 1;
    
    for(int i=0; i<grid.getSizeX(); i++){
        for(int j=0; j<grid.getSizeY(); j++){
            grid[i,j].rho = rho0;
            grid[i,j].vx = 0; grid[i,j].vy = 0; grid[i,j].vz = 0;
            
            double x = (i + 0.5) * 0.1 - 5;
            double y = (j + 0.5) * 0.1 - 5;
            double r = sqrt(x*x + y*y);
            grid[i,j].p = (r < r0) ? p_blast : p_ambient;
        }
    }
   
    for(int n = 0; n<1000; n++){
        if(n>0) grid.advance(0.01);
        
        auto numStr = std::string(n<100 ? "0" : "") + std::string(n<10 ? "0" : "") +  std::to_string(n);
        std::cout<<"Frame "<<numStr<<" computed\n";
    
        std::string filename = "/Users/bobbiemarkwick/DRAGON_OUT/frame-" + numStr  +".csv";
        std::ofstream out;
        out.open (filename);
        
        for(int i=0; i<grid.getSizeX(); i++){
            for(int j=0; j<grid.getSizeY(); j++){
                out<< grid[i,j].rho << (j+1 == grid.getSizeY() ? "\n" : ",");
            }
        }
        out.close();
    }
    
    
    
    return 0;
}
