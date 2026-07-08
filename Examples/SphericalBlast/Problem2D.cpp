//
//  Problem2D.cpp
//  DRAGON/Examples/SphericalBlast
//
//  Created by Bobbie Markwick on 24/06/2026.
//

#include "Problem.hpp"
#include "DistGrid.hpp"
#include <cmath>
#include "Constants.h"

typedef DistGrid2D MyGrid;//Choose the dimension of your grid here

constexpr double rho = 1.0;
constexpr double p_amb = 1e-5;
constexpr double E_blast = 1.0;
constexpr int n = 512;
constexpr double r0 = 12.0/n;


Grid& Problem::makeProblem(){
    //Construct your grid object. Don't worry about initial setup, you'll do that later
    auto grid = new MyGrid(n,n, 1.0/n,1.0/n);
    grid->boundary = Boundary::Outflow();
    return *grid;
}

double circleCellFraction(double x1, double x2, double y1, double y2, double r){
    const double area = (x2 - x1) * (y2 - y1);

    // Quick accept: all four corners inside / outside
    auto inside = [r](double x, double y) { return x*x + y*y <= r*r; };
    if (inside(x1,y1) && inside(x1,y2) && inside(x2,y1) && inside(x2,y2)) return 1.0;
    if (!inside(x1,y1) && !inside(x1,y2) && !inside(x2,y1) && !inside(x2,y2)) return 0.0;

    // 8-point Gauss-Legendre nodes/weights on [-1,1]
    static constexpr double xs[8] = {
        -0.9602898564975363,
        -0.7966664774136267,
        -0.5255324099163290,
        -0.1834346424956498,
         0.1834346424956498,
         0.5255324099163290,
         0.7966664774136267,
         0.9602898564975363
    };

    static constexpr double ws[8] = {
        0.1012285362903763,
        0.2223810344533745,
        0.3137066458778873,
        0.3626837833783620,
        0.3626837833783620,
        0.3137066458778873,
        0.2223810344533745,
        0.1012285362903763
    };

    double xm = 0.5 * (x1 + x2);
    double xr = 0.5 * (x2 - x1);

    double integral = 0.0;

    for (int q = 0; q < 8; q++) {
        double x = xm + xr * xs[q];

        double vertical_length = 0.0;

        if (std::abs(x) < r) {
            double yCircle = std::sqrt(r*r - x*x);

            double lo = std::max(y1, -yCircle);
            double hi = std::min(y2,  yCircle);

            if (hi > lo) {
                vertical_length = hi - lo;
            }
        }

        integral += ws[q] * vertical_length;
    }

    integral *= xr;

    return integral / area;
}

void Problem::initializeProblem(Grid& problem){
    MyGrid& grid = *dynamic_cast<MyGrid*>(&problem);

    const double p_blast = (_gamma - 1.0) * E_blast / (M_PI * r0 * r0);

    //Set up the ambient grid and calculate the size of the blast
    for(int i=0; i<n;i++){
        for(int j=0; j<n; j++){
            double x = (i + 0.5)/n  - 0.5;
            double y = (j + 0.5)/n - 0.5;
            double r = sqrt(x*x + y*y);

            grid[i,j].rho = rho;
            grid[i,j].v = {0,0,0};
            grid[i,j].p = p_amb;
            
            grid[i,j].p += p_blast * circleCellFraction(x-0.5/n, x+0.5/n, y-0.5/n, y+0.5/n, r0);
        }
    }
    //Inject the blast
    for(int i=0; i<n;i++){
        for(int j=0; j<n; j++){
            double x = (i + 0.5)/n  - 0.5;
            double y = (j + 0.5)/n - 0.5;
            double q = sqrt(x*x + y*y) / r0;

            if(q < 1) grid[i,j].p += p_blast ;
        }
    }
    
    
}


void Problem::cycleComplete(Grid& problem, int cycle){
    MyGrid& grid = *dynamic_cast<MyGrid*>(&problem);
    
    //Do any processing you need to do between cycles
    //This is called before the file write, and is called for the initial frame as well as subsequent updates
    
}
