//
//  AMRGrid.hpp
//  DRAGON
//
//  Created by Bobbie Markwick on 12/06/2026.
//

#ifndef AMRGrid_hpp
#define AMRGrid_hpp

#include "Grid.hpp"
#include "Config.h"
#include <vector>
#include <memory>


struct AMRGrid2D{
    
    AMRGrid2D(int nx, int ny, double dx, double dy, int ghosts=3);
    ~AMRGrid2D();
    AMRGrid2D(const Grid2D&) = delete; //No copying
    AMRGrid2D& operator=(const AMRGrid2D&) = delete;

    
    //Grid access
    PrimitiveState& operator[](int i,int j);
    const PrimitiveState& operator[](int i,int j) const;
    int getSizeX() const, getSizeY() const, getGhosts() const;

    //Boundary
    Boundary::BoundaryList boundary = Boundary::BoundaryList();
    
    //Advance Forward in time
    void advance(double dt, bool check_cfl = true);
    
private:
    Grid2D data;
    std::vector<std::unique_ptr<AMRGrid2D>> children;
    int ghosts, ncx, ncy;
    double size_x, size_y, dx, dy;
    
    void loadFromChildren();
    void pushToChildren();
    
    void debugCheckChildGhosts();

};


#endif /* AMRGrid_hpp */
