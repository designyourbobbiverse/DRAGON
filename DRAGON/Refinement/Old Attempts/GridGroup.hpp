//
//  GridGroup.hpp
//  DRAGON
//
//  Created by Bobbie Markwick on 11/06/2026.
//

#ifndef GridGroup_hpp
#define GridGroup_hpp

#include "Grid.hpp"
#include "GridLink.hpp"
#include <vector>

struct GridGroup2D{

    GridGroup2D(int nx, int ny, double dx, double dy, int ghosts=1);
    GridGroup2D(const GridGroup2D&) = delete; //No copying
    GridGroup2D& operator=(const GridGroup2D&) = delete;
    
    //Grid access
    PrimitiveState& operator[](int i,int j);
    const PrimitiveState& operator[](int i,int j) const;
    int getSizeX() const, getSizeY() const, getGhosts() const;
    
    //Boundary
    
    //Advance Forward in time
    void advance(double dt);
private:
    std::array<std::unique_ptr<Grid2D>,16> zones;
    void reloadLinks();
};


#endif /* GridGroup_hpp */
