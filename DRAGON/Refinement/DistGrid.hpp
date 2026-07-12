//
//  DistGrid.hpp
//  DRAGON/Refinement
//  User-Facing Header file
//
//  Created by Bobbie Markwick on 12/06/2026.
//

#ifndef AMRGrid_hpp
#define AMRGrid_hpp

#include "Grid.hpp"
#include "Config.h"
#include <vector>
#include <memory>

template <typename T>
struct DistGrid{
protected:
    std::vector<std::unique_ptr<T>> children;//The subgrids that actually do the processing
    //Synchronization with children
    virtual void loadFromChildren() = 0;
    virtual void pushToChildren() = 0;
    //Advance Forward in time
    void step(double dt);
};



struct DistGrid1D: public Grid1D, public DistGrid<DistGrid1D>{
    
    DistGrid1D(int nx, double dx, int ghosts=3, bool root=true);
    DistGrid1D(const DistGrid1D&) = delete; //No copying
    DistGrid1D& operator=(const DistGrid1D&) = delete;

    //Advance Forward in time
    void split_step(double dt) override;
    void unsplit_step(double dt) override;
    bool on_step_fail(const std::exception& e) override;
private:
    int ncx; //The number of children we have
    double size_x; //Total width
    //Synchronization with children
    void loadFromChildren() override;
    void pushToChildren() override;
};


struct DistGrid2D: public Grid2D, public DistGrid<DistGrid2D>{
    
    DistGrid2D(int nx, int ny, double dx, double dy, int ghosts=3, bool root=true);
    DistGrid2D(const DistGrid2D&) = delete; //No copying
    DistGrid2D& operator=(const DistGrid2D&) = delete;

    //Advance Forward in time
    void split_step(double dt) override;
    void unsplit_step(double dt) override;
    bool on_step_fail(const std::exception& e) override;
private:
    int ncx, ncy;//The number of children we have in each dimension
    double size_x, size_y;//Total physical size of the Grid
    //Synchronization with Children
    void loadFromChildren() override;
    void pushToChildren() override;
};


struct DistGrid3D: public Grid3D, public DistGrid<DistGrid3D>{
    
    DistGrid3D(int nx, int ny, int nz, double dx, double dy, double dz, int ghosts=3, bool root=true);
    DistGrid3D(const DistGrid3D&) = delete; //No copying
    DistGrid3D& operator=(const DistGrid3D&) = delete;
    
    
    //Advance Forward in time
    void split_step(double dt) override;
    void unsplit_step(double dt) override;
    bool on_step_fail(const std::exception& e) override;
private:
    int ncx, ncy, ncz;//The number of children we have in each dimension
    double size_x, size_y, size_z; //The total width of the grid

    //Synchronising with children
    void loadFromChildren() override;
    void pushToChildren() override;
};


#endif /* AMRGrid_hpp */
