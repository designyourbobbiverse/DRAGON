//
//  Boundary/GhostFill.hpp
//  DRAGON
//
//  Created by Bobbie Markwick on 28/05/2026.
//

#ifndef GridLink_hpp
#define GridLink_hpp

#include "GhostFill.hpp"
#include "Grid.hpp"
#include <vector>

using namespace Boundary;

class GridLink: public GhostFill {
public:
    GridLink(Grid1D* grid, int face);
    GridLink(Grid2D* grid, int face);
    GridLink(Grid3D* grid, int face);
    ~GridLink() = default;
    
    void apply(Grid1D& _target) override;
    void apply(Grid2D& _target) override;
    void apply(Grid3D& _target) override;
    
    void reload();
    
private:
    Grid1D* ref1D = nullptr;
    std::vector<PrimitiveState> ghosts1D;
    void reload1D();

    Grid2D* ref2D = nullptr;
    std::vector<std::vector<PrimitiveState>> ghosts2D;
    void reload2D();

    Grid3D* ref3D = nullptr;
    std::vector<std::vector<std::vector<PrimitiveState>>> ghosts3D;
    void reload3D();

};





#endif
