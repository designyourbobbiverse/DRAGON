//
//  ParentGridFill.hpp
//  DRAGON
//
//  Created by Bobbie Markwick on 12/06/2026.
//

#ifndef ParentGridFill_hpp
#define ParentGridFill_hpp

#include "GhostFill.hpp"
#include "AMRGrid.hpp"


class ParentGridFill : public GhostFill {
public:
    ParentGridFill(AMRGrid2D* parent, int x, int y);
    
    void apply(Grid1D& grid) override;
    void apply(Grid2D& grid) override;
    void apply(Grid3D& grid) override;
private:
    AMRGrid2D* parent;
    int x_offset, y_offset;
};


#endif /* ParentGridFill_hpp */
