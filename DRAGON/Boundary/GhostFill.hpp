//
//  GhostFill.hpp
//  DRAGON/Boundary
//
//  Created by Bobbie Markwick on 28/05/2026.
//

#ifndef GhostFill_hpp
#define GhostFill_hpp

class Grid;
class Grid1D;
class Grid2D;
class Grid3D;

class GhostFill {
public:
    virtual ~GhostFill() = default;
    int get_faces() const;
    
    //Fill the appropriate Ghost Cells within a given grid
    virtual void apply(Grid1D& grid) = 0;
    virtual void apply(Grid2D& grid) = 0;
    virtual void apply(Grid3D& grid) = 0;
    void apply(Grid& grid);
protected:
    int faces; //A bit mask desgnating to which faces this GhostFill shall be applied
    bool corners; //Designates whether this object should also fill corner cells of the respective face
    explicit GhostFill(int faces_, bool corner_ghosts) : faces(faces_), corners(corner_ghosts) {}
};





#endif
