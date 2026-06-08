//
//  Boundary/GhostFill.hpp
//  DRAGON
//
//  Created by Bobbie Markwick on 28/05/2026.
//

#ifndef GhostFill_hpp
#define GhostFill_hpp

class Grid1D;
class Grid2D;
class Grid3D;

class GhostFill {
public:
    virtual ~GhostFill() = default;
    int get_faces() const;
    
    virtual void apply(Grid1D& grid) = 0;
    virtual void apply(Grid2D& grid) = 0;
    virtual void apply(Grid3D& grid) = 0;
protected:
    int faces;
    bool corners;
    explicit GhostFill(int faces_, bool corner_ghosts) : faces(faces_), corners(corner_ghosts) {}
};





#endif
