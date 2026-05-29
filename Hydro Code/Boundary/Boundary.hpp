//
//  Boundary.hpp
//  Hydro Code
//
//  Created by Bobbie Markwick on 28/05/2026.
//

#ifndef Boundary_hpp
#define Boundary_hpp

#include "Godunov.hpp"
#include <vector>

namespace Boundary{

constexpr int X_negative = 1 << 0;
constexpr int X_positive = 1 << 1;
constexpr int X = X_negative | X_positive;

constexpr int Y_negative = 1 << 2;
constexpr int Y_positive = 1 << 3;
constexpr int Y = Y_negative | Y_positive;

constexpr int Z_negative = 1 << 4;
constexpr int Z_positive = 1 << 5;
constexpr int Z = Z_negative | Z_positive;



class BoundaryType {
public:
    explicit BoundaryType(int faces_, bool corner_ghosts) : faces(faces_), corners(corner_ghosts) {}
    virtual ~BoundaryType() = default;
    
    virtual void apply(Grid1D& grid) const = 0;
    virtual void apply(Grid2D& grid) const = 0;
    virtual void apply(Grid3D& grid) const = 0;
protected:
    int faces;
    bool corners;
};


class Fixed : public BoundaryType {
public:
    Fixed(PrimitiveState state);
    Fixed(PrimitiveState state, int faces);
    Fixed(PrimitiveState state, int faces, bool corner_ghosts);
    ~Fixed() = default;
    
    void apply(Grid1D& grid) const override;
    void apply(Grid2D& grid) const override;
    void apply(Grid3D& grid) const override;
protected:
    PrimitiveState state;
};

class Outflow : public BoundaryType {
public:
    Outflow();
    Outflow(int faces);
    Outflow(int faces, bool corner_ghosts);
    Outflow(int faces, bool corner_ghosts,bool gated);
    static Outflow Gated();
    static Outflow Gated(int faces);
    static Outflow Gated(int faces, bool corner_ghosts);
    
    void apply(Grid1D& grid) const override;
    void apply(Grid2D& grid) const override;
    void apply(Grid3D& grid) const override;

    bool gated; //Set to true if you want to ensure no accidental inflow occurs
};

class Reflective : public BoundaryType {
public:
    Reflective();
    Reflective(int faces);
    Reflective(int faces, bool corner_ghosts);

    
    
    void apply(Grid1D& grid) const override;
    void apply(Grid2D& grid) const override;
    void apply(Grid3D& grid) const override;
};

class Periodic : public BoundaryType {
public:
    Periodic();
    Periodic(int faces);
    Periodic(int faces, bool corner_ghosts);
    
    void apply(Grid1D& grid) const override;
    void apply(Grid2D& grid) const override;
    void apply(Grid3D& grid) const override;
};

};

#endif
