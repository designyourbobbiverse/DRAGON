//
//  Boundary/Boundary.hpp
//  DRAGON
//
//  Created by Bobbie Markwick on 28/05/2026.
//

#ifndef Boundary_hpp
#define Boundary_hpp

#include "Godunov.hpp"
#include <vector>
#include <memory>
#include <concepts>
#include <type_traits>
#include <utility>


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


//Intialization from a string
//X = X, Y = Y, Z = Z
//- after XYZ means only the negative, + after XYZ means only the positive
//Example: XYZ+ is everyting except negative Z
int face_mask(std::string s);

class Fixed : public GhostFill {
public:
    Fixed(PrimitiveState state, std::string faces, bool corner_ghosts=true);
    Fixed(PrimitiveState state, int faces = X|Y|Z, bool corner_ghosts=true);
    ~Fixed() = default;
    
    void apply(Grid1D& grid) override;
    void apply(Grid2D& grid) override;
    void apply(Grid3D& grid) override;
protected:
    PrimitiveState state;
};

class Outflow : public GhostFill {
public:
   
    Outflow(std::string faces, bool corner_ghosts=true, bool gated=false);
    Outflow(int faces = X|Y|Z, bool corner_ghosts=true, bool gated=false);
    static Outflow Gated(std::string faces, bool corner_ghosts=true);
    static Outflow Gated(int faces = X|Y|Z, bool corner_ghosts=true);
    
    void apply(Grid1D& grid) override;
    void apply(Grid2D& grid) override;
    void apply(Grid3D& grid) override;

    bool gated; //Set to true if you want to ensure no accidental inflow occurs
};

class Reflective : public GhostFill {
public:
    Reflective(std::string faces, bool corner_ghosts = true);
    Reflective(int faces = X|Y|Z, bool corner_ghosts = true);

    
    
    void apply(Grid1D& grid) override;
    void apply(Grid2D& grid) override;
    void apply(Grid3D& grid) override;
};


class Periodic : public GhostFill {
public:
    Periodic(std::string faces, bool corner_ghosts=true);
    Periodic(int faces = X|Y|Z, bool corner_ghosts=true);
    
    void apply(Grid1D& grid) override;
    void apply(Grid2D& grid) override;
    void apply(Grid3D& grid) override;
};



//A collection of boundary conditions
//The conditions will be applied in order, with later conditions overriding prior conditions for overlapping cells
//When initializing, Outflow(missing_faces) will be added at the beginning if any faces are missing.
class Boundaries : public GhostFill {
public:
    template<typename... Bs> Boundaries(Bs&&... bs);
    
    void apply(Grid1D& grid) override;
    void apply(Grid2D& grid) override;
    void apply(Grid3D& grid) override;
private:
    std::vector<std::unique_ptr<GhostFill>>  boundaries;
    template<typename B> void addBoundary(B&& b);
    void resetImplicit();
    
    Outflow implicits;
    bool stale = true;
    
};

}

#include "Boundaries.tpp"


#endif
