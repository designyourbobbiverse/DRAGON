//
//  Boundary/Boundary.hpp
//  DRAGON/Boundary
//  User-Facing Header file
//
//  Created by Bobbie Markwick on 28/05/2026.
//

#ifndef Boundary_hpp
#define Boundary_hpp

#include "GhostFill.hpp"

#include <string> //For face initialization
#include "FluidElement/FluidElement.hpp" //For Boundary::fixed

namespace DRAGON::Boundary{

//MARK: Boundary Composition
//The conditions will be applied in order, with later conditions overriding prior conditions for overlapping cells
//Outflow(missing_faces) will be added at the beginning if any faces are missing.
class BoundaryList; //See BoundaryList.hpp
//Supported Addition and Assignment operators
//List = Element + Element (Element is any GhostFill that isn't also another List)
//List = List + Element
//List += Element
//List = Element + List
//List = List + List
//List += List
//List = List
//List = Element


//MARK: Faces
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
//X = X, Y = Y, Z = Z (case insensitive)
//- after X Y or Z means only the negative-side boundary
//+ after X Y or Z means only the positive-side boundary
//Example: XYZ+ is everything except negative Z
int face_mask(std::string s);


//MARK: Predefined types

//Outflow boundary condition
//Sets each boundary cell to be equal to the nearest interior cell
class Outflow : public GhostFill {
public:
    Outflow(std::string faces, bool corner_ghosts=true, bool gated=false);
    Outflow(int faces = X|Y|Z, bool corner_ghosts=true, bool gated=false);
    
    //Variant which overrides any inflow component of velocity to be zero
    static Outflow Gated(std::string faces, bool corner_ghosts=true);
    static Outflow Gated(int faces = X|Y|Z, bool corner_ghosts=true);
    
    void apply(Grid1D& grid) override;
    void apply(Grid2D& grid) override;
    void apply(Grid3D& grid) override;
    
    bool gated; //If True, overrides any inflow component of velocity to be zero
};

//Fixed Boundary condition
//Sets each boundary cell to be equal to some specified state
class Fixed : public GhostFill {
public:
    Fixed(PrimitiveState state, std::string faces, bool corner_ghosts=true);
    Fixed(PrimitiveState state, int faces = X|Y|Z, bool corner_ghosts=true);
    
    void apply(Grid1D& grid) override;
    void apply(Grid2D& grid) override;
    void apply(Grid3D& grid) override;
protected:
    PrimitiveState state;
};

//Reflective Boundary Condition
//Sets the boundary cells to be a reflection of the physical cells
//The velocity component normal to the boundary face will be flipped, transverse velocities preserved
class Reflective : public GhostFill {
public:
    Reflective(std::string faces, bool corner_ghosts = true);
    Reflective(int faces = X|Y|Z, bool corner_ghosts = true);

    void apply(Grid1D& grid) override;
    void apply(Grid2D& grid) override;
    void apply(Grid3D& grid) override;
};

//Periodic Boundary Condition
//Sets the boundary cells to match the interior cells at the opposite side of the grid
class Periodic : public GhostFill {
public:
    //Including either X- or X+ in the face mask will automatically include the other, and similarly for other dimensions
    Periodic(std::string faces, bool corner_ghosts=true);
    Periodic(int faces = X|Y|Z, bool corner_ghosts=true);
    
    void apply(Grid1D& grid) override;
    void apply(Grid2D& grid) override;
    void apply(Grid3D& grid) override;
};

//Ignore Boundary Condition
//A placeholder which itself does nothing
//Useful in a BoundaryList if you need to manually fill cells, since any faces specified as part of an Ignore won't be subjected to an implicit Outflow
class Ignore : public GhostFill {
public:
    Ignore(std::string faces, bool corner_ghosts = true);
    Ignore(int faces = X|Y|Z, bool corner_ghosts = true);
    
    void apply(Grid1D& grid) override;
    void apply(Grid2D& grid) override;
    void apply(Grid3D& grid) override;
};

}

//Keep the implmentation of BoundaryList in a separate file
#define BOUNDARY_LIST_HPP
#include "BoundaryList.hpp"


#endif
