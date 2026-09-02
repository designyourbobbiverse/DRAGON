//
//  Hydro/Grid.hpp
//  DRAGON/Hydro
//  User-Facing Header file
//
//  Created by Bobbie Markwick on 12/05/2026.
//

#ifndef Grid_hpp
#define Grid_hpp

#include "Hydro/ExtendedArray/ArrayTypes.hpp"   //For ExtendedArray, PrimitiveState
#include "Boundary/Boundary.hpp"                //For Boundary::BoundaryList
#include "Passives/Passives.hpp"                //For Passives

#include "Config.h" //For enabling or disabling MHD declarations

namespace DRAGON{
class Grid{
public:
    //Advance a single timestep
    void advance_step(double dt); //Calls split_step or unsplit_step in accordance with Config.h
    virtual void split_step(double dt) = 0; //Should advance all dimensnions by dt
    virtual void unsplit_step(double dt) = 0;
    
    //Advance the grid by time dt. If check_cfl is true, the timestep will be CFL-limited and looped if necessary
    void advance(double dt, bool check_cfl=true); //Calls advance_split/advance_unsplit as determined in Config.h
    void advance_split(double dt, bool check_cfl=true); //Calls split_step one or more times
    void advance_unsplit(double dt, bool check_cfl=true); //Calls unsplit_step one or more times
    
    #ifdef MHD
    virtual void initialize_B_fields(){}
    #endif
    
    //Boundary
    Boundary::BoundaryList boundary = Boundary::BoundaryList();
    virtual ~Grid() = default;
};


class Grid1D: public Grid{
protected:
    ExtendedArray1D<PrimitiveState> w;
    PassiveArray1D q;
public:
    double dx; //Phsyical scale of a grid unit
    
    Grid1D(int size, double dx, int ghosts=1);
    Grid1D(const Grid1D&) = delete; //No copying
    Grid1D& operator=(const Grid1D&) = delete;
    ~Grid1D() = default;
    //Grid access
    //Can take inputs <0 or >= size to access ghost cells
    PrimitiveState& operator[](int k);
    const PrimitiveState& operator[](int k) const;
    int getSize() const, getGhosts() const;
    
    PassiveArray1D& passives(){return q;}
    const PassiveArray1D& passives() const {return q;}
    
    //Advance forward in time
    void split_step(double dt) override;
    void unsplit_step(double dt) override;
};

class Grid2D: public Grid{
protected:
    ExtendedArray2D<PrimitiveState> w;
    #ifdef MHD
    ExtendedArray2D<vec3> B;//B fields on the faces
    #endif
    PassiveArray2D q;
public:
    double dx, dy;
    
    Grid2D(int nx, int ny, double dx, double dy, int ghosts=1);
    Grid2D(const Grid2D&) = delete; //No copying
    Grid2D& operator=(const Grid2D&) = delete;
    ~Grid2D() = default;
    
    //Grid access
    //Can take inputs <0 or >= n to access ghost cells
    PrimitiveState& operator[](int i,int j);
    const PrimitiveState& operator[](int i,int j) const;
    int getSizeX() const, getSizeY() const, getGhosts() const;
    #ifdef MHD //Face-normal magnetic fields
    ExtendedArray2D<vec3>& _B(){return B;}
    const ExtendedArray2D<vec3>& _B() const {return B;}
    #endif
    //Passive Scalars
    PassiveArray2D& passives(){return q;}
    const PassiveArray2D& passives() const {return q;}
    
    //Advance Forward in time
    void split_step(double dt) override;
    void unsplit_step(double dt) override;
    #ifdef MHD
    void initialize_B_fields() override;
    #endif
protected:
    int sweep_step = 0;
    
    void advanceX(double dt); //Advance a single split step in X
    void advanceY(double dt); //Advance a single split step in Y
    #ifdef MHD
    void computeBodyAveragedFields(const ExtendedArray2D<vec3>& B);
    #endif
};

class Grid3D: public Grid{
protected:
    ExtendedArray3D<PrimitiveState>  w;
    #ifdef MHD
    ExtendedArray3D<vec3> B;//B fields on the faces
    #endif
    PassiveArray3D q;
public:
    double dx, dy, dz;
    
    Grid3D(int nx, int ny, int nz, double dx, double dy, double dz, int ghosts=2);
    Grid3D(const Grid3D&) = delete; //No copying
    Grid3D& operator=(const Grid3D&) = delete;
    ~Grid3D() = default;
    
    //Grid Access
    //Can take inputs <0 or >= n to access ghost cells
    PrimitiveState& operator[](int i,int j,int k);
    const PrimitiveState& operator[](int i,int j,int k) const;
    int getSizeX() const, getSizeY() const, getSizeZ() const, getGhosts() const;
    #ifdef MHD//Face-normal magnetic fields
    ExtendedArray3D<vec3>& _B(){return B;}
    const ExtendedArray3D<vec3>& _B() const {return B;}
    #endif
    //Passive Scalars
    PassiveArray3D& passives(){return q;}
    const PassiveArray3D& passives() const {return q;}
    
    //Advance Forward in time
    void split_step(double dt) override;
    void unsplit_step(double dt) override;
    #ifdef MHD
    void initialize_B_fields() override;
    #endif
protected:
    int sweep_step = 0;
    
    void advanceX(double dt); //Advance a single split step in X
    void advanceY(double dt); //Advance a single split step in Y
    void advanceZ(double dt); //Advance a single split step in Z
    #ifdef MHD
    void computeBodyAveragedFields(const ExtendedArray3D<vec3>& B);
    #endif
};
}

#endif
