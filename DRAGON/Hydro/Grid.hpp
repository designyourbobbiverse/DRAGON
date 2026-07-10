//
//  Grid.hpp
//  DRAGON/Hydro
//  User-Facing Header file
//
//  Created by Bobbie Markwick on 12/05/2026.
//

#ifndef Godunov_hpp
#define Godunov_hpp

#include "FluidElement.hpp"
#include "ExtendedArray.hpp"
#include "Boundary.hpp"
#include "Config.h"

class Grid{
public:
    virtual void advance(double dt, bool check_cfl=true) = 0;
    virtual ~Grid() = default;
    //Boundary
    Boundary::BoundaryList boundary = Boundary::BoundaryList();
};


class Grid1D: public Grid{
protected:
    ExtendedArray1D<PrimitiveState> w;
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
        
    //Advance forward in time
    void advance(double dt, bool check_cfl=true);
};

class Grid2D: public Grid{
protected:
    ExtendedArray2D<PrimitiveState> w;
#ifdef MHD
    ExtendedArray2D<vec3> A;//Magnetic Vector Potential
#endif
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
    #ifdef MHD
    //Access Edge Magnetic potentials. Only Az is used in 2D
    //A[i,j] is on the corner w[i-1/2,j-1/2]
    ExtendedArray2D<vec3>& _A(){return A;}
    #endif
    
    //Advance Forward in time
    void advance(double dt, bool check_cfl = true); //Automatically split or unsplit based whether DIMENSION_UNSPLIT is set in Config.h
    void advance_split(double dt, bool check_cfl = true);
    void advance_unsplit(double dt, bool check_cfl = true);

    #ifdef MHD
    void initialize_B_fields();
    #endif
protected:
    int sweep_step = 0;

    void advanceX(double dt); //Advance a single split step in X
    void advanceY(double dt); //Advance a single split step in Y
    void advanceXY(double dt); //Advance a single unsplit step
    #ifdef MHD
    void computeBodyAveragedFields(const ExtendedArray2D<vec3>& B);
    #endif
};

class Grid3D: public Grid{
protected:
    ExtendedArray3D<PrimitiveState>  w;
#ifdef MHD
    ExtendedArray3D<vec3> A;//Magnetic Vector Potential
#endif
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
    #ifdef MHD
    //Access Edge Magnetic potentials.
    //A[i,j,k] is the corner w[i-1/2,j-1/2,k-1/2] to each of the 3 adjacent corners of w[i,j,k]
    ExtendedArray3D<vec3>& _A(){return A;}
    #endif
    
    //Advance Forward in time
    void advance(double dt, bool check_cfl = true);//Automatically split or unsplit based whether DIMENSION_UNSPLIT is set in Config.h
    void advance_split(double dt, bool check_cfl = true);
    void advance_unsplit(double dt, bool check_cfl = true);
    
    #ifdef MHD
    void initialize_B_fields();
    #endif
protected:
    int sweep_step = 0;
    
    void advanceX(double dt); //Advance a single split step in X
    void advanceY(double dt); //Advance a single split step in Y
    void advanceZ(double dt); //Advance a single split step in Z
    void advanceXYZ(double dt); //Advance a single unsplit step
    #ifdef MHD
    void computeBodyAveragedFields(const ExtendedArray3D<vec3>& B);
    #endif
};


#endif
