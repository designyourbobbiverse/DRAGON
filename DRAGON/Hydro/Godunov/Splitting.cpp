//
//  Splitting.cpp
//  DRAGON/Hydro/Godunov
//
//  Created by Bobbie Markwick on 16/05/2026.
//  Implementation based mostly on Toro (2009). https://doi.org/10.1007/b79761
//

#include "Grid.hpp"
#include "Boundary.hpp"
#include "CFL.hpp"
#include <math.h>
#include <cassert>


static int validGhosts(int g){
#ifdef MUSCL_Hancock
    return std::max(g, 2);
#else
    return std::max(g, 1);
#endif
}
//MARK: Array Wrappers
Grid2D::Grid2D(int nx_, int ny_, double dx_, double dy_, int g_):  w(nx_, ny_,validGhosts(g_)),
#ifdef MHD
    A(nx_+1, ny_+1,w.getGhosts()),
#endif
    dx(dx_), dy(dy_) { }
PrimitiveState& Grid2D::operator[](int i, int j) { return w[i,j]; }
const PrimitiveState& Grid2D::operator[](int i, int j) const { return w[i,j]; }
int Grid2D::getSizeX() const { return w.getSizeX(); }
int Grid2D::getSizeY() const { return w.getSizeY(); }
int Grid2D::getGhosts() const { return w.getGhosts(); }


Grid3D::Grid3D(int nx_, int ny_, int nz_, double dx_, double dy_, double dz_, int g_): w(nx_, ny_, nz_, validGhosts(g_)) ,
#ifdef MHD
    A(nx_+1, ny_+1, nz_+1, w.getGhosts()),
#endif
    dx(dx_), dy(dy_), dz(dz_) {    }
PrimitiveState& Grid3D::operator[](int i, int j, int k) { return w[i,j,k]; }
const PrimitiveState& Grid3D::operator[](int i, int j, int k) const { return w[i,j,k]; }
int Grid3D::getSizeX() const { return w.getSizeX(); }
int Grid3D::getSizeY() const { return w.getSizeY(); }
int Grid3D::getSizeZ() const { return w.getSizeZ(); }
int Grid3D::getGhosts() const { return w.getGhosts(); }


//MARK: 2D Split
void Grid2D::advance_split(double dt, bool check_cfl){
    while(dt > CONFIG::Timestep_Tolerance){
        //CFL Time Constraint
        double t1 = check_cfl ? std::min(dt,CFL::cfl_time(*this)) : dt;
        dt -= t1;
        
        //Advance (Strang Split), alternating which step comes first
        if (sweep_step++ % 2 == 0) {
            advanceX(t1/2);
            advanceY(t1);
            advanceX(t1/2);
        } else {
            advanceY(t1/2);
            advanceX(t1);
            advanceY(t1/2);
        }
    }
}
//MARK: 3D Split
void Grid3D::advance_split(double dt, bool check_cfl){
    while(dt > CONFIG::Timestep_Tolerance){
        //CFL Time Constraint
        double t1 = check_cfl ? std::min(dt,CFL::cfl_time(*this)) : dt;
        dt -= t1;
        //Advance  (Strang Split), rotating step orders
        switch (sweep_step++ % 6) {
        case 0: //Cyclic XYZ
            advanceX(t1/2);
            advanceY(t1/2);
            advanceZ(t1);
            advanceY(t1/2);
            advanceX(t1/2);
            break;
        case 1: //Cyclic ZXY
            advanceZ(t1/2);
            advanceX(t1/2);
            advanceY(t1);
            advanceX(t1/2);
            advanceZ(t1/2);
            break;
        case 2: //Cyclic YZX
            advanceY(t1/2);
            advanceZ(t1/2);
            advanceX(t1);
            advanceZ(t1/2);
            advanceY(t1/2);
            break;
        case 3: //Anticyclic ZYX
            advanceZ(t1/2);
            advanceY(t1/2);
            advanceX(t1);
            advanceY(t1/2);
            advanceZ(t1/2);
            break;
        case 4: //Anticyclic XZY
            advanceX(t1/2);
            advanceZ(t1/2);
            advanceY(t1);
            advanceZ(t1/2);
            advanceX(t1/2);
            break;
        case 5: //Anticyclic YXZ
            advanceY(t1/2);
            advanceX(t1/2);
            advanceZ(t1);
            advanceX(t1/2);
            advanceY(t1/2);
            break;
        }
    }
}


//MARK: 2D Component Sweeps
void Grid2D::advanceX(double dt){
    boundary.apply(*this); //Apply Boundary Conditions before every sweep
    
    int nx = w.getSizeX(), ny = w.getSizeY(), ghosts = w.getGhosts();
    #ifdef PRESERVE_BUFFERS
    if(buffer_x == nullptr) buffer_x = new Grid1D(nx,dx,ghosts);
    buffer_x->dx = dx; //In case user changes dx for any reason
    buffer_x->boundary = Boundary::Ignore();
    Grid1D& _w = *buffer_x;
    #else
    Grid1D _w(nx, dx, ghosts);
    ExtendedArray1D<PrimitiveState> _B1(nx, ghosts), _B2(nx, ghosts);//Buffers
    #endif
    
    for(int j=-ghosts; j<ny+ghosts; j++){
        for(int i=-ghosts; i<nx+ghosts; i++) _w[i] = w[i,j]; //Copy to a 1D array

        #ifdef PRESERVE_BUFFERS
        _w.advance(dt, false); //Use _w's buffers
        #else
        _w.god_sweep(dt, _B1, _B2); //Sweep through the 1D array
        #endif
        
        for(int i=-ghosts; i<nx+ghosts; i++) w[i,j] = _w[i]; //Copy 1D array back to grid
    }
}
void Grid2D::advanceY(double dt){
    boundary.apply(*this); //Apply Boundary Conditions before every sweep

    int nx = w.getSizeX(), ny = w.getSizeY(), ghosts = w.getGhosts();
    #ifdef PRESERVE_BUFFERS
    if(buffer_y == nullptr) buffer_y = new Grid1D(ny,dy,ghosts);
    buffer_y->dx = dy; //In case user changes dy for any reason
    buffer_y->boundary = Boundary::Ignore();
    Grid1D& _w = *buffer_y;
    #else
    Grid1D _w(ny, dy, ghosts);
    ExtendedArray1D<PrimitiveState> _B1(ny, ghosts), _B2(ny, ghosts);
    #endif
    
    
    for(int i=-ghosts; i<nx+ghosts; i++){
        for(int j=-ghosts; j<ny+ghosts; j++)  _w[j] = w[i,j].swappedXY(); //Dimension swap + copy to a 1D array

        #ifdef PRESERVE_BUFFERS
        _w.advance(dt, false); //Use _w's buffers
        #else
        _w.god_sweep(dt, _B1, _B2); //Sweep through the 1D array
        #endif
        
        for(int j=-ghosts; j<ny+ghosts; j++)  w[i,j] = _w[j].swappedXY(); //Dimension swap back + copy back to grid
    }
}



//MARK: 3D Component Sweeps
void Grid3D::advanceX(double dt){
    boundary.apply(*this); //Apply Boundary Conditions before every sweep

    int nx = w.getSizeX(), ny = w.getSizeY(), nz = w.getSizeZ(), ghosts = w.getGhosts();
    #ifdef PRESERVE_BUFFERS
    if(buffer_x == nullptr) buffer_x = new Grid1D(nx,dx,ghosts);
    buffer_x->dx = dx; //In case user changes dx for any reason
    buffer_x->boundary = Boundary::Ignore();
    Grid1D& _w = *buffer_x;
    #else
    Grid1D _w(nx, dx, ghosts);
    ExtendedArray1D<PrimitiveState> _B1(nx, ghosts), _B2(nx, ghosts);//Buffers
    #endif
    
    for(int k=-ghosts; k<nz+ghosts; k++){
        for(int j=-ghosts; j<ny+ghosts; j++){
            for(int i=-ghosts; i<nx+ghosts; i++) _w[i] = w[i,j,k]; //Copy to a 1D array
            
            #ifdef PRESERVE_BUFFERS
            _w.advance(dt, false); //Use _w's buffers
            #else
            _w.god_sweep(dt, _B1, _B2); //Sweep through the 1D array
            #endif
            
            for(int i=-ghosts; i<nx+ghosts; i++) w[i,j,k] = _w[i]; //Copy back to grid
        }
    }
}
void Grid3D::advanceY(double dt){
    boundary.apply(*this); //Apply Boundary Conditions before every sweep

    int nx = w.getSizeX(), ny = w.getSizeY(), nz = w.getSizeZ(), ghosts = w.getGhosts();
    #ifdef PRESERVE_BUFFERS
    if(buffer_y == nullptr) buffer_y = new Grid1D(ny,dy,ghosts);
    buffer_y->dx = dy; //In case user changes dy for any reason
    buffer_y->boundary = Boundary::Ignore();
    Grid1D& _w = *buffer_y;
    #else
    Grid1D _w(ny, dy, ghosts);
    ExtendedArray1D<PrimitiveState> _B1(ny, ghosts), _B2(ny, ghosts);
    #endif
    
    for(int k=-ghosts; k<nz+ghosts; k++){
        for(int i=-ghosts; i<nx+ghosts; i++) {
            for(int j=-ghosts; j<ny+ghosts; j++) _w[j] = w[i,j,k].swappedXY(); //Dimension swap + copy to a 1D array

            #ifdef PRESERVE_BUFFERS
            _w.advance(dt, false); //Use _w's buffers
            #else
            _w.god_sweep(dt, _B1, _B2); //Sweep through the 1D array
            #endif
            
            for(int j=-ghosts; j<ny+ghosts; j++) w[i,j,k] = _w[j].swappedXY();  //Dimension swap back + copy back to grid
       }
    }
}
void Grid3D::advanceZ(double dt){
    boundary.apply(*this); //Apply Boundary Conditions before every sweep

    int nx = w.getSizeX(), ny = w.getSizeY(), nz = w.getSizeZ(), ghosts = w.getGhosts();
    #ifdef PRESERVE_BUFFERS
    if(buffer_z == nullptr) buffer_z = new Grid1D(nz,dz,ghosts);
    buffer_z->dx = dz; //In case user changes dz for any reason
    buffer_z->boundary = Boundary::Ignore();
    Grid1D& _w = *buffer_z;
    #else
    Grid1D _w(nz,dz,ghosts);
    ExtendedArray1D<PrimitiveState> _B1(nz, ghosts), _B2(nz, ghosts);
    #endif
    
    
    for(int i=-ghosts; i<nx+ghosts; i++) {
        for(int j=-ghosts; j<ny+ghosts; j++) {
            for(int k=-ghosts; k<nz+ghosts; k++) _w[k] = w[i,j,k].swappedXZ(); //Dimension swap + copy to a 1D array
            
            #ifdef PRESERVE_BUFFERS
            _w.advance(dt, false); //Use _w's buffers
            #else
            _w.god_sweep(dt, _B1, _B2); //Sweep through the 1D array
            #endif
            
            for(int k=-ghosts; k<nz+ghosts; k++) w[i,j,k] = _w[k].swappedXZ(); //Dimension swap back + copy back to grid
       }
    }
}
