//
//  Storage.hpp
//  DRAGON
//
//  Created by Bobbie Markwick on 16/06/2026.
//

#ifndef Storage_hpp
#define Storage_hpp
#include <cassert>

#pragma once

//MARK: 1D Array
template <class T>
struct ExtendedArray1D{
    
    ExtendedArray1D(int size, int ghosts=2): size(size), ghosts(ghosts) {
        x = new T[size+2*ghosts];
    }
    ~ExtendedArray1D(){ }// delete[] x; }
    ExtendedArray1D(const ExtendedArray1D&) = delete; //No cop.ying
    ExtendedArray1D& operator=(const ExtendedArray1D&) = delete;
    
    //Grid access
    T& operator[](int i){
    #ifdef TESTMODE
        assert(i + ghosts >= 0 && i < size+ghosts);
    #endif
        int m = i+ghosts;
        return x[m];
    }
    const T& operator[](int i) const{
    #ifdef TESTMODE
        assert(i + ghosts >= 0 && i < size+ghosts);
    #endif
        int m = i+ghosts;
        return x[m];
    }
    int getSize() const{ return size; }
    int getGhosts() const{ return ghosts; }
private:
    T* x;
    int ghosts, size;
};

//MARK: 2D Array
template <class T>
struct ExtendedArray2D{
    
    ExtendedArray2D(int nx_, int ny_, int ghosts=2): nx(nx_), ny(ny_), ghosts(ghosts) {
        x = new T[(nx+2*ghosts)*(ny+2*ghosts)];
    }
    ~ExtendedArray2D(){ }//delete[] x; }
    ExtendedArray2D(const ExtendedArray2D&) = delete; //No cop.ying
    ExtendedArray2D& operator=(const ExtendedArray2D&) = delete;
    
    //Grid access
    T& operator[](int i,int j){
    #ifdef TESTMODE
        assert(i + ghosts >= 0 && i < nx+ghosts);
        assert(j + ghosts >= 0 && j < ny+ghosts);
    #endif
        int m = (i+ghosts)*(ny+2*ghosts) + (j+ghosts);
        return x[m];
    }
    const T& operator[](int i,int j) const{
    #ifdef TESTMODE
        assert(i + ghosts >= 0 && i < nx+ghosts);
        assert(j + ghosts >= 0 && j < ny+ghosts);
    #endif
        int m = (i+ghosts)*(ny+2*ghosts) + (j+ghosts);
        return x[m];
    }
    int getSizeX() const{ return nx; }
    int getSizeY() const{ return ny; }
    int getGhosts() const{ return ghosts; }
private:
    T* x;
    int ghosts, nx, ny;
};


//MARK: 3D Array
template <class T>
struct ExtendedArray3D{
    
    ExtendedArray3D(int nx_, int ny_, int nz_, int ghosts=2): nx(nx_), ny(ny_), nz(nz_), ghosts(ghosts) {
        x = new T[(nx+2*ghosts)*(ny+2*ghosts)*(nz+2*ghosts)];
    }
    ~ExtendedArray3D(){ }//delete[] x; }
    ExtendedArray3D(const ExtendedArray3D&) = delete; //No cop.ying
    ExtendedArray3D& operator=(const ExtendedArray3D&) = delete;
    
    //Grid access
    T& operator[](int i,int j, int k){
    #ifdef TESTMODE
        assert(i + ghosts >= 0 && i < nx+ghosts);
        assert(j + ghosts >= 0 && j < ny+ghosts);
        assert(k + ghosts >= 0 && k < nz+ghosts);
    #endif
        int m = ((i+ghosts)*(ny+2*ghosts) + (j+ghosts)) * (nz+2*ghosts) + (k+ghosts);
        return x[m];
    }
    const T& operator[](int i,int j,int k) const{
    #ifdef TESTMODE
        assert(i + ghosts >= 0 && i < nx+ghosts);
        assert(j + ghosts >= 0 && j < ny+ghosts);
        assert(k + ghosts >= 0 && k < nz+ghosts);
    #endif
        int m = ((i+ghosts)*(ny+2*ghosts) + (j+ghosts)) * (nz+2*ghosts) + (k+ghosts);
        return x[m];
    }
    int getSizeX() const{ return nx; }
    int getSizeY() const{ return ny; }
    int getSizeZ() const{ return nz; }
    int getGhosts() const{ return ghosts; }
private:
    T* x;
    int ghosts, nx, ny, nz;
};

#endif /* Storage_hpp */
