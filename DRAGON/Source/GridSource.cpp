//
//  GridSource.cpp
//  DRAGON/Sources
//
//  Created by Bobbie Markwick on 08/09/2026.
//

#include "Grid.hpp"

using namespace DRAGON;


 void Grid1D::advanceSource(double dt, bool ghosts){
     if(sources.count() == 0) return;
     
     const int nx = getSize(), g = ghosts ? getGhosts() : 0;
     for(int i=-g; i<nx+g; i++){
         w[i] += sources.integrate(dt, w[i]);
     }
 }
 void Grid2D::advanceSource(double dt, bool ghosts){
     if(sources.count() == 0) return;
     
     const int nx = getSizeX(), ny = getSizeY(), g = ghosts ? getGhosts() : 0;
     for(int i=-g; i<nx+g; i++){
         for(int j=-g; j<ny+g; j++){
             w[i,j] +=  sources.integrate(dt, w[i,j]);
         }
     }
 }
 void Grid3D::advanceSource(double dt, bool ghosts){
     if(sources.count() == 0) return;
     
     const int nx = getSizeX(), ny = getSizeY(), nz = getSizeZ(), g = ghosts ? getGhosts() : 0;
     for(int i=-g; i<nx+g; i++){
         for(int j=-g; j<ny+g; j++){
             for(int k=-g; k<nz+g; k++){
                 w[i,j,k] +=  sources.integrate(dt, w[i,j,k]);
             }
         }
     }
 }
