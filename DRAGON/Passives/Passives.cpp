//
//  Passives.cpp
//  DRAGON/Passives
//
//  Created by Bobbie Markwick on 31/07/2026.
//

#include "Passives.hpp"
using namespace DRAGON;

const double advect_tol = 1e-18;

//MARK: Access
double& PassiveArray::lookup(PassiveSet& q, const std::string& key) const { return q[keys.at(key)]; }
const double& PassiveArray::lookup(const PassiveSet& q, const std::string& key) const { return q[keys.at(key)]; }
std::size_t PassiveArray::count() const { return keys.size(); }

//Get all the scalars in a given cell
PassiveSet& PassiveArray1D::operator[](int i){ return q[i]; }
const PassiveSet& PassiveArray1D::operator[](int i) const { return q[i]; }
PassiveSet& PassiveArray2D::operator[](int i, int j){ return q[i,j]; }
const PassiveSet& PassiveArray2D::operator[](int i, int j) const { return q[i,j]; }
PassiveSet& PassiveArray3D::operator[](int i, int j, int k){ return q[i,j,k]; }
const PassiveSet& PassiveArray3D::operator[](int i, int j, int k) const { return q[i,j,k]; }
//Get the nth scalar in a given cell
double& PassiveArray1D::operator[](int i, int n){ return q[i][n]; }
const double& PassiveArray1D::operator[](int i, int n) const { return q[i][n]; }
double& PassiveArray2D::operator[](int i, int j, int n){ return q[i,j][n]; }
const double& PassiveArray2D::operator[](int i, int j, int n) const { return q[i,j][n]; }
double& PassiveArray3D::operator[](int i, int j, int k, int n){ return q[i,j,k][n]; }
const double& PassiveArray3D::operator[](int i, int j, int k, int n) const { return q[i,j,k][n]; }
//Get the scalar assocaited with "key" in a given cell
double& PassiveArray1D::operator[](int i, const std::string& key){ return lookup(q[i], key); }
const double& PassiveArray1D::operator[](int i, const std::string& key) const { return lookup(q[i], key); }
double& PassiveArray2D::operator[](int i, int j, const std::string& key){ return lookup(q[i,j], key); }
const double& PassiveArray2D::operator[](int i, int j, const std::string& key) const { return lookup(q[i,j], key); }
double& PassiveArray3D::operator[](int i, int j, int k, const std::string& key){ return lookup(q[i,j,k], key); }
const double& PassiveArray3D::operator[](int i, int j, int k, const std::string& key) const { return lookup(q[i,j,k], key); }


//MARK: Cloning
void PassiveArray1D::clone(const PassiveArray1D& arr){
    keys = arr.keys;
    q.clone(arr.q);
}
void PassiveArray2D::clone(const PassiveArray2D& arr){
    keys = arr.keys;
    q.clone(arr.q);
}
void PassiveArray3D::clone(const PassiveArray3D& arr){
    keys = arr.keys;
    q.clone(arr.q);
}



//MARK: Adding Scalars
void PassiveArray1D::add(const std::string& key){
    if(keys.contains(key)) return; //scalar already exists
    keys[key] = count();
    const int nx = q.getSize(), g = q.getGhosts();
    for (int i=-g; i<nx+g; i++) {
        q[i].push_back(0);
    }
}
void PassiveArray2D::add(const std::string& key){
    if(keys.contains(key)) return; //scalar already exists
    keys[key] = count();
    const int nx = q.getSizeX(), ny = q.getSizeY(), g = q.getGhosts();
    for (int i=-g; i<nx+g; i++) {
        for (int j=-g; j<ny+g; j++) {
            q[i,j].push_back(0);
        }
    }
}
void PassiveArray3D::add(const std::string& key){
    if(keys.contains(key)) return; //scalar already exists
    keys[key] = count();
    const int nx = q.getSizeX(), ny = q.getSizeY(), nz = q.getSizeZ(), g = q.getGhosts();
    for (int i=-g; i<nx+g; i++) {
        for (int j=-g; j<ny+g; j++) {
            for (int k=-g; k<nz+g; k++) {
                q[i,j,k].push_back(0);
            }
        }
    }
}


//MARK: Removing Scalars
std::size_t PassiveArray::remove(const std::string& key){
    if(!keys.contains(key)) return  -1; //scalar doesn't exist
    auto idx = keys[key];
    keys.erase(key);
    for (auto& [key, val]: keys) { //Decrement anything greater than idx
        if(val > idx) val--;
    }
    return idx;
}
std::size_t PassiveArray1D::remove(const std::string& key){
    if(!keys.contains(key)) return -1; //scalar doesn't exist
    auto idx = PassiveArray::remove(key);
    
    const int nx = q.getSize(), g = q.getGhosts();
    for (int i=-g; i<nx+g; i++) {
        q[i].erase(q[i].begin() + idx);
    }
    
    return idx;
}
std::size_t PassiveArray2D::remove(const std::string& key){
    if(!keys.contains(key)) return  -1; //scalar doesn't exist
    auto idx = PassiveArray::remove(key);

    const int nx = q.getSizeX(), ny = q.getSizeY(), g = q.getGhosts();
    for (int i=-g; i<nx+g; i++) {
        for (int j=-g; j<ny+g; j++) {
            q[i,j].erase(q[i,j].begin()+idx);
        }
    }
    
    return idx;
}
std::size_t PassiveArray3D::remove(const std::string& key){
    if(!keys.contains(key)) return  -1; //scalar doesn't exist
    auto idx = PassiveArray::remove(key);
    
    const int nx = q.getSizeX(), ny = q.getSizeY(), nz = q.getSizeZ(), g = q.getGhosts();
    for (int i=-g; i<nx+g; i++) {
        for (int j=-g; j<ny+g; j++) {
            for (int k=-g; k<nz+g; k++) {
                q[i,j,k].erase(q[i,j,k].begin()+idx);
            }
        }
    }
    
    return idx;
}

//MARK: Advection

std::unique_ptr<PassiveArray1D> PassiveArray1D::advected(const ExtendedArray1D<ConservativeState>& F, const FluidArray1D& w_old, const FluidArray1D& w_new, double dt_dx){
    const int nx = q.getSize(), g = q.getGhosts();
    //Make a copy
    auto advected = std::make_unique<PassiveArray1D>(nx, g);
    advected->clone(*this);
    //Convert to conservative
    for(int i=-g; i<nx+g; i++){
        for(unsigned int n = 0; n < count(); n++){
            (*advected)[i,n] *= w_old[i].rho;
        }
    }
    //Apply fluxes
    for(int i=-g+1; i<nx+g; i++){
        for(unsigned int n = 0; n < count(); n++){
            auto f = F[i].rho * dt_dx;
            if(F[i].rho > advect_tol) { //Left to Right
                (*advected)[i,n] += f * q[i-1][n];
                (*advected)[i-1,n] -= f * q[i-1][n];
            } else if(F[i].rho < -advect_tol){ //Right to Left
                (*advected)[i,n] += f * q[i][n];
                (*advected)[i-1,n] -= f * q[i][n];
            }
        }
    }
    //Convert to primitive
    for(int i=-g; i<nx+g; i++){
        for(unsigned int n = 0; n < count(); n++){
            (*advected)[i,n] /= w_new[i].rho;
        }
    }
    return advected;
}
std::unique_ptr<PassiveArray2D> PassiveArray2D::advected(const FluxArray2D& F_X, const FluxArray2D& F_Y, const FluidArray2D& w_old, const FluidArray2D& w_new, double dt_dx, double dt_dy){
    const int nx = q.getSizeX(), ny = q.getSizeY(), g = q.getGhosts();
    //Make a copy
    auto advected = std::make_unique<PassiveArray2D>(nx, ny, g);
    advected->clone(*this);
    //Convert to conservative
    for(int i=-g; i<nx+g; i++){
        for (int j=-g; j<ny+g; j++) {
            for(unsigned int n = 0; n < count(); n++){
                (*advected)[i,j,n] *= w_old[i,j].rho;
            }
        }
    }
    //Apply fluxes
    for(int i=-g+1; i<nx+g; i++){
        for (int j=-g+1; j<ny+g; j++) {
            for(unsigned int n = 0; n < count(); n++){
                //X fluxes
                double fx = F_X[i,j].rho * dt_dx;
                if(F_X[i,j].rho > advect_tol) { //Left to Right
                    (*advected)[i,j,n] += fx * q[i-1,j][n];
                    (*advected)[i-1,j,n] -= fx * q[i-1,j][n];
                } else if(F_X[i,j].rho < -advect_tol){ //Right to Left
                    (*advected)[i,j,n] += fx * q[i,j][n];
                    (*advected)[i-1,j,n] -= fx * q[i,j][n];
                }
                //Y fluxes
                double fy = F_Y[i,j].rho * dt_dy;
                if(F_Y[i,j].rho > advect_tol) { //Left to Right
                    (*advected)[i,j,n] += fy * q[i,j-1][n];
                    (*advected)[i,j-1,n] -= fy * q[i,j-1][n];
                } else if(F_Y[i,j].rho < -advect_tol){ //Right to Left
                    (*advected)[i,j,n] += fy * q[i,j][n];
                    (*advected)[i,j-1,n] -= fy * q[i,j][n];
                }
            }
        }
    }
    //Convert to primitive
    for(int i=-g; i<nx+g; i++){
        for (int j=-g; j<ny+g; j++) {
            for(unsigned int n = 0; n < count(); n++){
                (*advected)[i,j,n] /= w_new[i,j].rho;
            }
        }
    }
    return advected;
}
std::unique_ptr<PassiveArray3D> PassiveArray3D::advected(const FluxArray3D& F_X, const FluxArray3D& F_Y, const FluxArray3D& F_Z, const FluidArray3D& w_old, const FluidArray3D& w_new, double dt_dx, double dt_dy, double dt_dz){
    const int nx = q.getSizeX(), ny = q.getSizeY(), nz = q.getSizeZ(), g = q.getGhosts();
    //Make a copy
    auto advected = std::make_unique<PassiveArray3D>(nx, ny, nz, g);
    advected->clone(*this);
    //Convert to conservative
    for(int i=-g; i<nx+g; i++){
        for (int j=-g; j<ny+g; j++) {
            for (int k=-g; k<nz+g; k++) {
                for(unsigned int n = 0; n < count(); n++){
                    (*advected)[i,j,k,n] *= w_old[i,j,k].rho;
                }
            }
        }
    }
    //Apply fluxes
    for(int i=-g+1; i<nx+g; i++){
        for (int j=-g+1; j<ny+g; j++) {
            for (int k=-g+1; k<nz+g; k++) {
                for(unsigned int n = 0; n < count(); n++){
                    //X fluxes
                    double fx = F_X[i,j,k].rho * dt_dx;
                    if(F_X[i,j,k].rho > advect_tol) { //Left to Right
                        (*advected)[i,j,k,n] += fx * q[i-1,j,k][n];
                        (*advected)[i-1,j,k,n] -= fx * q[i-1,j,k][n];
                    } else if(F_X[i,j,k].rho < -advect_tol){ //Right to Left
                        (*advected)[i,j,k,n] += fx * q[i,j,k][n];
                        (*advected)[i-1,j,k,n] -= fx * q[i,j,k][n];
                    }
                    //Y fluxes
                    double fy = F_Y[i,j,k].rho * dt_dy;
                    if(F_Y[i,j,k].rho > advect_tol) { //Left to Right
                        (*advected)[i,j,k,n] += fy * q[i,j-1,k][n];
                        (*advected)[i,j-1,k,n] -= fy * q[i,j-1,k][n];
                    } else if(F_Y[i,j,k].rho < -advect_tol){ //Right to Left
                        (*advected)[i,j,k,n] += fy * q[i,j,k][n];
                        (*advected)[i,j-1,k,n] -= fy * q[i,j,k][n];
                    }
                    //Z fluxes
                    double fz = F_Z[i,j,k].rho * dt_dz;
                    if(F_Z[i,j,k].rho > advect_tol) { //Left to Right
                        (*advected)[i,j,k,n] += fz * q[i,j,k-1][n];
                        (*advected)[i,j,k-1,n] -= fz * q[i,j,k-1][n];
                    } else if(F_Z[i,j,k].rho < -advect_tol){ //Right to Left
                        (*advected)[i,j,k,n] += fz * q[i,j,k][n];
                        (*advected)[i,j,k-1,n] -= fz * q[i,j,k][n];
                    }
                }
            }
        }
    }
    //Convert to primitive
    for(int i=-g; i<nx+g; i++){
        for (int j=-g; j<ny+g; j++) {
            for (int k=-g; k<nz+g; k++) {
                for(unsigned int n = 0; n < count(); n++){
                    (*advected)[i,j,k,n] /= w_new[i,j,k].rho;
                }
            }
        }
    }
    return advected;
}
