//
//  Passives.hpp
//  DRAGON/Passives
//
//  Created by Bobbie Markwick on 31/07/2026.
//

#ifndef Passives_hpp
#define Passives_hpp


#include <vector>
#include <map>
#include <string>
#include "Hydro/ExtendedArray/ArrayTypes.hpp"


namespace DRAGON{
typedef std::vector<double> PassiveSet;

struct PassiveArray{
protected:
    std::map<std::string, std::size_t> keys;
public:
    double& lookup(PassiveSet& q, const std::string& key) const;
    const double& lookup(const PassiveSet& q, const std::string& key) const;

    virtual void add(const std::string& key) = 0;
    virtual std::size_t remove(const std::string& key); //Returns index of removed key
    std::size_t count() const;
    
    virtual ~PassiveArray() = default;
};

struct PassiveArray1D: public PassiveArray{
private:
    ExtendedArray1D<PassiveSet> q;
    friend class DistGrid1D; //For syncing keys between parent/child
public:
    PassiveArray1D(int nx, int g=1): q(nx, g) {}
    void clone(const PassiveArray1D& arr);

    void add(const std::string& key) override;
    std::size_t remove(const std::string& key) override;
    
    PassiveSet& operator[](int i);
    const PassiveSet& operator[](int i) const;
    double& operator[](int i, int n);
    const double& operator[](int i, int n) const;
    double& operator[](int i, const std::string& key);
    const double& operator[](int i, const std::string& key) const;
    
    
    std::unique_ptr<PassiveArray1D> advected(const ExtendedArray1D<ConservativeState>& F, const FluidArray1D& w_old, const FluidArray1D& w_new, double dt_dx);
};

struct PassiveArray2D: public PassiveArray{
private:
    ExtendedArray2D<PassiveSet> q;
    friend class DistGrid2D; //For syncing keys between parent/child
public:
    PassiveArray2D(int nx, int ny, int g=1): q(nx, ny, g) {}
    void clone(const PassiveArray2D& arr);

    void add(const std::string& key) override;
    std::size_t remove(const std::string& key) override;

    
    PassiveSet& operator[](int i, int j);
    const PassiveSet& operator[](int i, int j) const;
    double& operator[](int i, int j, int n);
    const double& operator[](int i, int j, int n) const;
    double& operator[](int i, int j, const std::string& key);
    const double& operator[](int i, int j, const std::string& key) const;
    
    std::unique_ptr<PassiveArray2D> advected(const FluxArray2D& F_X, const FluxArray2D& F_Y, const FluidArray2D& w_old, const FluidArray2D& w_new, double dt_dx, double dt_dy);

};

struct PassiveArray3D: public PassiveArray{
private:
    ExtendedArray3D<PassiveSet> q;
    friend class DistGrid3D; //For syncing keys between parent/child
public:
    PassiveArray3D(int nx, int ny, int nz, int g=1): q(nx, ny, nz, g) {}
    void clone(const PassiveArray3D& arr);

    void add(const std::string& key) override;
    std::size_t remove(const std::string& key) override;

    PassiveSet& operator[](int i, int j, int k);
    const PassiveSet& operator[](int i, int j, int k) const;
    double& operator[](int i, int j, int k, int n);
    const double& operator[](int i, int j, int k, int n) const;
    double& operator[](int i, int j, int k, const std::string& key);
    const double& operator[](int i, int j, int k, const std::string& key) const;
    
    std::unique_ptr<PassiveArray3D> advected(const FluxArray3D& F_X, const FluxArray3D& F_Y, const FluxArray3D& F_Z, const FluidArray3D& w_old, const FluidArray3D& w_new, double dt_dx, double dt_dy, double dt_dz);

};

}



#endif
