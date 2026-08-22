//
//  ArrayGuard.hpp
//  DRAGONWING/Memory_Management
//
//  Created by Bobbie Markwick on 10/07/2026.
//

#pragma once
#include <vector>

#include <concepts>     //for ValidGridType template
#include <type_traits>  //for ValidGridType template

#include "Hydro/ExtendedArray/ArrayTypes.hpp"

namespace DRAGONWING {
//ArrayGuard supports ExtendedArrays (1D, 2D, or 3D) containing any of the types defined in FluidElement/FluidElement.hpp
template<typename X> concept GridValueType =
    std::same_as<X, DRAGON::PrimitiveState> || std::same_as<X, DRAGON::ConservativeState> || std::same_as<X, DRAGON::vec3>;

template<typename T> struct is_allowed_grid : std::false_type {};
template<GridValueType X> struct is_allowed_grid<DRAGON::ExtendedArray1D<X>> : std::true_type {};
template<GridValueType X> struct is_allowed_grid<DRAGON::ExtendedArray2D<X>> : std::true_type {};
template<GridValueType X> struct is_allowed_grid<DRAGON::ExtendedArray3D<X>> : std::true_type {};
template<typename T> concept ValidGridType = is_allowed_grid<std::remove_cvref_t<T>>::value;

//Release an array back to the shared memory pool. These generally should not be called outside of DRAGONWING
namespace Internal{
void releaseArray(DRAGON::ExtendedArray1D<DRAGON::PrimitiveState>* arr);
void releaseArray(DRAGON::ExtendedArray2D<DRAGON::PrimitiveState>* arr);
void releaseArray(DRAGON::ExtendedArray3D<DRAGON::PrimitiveState>* arr);
void releaseArray(DRAGON::ExtendedArray1D<DRAGON::ConservativeState>* arr);
void releaseArray(DRAGON::ExtendedArray2D<DRAGON::ConservativeState>* arr);
void releaseArray(DRAGON::ExtendedArray3D<DRAGON::ConservativeState>* arr);
void releaseArray(DRAGON::ExtendedArray1D<DRAGON::vec3>* arr);
void releaseArray(DRAGON::ExtendedArray2D<DRAGON::vec3>* arr);
void releaseArray(DRAGON::ExtendedArray3D<DRAGON::vec3>* arr);
}

//ArrayGuard provides a batch of ExtendedArrays.
//When ArrayGuard goes out of scope (or earlier if release() is called), these arrays are returned to the shared memory pool
template <ValidGridType T>
class ArrayGuard{
private:
    std::vector<T*> arrays;
public:
    ArrayGuard(std::vector<T*> arrs): arrays(arrs) {}
    ArrayGuard(const ArrayGuard&) = delete; //No copying
    ArrayGuard& operator=(const ArrayGuard&) = delete;
    
    //Access the arrays in the batch. Calling this (or using the arrays) after release() is undefined
    T* get(int i){ return arrays[i]; }
    T* operator[](int i){ return arrays[i]; }
    
    //Return all of the arrays in the batch to the shared memory pool
    void release(){
        while (!arrays.empty()) {
            T* arr = arrays.back();
            arrays.pop_back();
            DRAGONWING::Internal::releaseArray(arr);
        }
    }
    ~ArrayGuard(){ release(); }
};
}
