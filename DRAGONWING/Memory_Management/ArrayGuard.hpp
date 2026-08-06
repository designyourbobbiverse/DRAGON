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

#include "ArrayTypes.hpp"

template<typename X>
concept GridValueType = std::same_as<X, DRAGON::PrimitiveState> || std::same_as<X, DRAGON::ConservativeState> || std::same_as<X, DRAGON::vec3>;

template<typename T> struct is_allowed_grid : std::false_type {};
template<GridValueType X> struct is_allowed_grid<DRAGON::ExtendedArray1D<X>> : std::true_type {};
template<GridValueType X> struct is_allowed_grid<DRAGON::ExtendedArray2D<X>> : std::true_type {};
template<GridValueType X> struct is_allowed_grid<DRAGON::ExtendedArray3D<X>> : std::true_type {};

template<typename T>
concept ValidGridType = is_allowed_grid<std::remove_cvref_t<T>>::value;


namespace DW_Internal{
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

namespace DRAGONWING {
template <ValidGridType T>
class ArrayGuard{
private:
    std::vector<T*> arrays;
public:
    ArrayGuard(std::vector<T*> arrs): arrays(arrs) {}
    ArrayGuard(const ArrayGuard&) = delete;
    ArrayGuard& operator=(const ArrayGuard&) = delete;
    
    T* get(int i){ return arrays[i]; }
    T* operator[](int i){ return arrays[i]; }
    
    void release(){
        while (!arrays.empty()) {
            T* arr = arrays.back();
            arrays.pop_back();
            DW_Internal::releaseArray(arr);
        }
    }
    ~ArrayGuard(){ release(); }
};
}
