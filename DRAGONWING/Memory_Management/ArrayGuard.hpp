//
//  ForwardDecls.h
//  DRAGON
//
//  Created by Bobbie Markwick on 10/07/2026.
//

#pragma once
#include <concepts>
#include <type_traits>


struct PrimitiveState;
struct ConservativeState;
struct vec3;
template <class T> struct ExtendedArray1D;
template <class T> struct ExtendedArray2D;
template <class T> struct ExtendedArray3D;

template<typename X>
concept GridValueType = std::same_as<X, PrimitiveState> || std::same_as<X, ConservativeState> || std::same_as<X, vec3>;

template<typename T> struct is_allowed_grid : std::false_type {};
template<GridValueType X> struct is_allowed_grid<ExtendedArray1D<X>> : std::true_type {};
template<GridValueType X> struct is_allowed_grid<ExtendedArray2D<X>> : std::true_type {};
template<GridValueType X> struct is_allowed_grid<ExtendedArray3D<X>> : std::true_type {};

template<typename T>
concept ValidGridType = is_allowed_grid<std::remove_cvref_t<T>>::value;


namespace DW_Internal{
void releaseArray(ExtendedArray1D<PrimitiveState>* arr);
void releaseArray(ExtendedArray2D<PrimitiveState>* arr);
void releaseArray(ExtendedArray3D<PrimitiveState>* arr);
void releaseArray(ExtendedArray1D<ConservativeState>* arr);
void releaseArray(ExtendedArray2D<ConservativeState>* arr);
void releaseArray(ExtendedArray3D<ConservativeState>* arr);
void releaseArray(ExtendedArray1D<vec3>* arr);
void releaseArray(ExtendedArray2D<vec3>* arr);
void releaseArray(ExtendedArray3D<vec3>* arr);
}

namespace DRAGONWING {
template <ValidGridType T>
class ArrayGuard{
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
private:
    std::vector<T*> arrays;
};
}
