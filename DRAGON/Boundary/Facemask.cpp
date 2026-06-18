//
//  Boundary/Facemask.cpp
//  DRAGON
//
//  Created by Bobbie Markwick on 29/05/2026.
//

#include "Boundary.hpp"

int ::GhostFill::get_faces() const { return faces; }

//Convert a String such as "X-" or "YZ" into a bitmask
int Boundary::face_mask(std::string s){
    int mask = 0;
    s += " ";
    for(auto i = s.begin(); i+1 < s.end(); i++){
        char current = *i, next = *(i + 1);
        if(current == 'X' || current == 'x'){
            switch(next){
            case '+': mask |= X_positive; break;
            case '-': mask |= X_negative; break;
            default: mask |= X;
            }
        } else if(current == 'Y' || current == 'y'){
            switch(next){
            case '+': mask |= Y_positive; break;
            case '-': mask |= Y_negative; break;
            default: mask |= Y;
            }
        } else if(current == 'Z' || current == 'z'){
            switch(next){
            case '+': mask |= Z_positive; break;
            case '-': mask |= Z_negative; break;
            default: mask |= Z;
            }
        } else if(current == 'L' || current == 'l'){
            mask |= X_negative;
        } else if(current == 'R' || current == 'r'){
            mask |= X_positive;
        } else if(current == 'U' || current == 'u'){
            mask |= Y_positive;
        } else if(current == 'D' || current == 'd'){
            mask |= Y_negative;
        }
    }
    return mask;
}
