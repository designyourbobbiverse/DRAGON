//
//  Constants.h
//  DRAGON
//
//  Created by Bobbie Markwick on 28/04/2026.
//

#ifndef Constants_h
#define Constants_h

#include <cmath>

//MARK: Physical Constants

constexpr double _gamma = 5.0/3.0;



//MARK: Gamma-dependent constants
constexpr double  _Gm1_2G = (_gamma - 1.0)/(2.0 * _gamma);
constexpr double _Gp1_2G = (_gamma + 1.0)/(2.0*_gamma);
constexpr double _2G_Gm1 = 2.0*_gamma/(_gamma - 1.0);
constexpr double _G_Gm1 = _gamma / (_gamma - 1.0);
constexpr double _Gm1_2 = (_gamma - 1.0)/2.0;
constexpr double _2_Gm1 = 2.0/(_gamma - 1.0);
constexpr double _2_Gp1 = 2.0/(_gamma + 1.0);
constexpr double _Gm1_Gp1 = (_gamma - 1.0)/(_gamma + 1.0);
constexpr double _Ginv = 1.0/_gamma;


#endif /* Constants_h */
