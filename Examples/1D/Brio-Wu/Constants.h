//
//  Constants.h
//  DRAGON/Examples/1D/Brio-Wu
//
//  Created by Bobbie Markwick on 28/04/2026.
//

#ifndef Constants_h
#define Constants_h

#include <cmath>

//MARK: Physical Constants

constexpr double _gamma = 2.0;



//******************************************************************//
//*DO NOT TOUCH anything below this line (but feel free to read it)*//
//******************************************************************//

//MARK: Pi-related constants
constexpr double _2pi = 0.5/M_PI;
constexpr double _4pi = 0.25/M_PI;
constexpr double _8pi = 0.125/M_PI;
const double sq2pi = sqrt(2.0*M_PI);
const double sq4pi = sqrt(4.0*M_PI);
const double sq8pi = sqrt(8.0*M_PI);
const double inv_sq2pi = 1.0/sq2pi;
const double inv_sq4pi = 1.0/sq4pi;
const double inv_sq8pi = 1.0/sq8pi;

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
