//
//  Constants.h
//  DRAGON/Examples/MHD/MHDRotor
//
//  Created by Bobbie Markwick on 28/04/2026.
//

#ifndef Constants_h
#define Constants_h
#include <cmath>





//MARK: Physical Constants
namespace DRAGON{

constexpr double gamma = 1.4;


//******************************************************************//
//*DO NOT TOUCH anything below this line (but feel free to read it)*//
//******************************************************************//

//MARK: Pi-related constants
constexpr double _1_4pi = 0.25/M_PI;
constexpr double _1_8pi = 0.125/M_PI;
const double sq4pi = std::sqrt(4.0*M_PI);
const double sq8pi = std::sqrt(8.0*M_PI);

//MARK: Gamma-dependent constants
constexpr double  _Gm1_2G = (gamma - 1.0)/(2.0 * gamma);
constexpr double _Gp1_2G = (gamma + 1.0)/(2.0*gamma);
constexpr double _2G_Gm1 = 2.0*gamma/(gamma - 1.0);
constexpr double _G_Gm1 = gamma / (gamma - 1.0);
constexpr double _Gm1_2 = (gamma - 1.0)/2.0;
constexpr double _2_Gm1 = 2.0/(gamma - 1.0);
constexpr double _2_Gp1 = 2.0/(gamma + 1.0);
constexpr double _Gm1_Gp1 = (gamma - 1.0)/(gamma + 1.0);
constexpr double _Ginv = 1.0/gamma;
}


#endif /* Constants_h */
