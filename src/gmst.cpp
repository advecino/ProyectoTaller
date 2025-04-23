/*%--------------------------------------------------------------------------
%
% Purpose:
%   Greenwich Mean Sidereal Time
%
% Input:
%   Mjd_UT1    Modified Julian Date UT1
%
% Output:
%   gmstime	   GMST in [rad]
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------*/

/**
 * @file GMST.cpp
 * @brief Implementación del cálculo del Tiempo Sidéreo Medio de Greenwich
 */

#include "../include/gmst.h"

#include <cmath>

const double Secs = 86400.0;         // Seconds per day
const double MJD_J2000 = 51544.5;    // Modified Julian Date of J2000

double gmst(double Mjd_UT1) {
    double Mjd_0 = floor(Mjd_UT1);
    double UT1 = Secs * (Mjd_UT1 - Mjd_0);    // [s]
    double T_0 = (Mjd_0 - MJD_J2000) / 36525.0;
    double T = (Mjd_UT1 - MJD_J2000) / 36525.0;

    double gmst_sec = 24110.54841 + 8640184.812866 * T_0
                      + 1.002737909350795 * UT1
                        + (0.093104-6.2e-6*T)*T*T;  // [s]

    return 2 * M_PI * Frac(gmst_sec / Secs);  // [rad], 0..2pi
}