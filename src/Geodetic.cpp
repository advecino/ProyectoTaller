//
// Created by adria on 11/05/2025.
//

#include <limits>
#include <stdexcept>
#include "../include/Geodetic.h"
#include <cmath>

/*
%--------------------------------------------------------------------------
%
% Geodetic.m
%
% Purpose:
%   geodetic coordinates (Longitude [rad], latitude [rad], altitude [m])
%   from given position vector (r [m])
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------*/
GeodeticCoords Geodetic(const Matrix& r) {
    const double eps = std::numeric_limits<double>::epsilon();
    const double epsRequ = eps * R_Earth; //Convergence criterion
    const double e2 = f_Earth * (2.0 - f_Earth); //Square of eccentricity


    double X = r(1,1);  //Cartesian coordinates
    double Y = r(2,1);
    double Z = r(3,1);
    double rho2 = X*X + Y*Y;   //quare of distance from z-axis

    if (r.norm() == 0.0) {
        return GeodeticCoords(0.0, 0.0, -R_Earth);
    }

    // Iteration
    double dZ = e2 * Z;
    double ZdZ, Nh, SinPhi, N;

    while(true) {
        ZdZ = Z + dZ;
        Nh = sqrt(rho2 + ZdZ*ZdZ);
        SinPhi = ZdZ / Nh;                //Sine of geodetic latitude
        N = R_Earth / sqrt(1.0 - e2*SinPhi*SinPhi);
        double dZ_new = N * e2 * SinPhi;

        if (std::abs(dZ - dZ_new) < epsRequ) {
            dZ = dZ_new;
            break;
        }
        dZ = dZ_new;
    }

    ZdZ = Z + dZ;
    Nh  = std::sqrt(rho2 + ZdZ*ZdZ);
    SinPhi = ZdZ / Nh;
    N  = R_Earth / std::sqrt(1.0 - e2 * SinPhi * SinPhi);

    double lon = std::atan2(Y, X);
    double lat = std::atan2(ZdZ, std::sqrt(rho2));
    double h   = Nh - N;

    return GeodeticCoords(lon, lat, h);
}
