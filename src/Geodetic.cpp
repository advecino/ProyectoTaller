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
        throw std::invalid_argument("invalid input in Geodetic constructor");
    }

    // Iteration
    double dZ = e2 * Z;
    double ZdZ, Nh, SinPhi, N, dZ_new;

    while(true) {
        ZdZ = Z + dZ;
        Nh = sqrt(rho2 + ZdZ*ZdZ);
        SinPhi = ZdZ / Nh;                //Sine of geodetic latitude
        N = R_Earth / sqrt(1.0 - e2*SinPhi*SinPhi);
        dZ_new = N * e2 * SinPhi;

        if (std::abs(dZ - dZ_new) < epsRequ) {
            break;
        }
        dZ = dZ_new;
    }

    return GeodeticCoords(
            atan2(Y, X),       // Longitud
            atan2(ZdZ, sqrt(rho2)), // Latitud
            Nh - N             // Altitud
    );
}
