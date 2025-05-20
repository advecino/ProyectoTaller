#include "../include/AzElPa.h"
#include <cmath>
#include <stdexcept>

/*
%--------------------------------------------------------------------------
%
% Purpose:
%  Computes azimuth, elevation and partials from local tangent coordinates
%
% Input:
%   s      Topocentric local tangent coordinates (East-North-Zenith frame)
%
% Outputs:
%   A      Azimuth [rad]
%   E      Elevation [rad]
%   dAds   Partials of azimuth w.r.t. s
%   dEds   Partials of elevation w.r.t. s
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------*/


const double pi2 = 2.0 * M_PI;

AzElPaResult AzElPa(const Matrix& s) {
    AzElPaResult result;

    if (s.getFilas() != 3 || s.getColumnas() != 1) {
        throw std::invalid_argument("Input s must be a 3x1 matrix");
    }

    double s1 = s(1,1);
    double s2 = s(2,1);
    double s3 = s(3,1);

    double rho = sqrt(s1*s1 + s2*s2);


    result.Az = atan2(s1, s2);
    if (result.Az < 0.0) {
        result.Az += pi2;
    }


    result.El = atan(s3 / rho);


    result.dAds = Matrix(1, 3);
    result.dAds(1,1) = s2 / (rho*rho);
    result.dAds(1,2) = -s1 / (rho*rho);
    result.dAds(1,3) = 0.0;

    double s_dot_s = s1*s1 + s2*s2 + s3*s3;
    result.dEds = Matrix(1, 3);
    result.dEds(1,1) = (-s1*s3/rho) / s_dot_s;
    result.dEds(1,2) = (-s2*s3/rho) / s_dot_s;
    result.dEds(1,3) = rho / s_dot_s;

    return result;
}

