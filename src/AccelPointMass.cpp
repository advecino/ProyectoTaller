#include "../include/AccelPointMass.h"
#include <cmath>
/*
%--------------------------------------------------------------------------
%
% AccelPointMass: Computes the perturbational acceleration due to a point
%				  mass
%
% Inputs:
%   r           Satellite position vector
%   s           Point mass position vector
%   GM          Gravitational coefficient of point mass
%
% Output:
%   a    		Acceleration (a=d^2r/dt^2)
%
% Last modified:   2018/01/27   M. Mahooti
%
%--------------------------------------------------------------------------*/

Matrix AccelPointMass( Matrix& r, const Matrix& s, double GM)
{
    // Relative position vector of satellite w.r.t. point mass
    Matrix d = r - s;

    double norm_d = 0.0;
    double norm_s = 0.0;

    for(int i = 1; i <= 3; ++i) {
        norm_d += d(i,1) * d(i,1);
        norm_s += s(i,1) * s(i,1);
    }

    norm_d = sqrt(norm_d);
    norm_s = sqrt(norm_s);

    double term1 = 1.0 / (norm_d * norm_d * norm_d);
    double term2 = 1.0 / (norm_s * norm_s * norm_s);

    Matrix a(3, 1);
    for(int i = 1; i <= 3; ++i) {
        a(i,1) = -GM * (d(i,1) * term1 + s(i,1) * term2);
    }

    return a;
}
