#include "../include/AccelHarmonic.h"
#include "../include/global.h"
#include "../include/Legendre.h"
#include <cmath>
#include <stdexcept>
#include <iostream>

/*%--------------------------------------------------------------------------
%
% AccelHarmonic.m
%
% Purpose:
%   Computes the acceleration due to the harmonic gravity field of the
%   central body
%
% Inputs:
%   r           Satellite position vector in the inertial system
%   E           Transformation matrix to body-fixed system
%   n_max       Maximum degree
%   m_max       Maximum order (m_max<=n_max; m_max=0 for zonals, only)
%
% Output:
%   a           Acceleration (a=d^2r/dt^2)
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------*/


Matrix AccelHarmonic(Matrix r, Matrix E, int n_max, int m_max) {
    const double r_ref = 6378.1363e3;   // Earth's radius [m]
    const double gm    = 398600.4415e9; // Standard gravitational parameter [m^3/s^2]

    if (r.getFilas() != 3 || r.getColumnas() != 1)
        throw std::runtime_error("El vector 'r' debe ser 3x1");
    if (E.getFilas() != 3 || E.getColumnas() != 3)
        throw std::runtime_error("La matriz 'E' debe ser 3x3");
    if (n_max < 0 || m_max < 0 || m_max > n_max)
        throw std::runtime_error("Valores de n_max/m_max inválidos");

    // Transform to body-fixed frame
    Matrix r_bf = E * r;
    double d   = r_bf.norm();
    double lat = std::asin(r_bf(3,1) / d);
    double lon = std::atan2(r_bf(2,1), r_bf(1,1));

    // Prepare matrices pnm and dpnm with correct dimensions
    Matrix pnm(n_max+2, m_max+2);
    Matrix dpnm(n_max+2, m_max+2);
    Legendre(n_max, m_max, lat, pnm, dpnm);
    std::cout<<"hola"<<std::endl;
    pnm.print();
    dpnm.print();

    // Compute partial derivatives of potential
    double dUdr = 0.0, dUdlat = 0.0, dUdlon = 0.0;
    for (int n = 0; n <= n_max; ++n) {
        double r_ratio_n = std::pow(r_ref / d, n);
        double b1 = -gm * (n+1) * r_ratio_n / (d * d);
        double b2 =  gm * r_ratio_n / d;
        for (int m = 0; m <= m_max; ++m) {
            double C = Cnm[n+1][m+1];
            double S = Snm[n+1][m+1];
            double cosml = std::cos(m * lon);
            double sinml = std::sin(m * lon);
            double P = pnm(n+1, m+1);
            double dP = dpnm(n+1, m+1);
            dUdr   += P  * (C * cosml + S * sinml) * b1;
            dUdlat += dP * (C * cosml + S * sinml) * b2;
            dUdlon += m  * P  * (S * cosml - C * sinml) * b2;
        }
    }

    // Body-fixed acceleration
    double x = r_bf(1,1), y = r_bf(2,1), z = r_bf(3,1);
    double r2xy = x*x + y*y;

    double ax = ((dUdr / d) - (z * dUdlat) / (d*d * std::sqrt(r2xy))) * x
                - (dUdlon / r2xy) * y;
    double ay = ((dUdr / d) - (z * dUdlat) / (d*d * std::sqrt(r2xy))) * y
                + (dUdlon / r2xy) * x;
    double az = (dUdr / d) * z
                + (std::sqrt(r2xy) / (d*d)) * dUdlat;

    Matrix a_bf(3,1);
    a_bf(1,1) = ax;
    a_bf(2,1) = ay;
    a_bf(3,1) = az;

    // Convert back to inertial frame
    return E.transpuesta() * a_bf;
}

