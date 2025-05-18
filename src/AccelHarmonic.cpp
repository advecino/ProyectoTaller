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


Matrix AccelHarmonic(
        Matrix& r,
        Matrix& E,
        int n_max,
        int m_max
) {
    // Physical constants
    const double r_ref = 6378.1363e3;   // [m]
    const double gm    = 398600.4415e9; // [m^3/s^2]

    // Sanity checks
    if (r.getFilas()!=3 || r.getColumnas()!=1)
        throw std::invalid_argument("r must be 3×1");
    if (E.getFilas()!=3 || E.getColumnas()!=3)
        throw std::invalid_argument("E must be 3×3");
    if (n_max<0 || m_max<0 || m_max>n_max)
        throw std::invalid_argument("Invalid n_max/m_max");

    // Body‐fixed position
    Matrix r_bf = E * r;

    double d = r_bf.norm();
    if (d < 1e-16) {
        return Matrix(3,1);  // aceleración cero
    }

    if (n_max==0 && m_max==0) {
        double inv_r3 = -gm/(d*d*d);
        Matrix a_bf(3,1);
        for(int i=1;i<=3;++i) a_bf(i,1) = r_bf(i,1)*inv_r3;
        return E.transpuesta() * a_bf;
    }
    // Compute geocentric latitude and longitude
    double x = r_bf(1,1), y = r_bf(2,1), zc = r_bf(3,1);
    double latgc = std::asin(zc/d);
    double lon   = std::atan2(y, x);

    // Compute associated Legendre functions and their derivatives
    Matrix pnm(0,0), dpnm(0,0);
    Legendre(n_max, m_max, latgc, pnm, dpnm);

    // Accumulate potential derivatives
    double dUdr = 0.0, dUdlat = 0.0, dUdlon = 0.0;
    for (int n = 0; n <= n_max; ++n) {
        double factor = std::pow(r_ref/d, n);
        double b1 = -gm/(d*d) * factor * (n+1);
        double b2 =  gm/d       * factor;
        double b3 =  b2;
        double q1 = 0.0, q2 = 0.0, q3 = 0.0;
        for (int m = 0; m <= m_max; ++m) {
            double C = Cnm[n+1][m+1];
            double S = Snm[n+1][m+1];
            double P = pnm(n+1,m+1);
            double dP= dpnm(n+1,m+1);
            double cosml = std::cos(m*lon);
            double sinml = std::sin(m*lon);
            q1 += P * (C*cosml + S*sinml);
            q2 += dP* (C*cosml + S*sinml);
            q3 += m * P * (S*cosml - C*sinml);
        }
        dUdr     += q1 * b1;
        dUdlat   += q2 * b2;
        dUdlon   += q3 * b3;
    }

    // Body‐fixed acceleration components
    double r2xy = x*x + y*y;
    double common = (1.0/d)*dUdr - (zc/(d*d*std::sqrt(r2xy)))*dUdlat;
    double ax_bf = common * x - (1.0/r2xy*dUdlon) * y;
    double ay_bf = common * y + (1.0/r2xy*dUdlon) * x;
    double az_bf = (1.0/d)*dUdr * zc + (std::sqrt(r2xy)/(d*d)) * dUdlat;

    // Assemble body‐fixed acceleration vector
    Matrix a_bf(3,1);
    a_bf(1,1)=ax_bf;
    a_bf(2,1)=ay_bf;
    a_bf(3,1)=az_bf;


    // Transform back to inertial frame: a = E' * a_bf
    return E.transpuesta() * a_bf;
}

