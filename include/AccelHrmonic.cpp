#include "AccelHrmonic.h"
#include "global.h"
#include <cmath>

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



Matrix AccelHarmonic(const Matrix& r, Matrix E, int n_max, int m_max) {
    const double r_ref = 6378.1363e3;   //Earth's radius [m]; GGM03S
    const double gm = 398600.4415e9;    // [m^3/s^2]; GGM03S

    // Body-fixed position
    Matrix r_bf = E.operator*(r);

    // Cantidades auxiliares
    double d = r_bf.norm();                     // distance
    double latgc = asin(r_bf(2,0)/d);
    double lon = atan2(r_bf(1,0), r_bf(0,0));

    Matrix dpnm = Matrix(0, 0);
    Matrix pnm = Matrix(0, 0);
    Legendre(n_max, m_max, latgc,pnm,dpnm);

    double dUdr = 0;
    double dUdlatgc = 0;
    double dUdlon = 0;
    double q3 = 0, q2 = q3, q1 = q2;

    for(int n = 0; n <= n_max; n++) {
        double b1 = (-gm/pow(d,2)) * pow(r_ref/d,n) * (n+1);
        double b2 = (gm/d) * pow(r_ref/d,n);
        double b3 = (gm/d) * pow(r_ref/d,n);;

        for(int m = 0; m <= m_max; m++) {
            q1 += pnm(n+1,m+1) * (Cnm[n+1][m+1]*cos(m*lon) + Snm[n+1][m+1]*sin(m*lon));
            q2 += dpnm(n+1,m+1) * (Cnm[n+1][m+1]*cos(m*lon) + Snm[n+1][m+1]*sin(m*lon));
            q3 += m * pnm(n+1,m+1) * (Snm[n+1][m+1]*cos(m*lon) - Cnm[n+1][m+1]*sin(m*lon));
        }

        dUdr += q1 * b1;
        dUdlatgc += q2 * b2;
        dUdlon += q3 * b3;

        q3 = q2 = q1 = 0;
    }

    // Body-fixed acceleration
    double r2xy = pow(r_bf(0,0), 2) + pow(r_bf(1,0), 2);

    double ax = (1/d*dUdr - r_bf(2,0)/(pow(d,2)*sqrt(r2xy))*dUdlatgc)*r_bf(0,0) - (1/r2xy*dUdlon)*r_bf(1,0);
    double ay = (1/d*dUdr - r_bf(2,0)/(pow(d,2)*sqrt(r2xy))*dUdlatgc)*r_bf(1,0) + (1/r2xy*dUdlon)*r_bf(0,0);
    double az = 1/d*dUdr*r_bf(2,0) + sqrt(r2xy)/pow(d,2)*dUdlatgc;

    Matrix a_bf(3, 1);
    a_bf(0,0) = ax;
    a_bf(1,0) = ay;
    a_bf(2,0) = az;

    // Inertial acceleration
    return E.transpuesta() * a_bf;
}



