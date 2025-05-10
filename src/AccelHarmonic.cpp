#include "../include/AccelHarmonic.h"
#include "../include/global.h"
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


Matrix AccelHarmonic(const Matrix r, Matrix E, int n_max, int m_max) {
    const double r_ref = 6378.1363e3;   //Earth's radius [m]; GGM03S
    const double gm = 398600.4415e9;    // [m^3/s^2]; GGM03S

    if (r.getFilas() != 3 || r.getColumnas() != 1) {
        throw std::runtime_error("El vector 'r' debe ser 3x1");
    }
    if (E.getFilas() != 3 || E.getColumnas() != 3) {
        throw std::runtime_error("La matriz 'E' debe ser 3x3");
    }
    if (n_max < 0 || m_max < 0 || m_max > n_max) {
        throw std::runtime_error("n_max y m_max deben ser >= 0, y m_max <= n_max");
    }

    // Body-fixed position
    Matrix r_bf = E * r;


    if(r_bf.getFilas() < 3 || r_bf.getColumnas() < 1) {
        throw std::runtime_error("r_bf debe ser al menos 3x1");
    }

    double d = r_bf.norm();// distance
    double latgc = asin(r_bf(3,1)/d);
    double lon = atan2(r_bf(2,1), r_bf(1,1));

    Matrix dpnm = Matrix(n_max+1, m_max+1);
    Matrix pnm = Matrix(n_max+1, m_max+1);
    Legendre(n_max, m_max, latgc,pnm,dpnm);

    if (pnm.getFilas() < n_max+1 || pnm.getColumnas() < m_max+1 ||
        dpnm.getFilas() < n_max+1 || dpnm.getColumnas() < m_max+1) {
        throw std::runtime_error("Las matrices pnm/dpnm son demasiado pequeñas");
    }

    double dUdr = 0;
    double dUdlatgc = 0;
    double dUdlon = 0;
    double q3 = 0, q2 = q3, q1 = q2;

    for(int n = 0; n <= n_max; n++) {
        double b1 = (-gm/pow(d,2)) * pow(r_ref/d,n) * (n+1);
        double b2 = (gm/d) * pow(r_ref/d,n);
        double b3 = (gm/d) * pow(r_ref/d,n);

        for(int m = 0; m <= m_max; m++) {
            q1 += pnm(n+1,m+1) * (Cnm[n+1][m+1]*cos(m*lon) + Snm[n+1][m+1]*sin(m*lon));
            q2 += dpnm(n+1,m+1) * (Cnm[n+1][m+1]*cos(m*lon) + Snm[n+1][m+1]*sin(m*lon));
            q3 += m * pnm(n+1,m+1) * (Snm[n+1][m+1]*cos(m*lon) - Cnm[n+1][m+1]*sin(m*lon));
        }

        dUdr     = dUdr     + q1*b1;
        dUdlatgc = dUdlatgc + q2*b2;
        dUdlon   = dUdlon   + q3*b3;
        q3 = 0; q2 = q3; q1 = q2;
    }

    // Body-fixed acceleration
    double r2xy = r_bf(1, 1) * r_bf(1,1) + r_bf(2, 1) * r_bf(2, 1);


    double ax = (1/d*dUdr - r_bf(3,1)/(pow(d,2)*sqrt(r2xy))*dUdlatgc)*r_bf(1,1) - (1/r2xy*dUdlon)*r_bf(2,1);
    double ay = (1/d*dUdr - r_bf(3,1)/(pow(d,2)*sqrt(r2xy))*dUdlatgc)*r_bf(2,1) + (1/r2xy*dUdlon)*r_bf(1,1);
    double az = 1/d*dUdr*r_bf(3,1) + sqrt(r2xy)/pow(d,2)*dUdlatgc;


    Matrix a_bf(3, 1);
    a_bf(1,1) = ax;
    a_bf(2,1) = ay;
    a_bf(3,1) = az;

    // Inertial acceleration
    Matrix a_final(E.transpuesta());
    Matrix a_fin = a_final*a_bf;
    return a_fin;
}