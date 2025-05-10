//
// Created by adria on 10/05/2025.
//

#include "elements.h"
/*
%--------------------------------------------------------------------------
%
% Purpose:
%   Computes the osculating Keplerian elements from the satellite state
%   vector for elliptic orbits
%
% Input:
%    y        State vector (x,y,z,vx,vy,vz)
%
% Outputs:
%    p        semilatus rectum [m]
%    a        Semimajor axis
%    e        Eccentricity
%    i        Inclination [rad]
%    Omega    Longitude of the ascending node [rad]
%    omega    Argument of pericenter [rad]
%    M        Mean anomaly [rad]
%
% Notes:
%   The function cannot be used with state vectors describing a circular
%   or non-inclined orbit.
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------*/
KeplerianElements elements(const Matrix& y) {
    const double pi2 = 2.0 * M_PI;
    KeplerianElements el;

    // Extract position and velocity from state vector
    Matrix r = y.getSubMatrix(1, 3, 1, 1);  // First 3 elements are position
    Matrix v = y.getSubMatrix(4, 6, 1, 1);  // Next 3 elements are velocity

    // Compute areal velocity (h = r × v)
    Matrix h = Matrix::cross(r, v);
    double magh = h.norm();
    el.p = magh * magh / GM_Earth;

    // Longitude of ascending node (Ω)
    el.Omega = atan2(h(1,1), -h(2,1));
    el.Omega = fmod(el.Omega, pi2);

    // Inclination (i)
    el.i = atan2(sqrt(h(1,1)*h(1,1) + h(2,1)*h(2,1)), h(3,1));

    // Argument of latitude (u)
    double H = h.norm();
    double u = atan2(r(3,1)*H, -r(1,1)*h(2,1) + r(2,1)*h(1,1));

    // Distance (R)
    double R = r.norm();

    // Semi-major axis (a)
    el.a = 1.0 / (2.0/R - Matrix::dot(v,v)/GM_Earth);

    // Eccentricity (e)
    double eCosE = 1.0 - R/el.a;
    double eSinE = Matrix::dot(r,v)/sqrt(GM_Earth * el.a);
    double e2 = eCosE*eCosE + eSinE*eSinE;
    el.e = sqrt(e2);

    // Eccentric anomaly (E)
    double E = atan2(eSinE, eCosE);

    // Mean anomaly (M)
    el.M = fmod(E - eSinE, pi2);

    // True anomaly (ν)
    double nu = atan2(sqrt(1.0 - e2)*eSinE, eCosE - e2);

    // Argument of pericenter (ω)
    el.omega = fmod(u - nu, pi2);

    return el;
}