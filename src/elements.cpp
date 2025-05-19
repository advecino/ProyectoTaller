//
// Created by adria on 10/05/2025.
//

#include <stdexcept>
#include "../include/elements.h"

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


KeplerianElements elements(const Matrix& r, const Matrix& v) {

    const double TWO_PI = 2 * M_PI;

    // 1) angular momentum
    Matrix h = Matrix::cross(r, v);
    double H = h.norm();
    double p = H*H / GM_Earth;

    // check non-circular, non-equatorial
    if (p <= 0) throw std::invalid_argument("elements: degenerate orbit");
    double hx = h(1,1), hy = h(2,1), hz = h(3,1);
    double n2 = hx*hx + hy*hy;
    if (n2 < 1e-12) throw std::invalid_argument("elements: equatorial orbit");

    // 2) node vector and Omega
    double Omega = std::atan2(hx, -hy);
    Omega = fmod(Omega + TWO_PI, TWO_PI);

    // 3) inclination
    double inc = std::atan2(std::sqrt(n2), hz);

    // 4) semimajor axis
    double R = r.norm();
    double V2 = Matrix::dot(v,v);
    double a = 1.0 / (2.0/R - V2/GM_Earth);

    // 5) eccentricity vector
    // e cos E = 1 - R/a ; e sin E = (r·v)/sqrt(GM a)
    double eCosE = 1 - R/a;
    double eSinE = Matrix::dot(r,v)/std::sqrt(GM_Earth*a);
    double e2 = eCosE*eCosE + eSinE*eSinE;
    double ecc = std::sqrt(e2);
    if (ecc < 1e-8) throw std::invalid_argument("elements: circular orbit");

    // 6) eccentric and mean anomalies
    double E = std::atan2(eSinE, eCosE);
    double M = fmod(E - eSinE + TWO_PI, TWO_PI);

    // 7) argument of latitude u and true anomaly nu
    double u = std::atan2(r(3,1)*H, -r(1,1)*hy + r(2,1)*hx);
    double nu = std::atan2(std::sqrt(1-e2)*eSinE, eCosE-e2);

    // 8) argument of pericenter
    double omega = fmod(u - nu + TWO_PI, TWO_PI);

    return { p, a, ecc, inc, Omega, omega, M };
}