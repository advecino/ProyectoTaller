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

    Matrix h = Matrix::cross(r, v);
    double H = h.norm();
    double p = H*H / GM_Earth;

    if (p <= 0) throw std::invalid_argument("elements: degenerate orbit");
    double hx = h(1,1), hy = h(2,1), hz = h(3,1);
    double n2 = hx*hx + hy*hy;
    if (n2 < 1e-12) throw std::invalid_argument("elements: equatorial orbit");

    double Omega = std::atan2(hx, -hy);
    Omega = fmod(Omega + TWO_PI, TWO_PI);

    double inc = std::atan2(std::sqrt(n2), hz);

    double R = r.norm();
    double V2 = Matrix::dot(v,v);
    double a = 1.0 / (2.0/R - V2/GM_Earth);


    double eCosE = 1 - R/a;
    double eSinE = Matrix::dot(r,v)/std::sqrt(GM_Earth*a);
    double e2 = eCosE*eCosE + eSinE*eSinE;
    double ecc = std::sqrt(e2);
    if (ecc < 1e-8) throw std::invalid_argument("elements: circular orbit");

    double E = std::atan2(eSinE, eCosE);
    double M = fmod(E - eSinE + TWO_PI, TWO_PI);

    double u = std::atan2(r(3,1)*H, -r(1,1)*hy + r(2,1)*hx);
    double nu = std::atan2(std::sqrt(1-e2)*eSinE, eCosE-e2);

    double omega = fmod(u - nu + TWO_PI, TWO_PI);

    return { p, a, ecc, inc, Omega, omega, M };
}