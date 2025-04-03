//
// Created by advecino on 03/04/2025.
//

#include "../include/Accel.h"
/*
%--------------------------------------------------------------------------
%
% Accel.m
%
% Purpose:
%   Computes the acceleration of an Earth orbiting satellite due to
%    - the Earth's harmonic gravity field,
%    - the gravitational perturbations of the Sun and Moon
%    - the solar radiation pressure and
%    - the atmospheric drag
%
% Inputs:
%   Mjd_TT      Terrestrial Time (Modified Julian Date)
%   Y           Satellite state vector in the ICRF/EME2000 system
%
% Output:
%   dY		    Acceleration (a=d^2r/dt^2) in the ICRF/EME2000 system
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------*/
#include <cmath>
#include <iostream>


Constants consts;
AuxParams AuxParam;
std::vector<double> eopdata;

std::vector<double> Accel(double x, const std::vector<double>& Y) {
    std::vector<double> dY(6);

    double mjd = AuxParam.Mjd_UTC + x / 86400.0;
    auto [x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC] = IERS(eopdata, mjd, 'l');
    auto [UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC] = timediff(UT1_UTC, TAI_UTC);

    double Mjd_UT1 = mjd + UT1_UTC / 86400.0;
    double Mjd_TT = mjd + TT_UTC / 86400.0;


    std::vector<double> a = AccelHarmonic({Y[0], Y[1], Y[2]}, {}, AuxParam.n, AuxParam.m);


    if (AuxParam.sun) {
        std::vector<double> r_Sun; // Debe ser obtenido de JPL_Eph_DE430
        a = AccelPointMass({Y[0], Y[1], Y[2]}, r_Sun, consts.GM_Sun);
    }
    if (AuxParam.moon) {
        std::vector<double> r_Moon; // Debe ser obtenido de JPL_Eph_DE430
        a = AccelPointMass({Y[0], Y[1], Y[2]}, r_Moon, consts.GM_Moon);
    }
    if (AuxParam.planets) {
        std::vector<std::vector<double>> r_planets; // Debe ser obtenido de JPL_Eph_DE430
        std::vector<double> GM_planets = {consts.GM_Mercury, consts.GM_Venus, consts.GM_Mars,
                                          consts.GM_Jupiter, consts.GM_Saturn, consts.GM_Uranus,
                                          consts.GM_Neptune, consts.GM_Pluto};
        for (size_t i = 0; i < r_planets.size(); i++) {
            a = AccelPointMass({Y[0], Y[1], Y[2]}, r_planets[i], GM_planets[i]);
        }
    }

    dY[0] = Y[3];
    dY[1] = Y[4];
    dY[2] = Y[5];
    dY[3] = a[0];
    dY[4] = a[1];
    dY[5] = a[2];

    return dY;

}



