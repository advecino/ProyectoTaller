#include "../include/Accel.h"
#include "../include/IERS.h"
#include "../include/TimeDiff.h"
#include "../include/Mjday_TDB.h"
#include "../include/JPL_Eph_DE430.h"
#include "../include/AccelHarmonic.h"
#include "../include/AccelPointMass.h"
#include "../include/PrecMatrix.h"
#include "../include/NutMatrix.h"
#include "../include/PoleMatrix.h"
#include "../include/GHAMatrix.h"
#include "../include/Sat_const.h"
#include "../include/AuxParam.h"
#include "../include/global.h"
#include <stdexcept>

/*
 * %--------------------------------------------------------------------------
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


Matrix Accel(
        double x,
        const Matrix& Y,
        AuxParam& params,
        Matrix& eop
) {
    if (Y.getFilas() != 6 || Y.getColumnas() != 1) {
        throw std::invalid_argument("Accel: Y must be 6×1");
    }

    double Mjd_UTC = params.Mjd_UTC + x/86400.0;
    IERSResult ires = IERS(eop, Mjd_UTC, 'l');
    double x_pole  = ires.x_pole;
    double y_pole  = ires.y_pole;
    double UT1_UTC = ires.UT1_UTC;
    double TAI_UTC = ires.TAI_UTC;

    TimeDiffs td = timediff(UT1_UTC, TAI_UTC);
    double TT_UTC = td.TT_UTC;

    double Mjd_UT1 = params.Mjd_UTC + x/86400.0 + UT1_UTC/86400.0;
    double Mjd_TT  = params.Mjd_UTC + x/86400.0 + TT_UTC /86400.0;

    Matrix P = PrecMatrix(MJD_J2000, Mjd_TT);
    Matrix N = NutMatrix(    Mjd_TT);
    Matrix Tmat = N * P;
    Matrix E = PoleMatrix(x_pole, y_pole)
               * GHAMatrix(Mjd_UT1)
               * Tmat;

    Matrix r(3,1), v(3,1);
    for (int i = 1; i <= 3; ++i) {
        r(i,1) = Y(i,   1);
        v(i,1) = Y(i+3, 1);
    }

    Matrix a = AccelHarmonic(r, E, params.n, params.m);

    if (params.sun || params.moon || params.planets) {
        double Mjd_TDB = Mjday_TDB(Mjd_TT);
        PlanetaryPositions jr = JPL_Eph_DE430(Mjd_TDB);
        if (params.sun)   a = a + AccelPointMass(r, jr.r_Sun,     GM_Sun);
        if (params.moon)  a = a + AccelPointMass(r, jr.r_Moon,    GM_Moon);
        if (params.planets) {
            a = a + AccelPointMass(r, jr.r_Mercury, GM_Mercury);
            a = a + AccelPointMass(r, jr.r_Venus,   GM_Venus);
            a = a + AccelPointMass(r, jr.r_Mars,    GM_Mars);
            a = a + AccelPointMass(r, jr.r_Jupiter, GM_Jupiter);
            a = a + AccelPointMass(r, jr.r_Saturn,  GM_Saturn);
            a = a + AccelPointMass(r, jr.r_Uranus,  GM_Uranus);
            a = a + AccelPointMass(r, jr.r_Neptune, GM_Neptune);
            a = a + AccelPointMass(r, jr.r_Pluto,   GM_Pluto);
        }
    }

    Matrix dY(6,1);
    for (int i = 1; i <= 3; ++i) {
        dY(i,   1) = v(i,1);
        dY(i+3, 1) = a(i,1);
    }
    return dY;
}








