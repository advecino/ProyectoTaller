

#include "Accel.h"
#include "IERS.h"
#include "TimeDiff.h"
#include "Mjday_TDB.h"
#include "JPL_Eph_DE430.h"
#include "AccelHarmonic.h"
#include "AccelPointMass.h"
#include "PrecMatrix.h"
#include "NutMatrix.h"
#include "PoleMatrix.h"
#include "GHAMatrix.h"
#include "Sat_const.h"
#include "AuxParam.h"
#include "global.h"

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



Matrix Accel(double x, const Matrix& Y) {
    double Mjd_UTC = AuxParamGlob.Mjd_UTC + x/86400.0;
    auto iers = IERS(eopdata, Mjd_UTC, 'l');

    auto td   = timediff(iers.UT1_UTC, iers.TAI_UTC);
    double Mjd_UT1 = Mjd_UTC + iers.UT1_UTC/86400.0;
    double Mjd_TT  = Mjd_UTC + td.TT_UTC   /86400.0;

    // 2) Construir matriz E de transformación a cuerpo fijo
    Matrix P = PrecMatrix(MJD_J2000, Mjd_TT);
    Matrix N = NutMatrix(Mjd_TT);
    Matrix T = N * P;
    Matrix E = PoleMatrix(iers.x_pole, iers.y_pole)
               * GHAMatrix(Mjd_UT1)
               * T;

    // 3) Posiciones planetarias
    double Mjd_TDB = Mjday_TDB(Mjd_TT);
    auto pp = JPL_Eph_DE430(Mjd_TDB);

    // 4) Aceleración armónica
    Matrix r = Y.getSubMatrix(1,3,1,1);
    Matrix a = AccelHarmonic(r, E, AuxParamGlob.n, AuxParamGlob.m);

    // 5) Perturbaciones luni–solares
    if (AuxParamGlob.sun)  a = a + AccelPointMass(r, pp.r_Sun,   GM_Sun);
    if (AuxParamGlob.moon) a = a + AccelPointMass(r, pp.r_Moon,  GM_Moon);

    // 6) Perturbaciones planetarias
    if (AuxParamGlob.planets) {
        a = a + AccelPointMass(r, pp.r_Mercury, GM_Mercury);
        a = a + AccelPointMass(r, pp.r_Venus,   GM_Venus);
        a = a + AccelPointMass(r, pp.r_Mars,    GM_Mars);
        a = a + AccelPointMass(r, pp.r_Jupiter, GM_Jupiter);
        a = a + AccelPointMass(r, pp.r_Saturn,  GM_Saturn);
        a = a + AccelPointMass(r, pp.r_Uranus,  GM_Uranus);
        a = a + AccelPointMass(r, pp.r_Neptune, GM_Neptune);
        a = a + AccelPointMass(r, pp.r_Pluto,   GM_Pluto);
    }

    // 7) Devolver [v; a]
    Matrix dY(6,1);
    for (int i = 1; i <= 3; ++i) {
        dY(i,1)   = Y(i+3,1);  // componente velocidad
        dY(i+3,1) = a(i,1);    // componente aceleración
    }
    return dY;
}