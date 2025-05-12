//
// Created by adria on 12/05/2025.
//

#include "Accel.h"
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
/*

Matrix Accel(double x, const Matrix& Y, const Matrix& eopdata)
    {
    double Mjd_UTC = 59000.5;
        // Extraer parámetros de IERS
        IERSResult iersResult = IERS(eopdata, Mjd_UTC + x/86400.0, "l");

        // Diferencias de tiempo
        TimeDiffs timeDiffResult = timediff(iersResult.UT1_UTC, iersResult.TAI_UTC);


        // Cálculos de fechas
        double Mjd_UT1 = Mjd_UTC + x/86400.0 + iersResult.UT1_UTC/86400.0;
        double Mjd_TT = Mjd_UTC + x/86400.0 + timeDiffResult.TT_UTC/86400.0;

        // Matrices de transformación
        Matrix P = PrecMatrix(MJD_J2000, Mjd_TT);
        Matrix N = NutMatrix(Mjd_TT);
        Matrix T = N * P;
        Matrix E = PoleMatrix(iersResult.x_pole, iersResult.y_pole) * GHAMatrix(Mjd_UT1) * T;

        // Efemérides planetarias
        double MJD_TDB = Mjday_TDB(Mjd_TT);
        PlanetaryPositions jplResult = JPL_Eph_DE430(MJD_TDB, PC);


        // Aceleración debido al campo gravitatorio armónico de la Tierra
        Matrix position = Y.getSubMatrix(1,3,1,1);
        Matrix a = AccelHarmonic(position, E, 20, 20);

        // Perturbaciones lunares y solares
        if (auxParam.sun) {
            a = a + AccelPointMass(position, R_Sun, GM_Sun);
        }

        if (auxParam.moon) {
            a = a + AccelPointMass(position, R_Moon, GM_Moon);
        }

        // Perturbaciones planetarias
        if (auxParam.planets) {
            a = a + AccelPointMass(position, r_Mercury, consts.GM_Mercury);
            a = a + AccelPointMass(position, r_Venus, consts.GM_Venus);
            a = a + AccelPointMass(position, r_Mars, consts.GM_Mars);
            a = a + AccelPointMass(position, r_Jupiter, consts.GM_Jupiter);
            a = a + AccelPointMass(position, r_Saturn, consts.GM_Saturn);
            a = a + AccelPointMass(position, r_Uranus, consts.GM_Uranus);
            a = a + AccelPointMass(position, r_Neptune, consts.GM_Neptune);
            a = a + AccelPointMass(position, r_Pluto, consts.GM_Pluto);
        }

        // Construir vector de estado derivado
        Matrix velocity = Y.getSubMatrix(3, 5, 0, 0);
        Matrix dY(6, 1);
        dY.setSubMatrix(0, 0, velocity);
        dY.setSubMatrix(3, 0, a);

        return dY;
    }
*/