#include "JPL_Eph_DE430.h"

/*
%--------------------------------------------------------------------------
%
% JPL_Eph_DE430: Computes the sun, moon, and nine major planets' equatorial
%                position using JPL Ephemerides
%
% Inputs:
%   Mjd_TDB         Modified julian date of TDB
%
% Output:
%   r_Earth(solar system barycenter (SSB)),r_Mars,r_Mercury,r_Venus,
%   r_Jupiter,r_Saturn,r_Uranus,r_Neptune,r_Pluto,r_Moon,
%   r_Sun(geocentric equatorial position ([m]) referred to the
%   International Celestial Reference Frame (ICRF))
%
% Notes: Light-time is already taken into account
%
% Last modified:   2018/01/11   M. Mahooti
%
%--------------------------------------------------------------------------*/


PlanetaryPositions JPL_Eph_DE430(double Mjd_TDB,  const Matrix& PC) {
    PlanetaryPositions positions;


    double JD = Mjd_TDB + 2400000.5;

    // Find the appropriate interval in the ephemeris data
    int i = 1;
    for (int col = 1; col <= PC.getColumnas(); ++col) {
        if (PC(1, col) <= JD && JD <= PC(2, col)) {
            i = col;
            break;
        }
    }

    if (i == -1) {
        throw std::runtime_error("Could not find MJD in PC data");
    }

    Matrix PCtemp = PC.getColumn(i);

    double t1 = PCtemp(1,1) - 2400000.5; // MJD at start of interval
    double dt = Mjd_TDB - t1;

    // Earth calculation
    int temp_start = 231;
    Matrix Cx_Earth = PCtemp.getSubMatrix(temp_start, temp_start + 12, 1, 1);
    Matrix Cy_Earth = PCtemp.getSubMatrix(temp_start + 13, temp_start + 25, 1, 1);
    Matrix Cz_Earth = PCtemp.getSubMatrix(temp_start + 26, temp_start + 38, 1, 1);

    temp_start += 39;
    Matrix Cx = PCtemp.getSubMatrix(temp_start, temp_start + 12, 1,1);
    Matrix Cy = PCtemp.getSubMatrix(temp_start + 13, temp_start + 25, 1,1);
    Matrix Cz = PCtemp.getSubMatrix(temp_start + 26, temp_start + 38, 1,1);

    Cx_Earth = Cx_Earth.concatenate(Cx, 0);
    Cy_Earth = Cy_Earth.concatenate(Cy, 0);
    Cz_Earth = Cz_Earth.concatenate(Cz, 0);

    int j = 0;
    double Mjd0 = t1;
    if (dt > 16 && dt <= 32) {
        j = 1;
        Mjd0 = t1 + 16 * j;
    }

    positions.r_Earth = Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0 + 16,
                               Cx_Earth.getSubMatrix(13*j, 13*j + 12, 1, 1),
                               Cy_Earth.getSubMatrix(13*j, 13*j + 12, 1, 1),
                               Cz_Earth.getSubMatrix(13*j, 13*j + 12, 1, 1)) * 1e3;

    // Moon calculation
    temp_start = 441 - 1;
    Matrix Cx_Moon = PCtemp.getSubMatrix(temp_start, temp_start + 12, 1, 1);
    Matrix Cy_Moon = PCtemp.getSubMatrix(temp_start + 13, temp_start + 25, 1, 1);
    Matrix Cz_Moon = PCtemp.getSubMatrix(temp_start + 26, temp_start + 38, 1, 1);

    for (int k = 1; k <= 7; k++) {
        temp_start += 39;
        Cx = PCtemp.getSubMatrix(temp_start, temp_start + 12, 1, 1);
        Cy = PCtemp.getSubMatrix(temp_start + 13, temp_start + 25, 1, 1);
        Cz = PCtemp.getSubMatrix(temp_start + 26, temp_start + 38, 1, 1);

        Cx_Moon = Cx_Moon.concatenate(Cx, 0);
        Cy_Moon = Cy_Moon.concatenate(Cy, 0);
        Cz_Moon = Cz_Moon.concatenate(Cz, 0);
    }

    j = static_cast<int>(dt / 4);
    if (j > 7) j = 7;
    Mjd0 = t1 + 4 * j;

    positions.r_Moon = Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0 + 4,
                              Cx_Moon.getSubMatrix(13*j, 13*j + 12, 1, 1),
                              Cy_Moon.getSubMatrix(13*j, 13*j + 12, 1, 1),
                              Cz_Moon.getSubMatrix(13*j, 13*j + 12, 1, 1)) * 1e3;

    // Sun calculation
    temp_start = 753 - 1;
    Matrix Cx_Sun = PCtemp.getSubMatrix(temp_start, temp_start + 10, 1, 1);
    Matrix Cy_Sun = PCtemp.getSubMatrix(temp_start + 11, temp_start + 21, 1, 1);
    Matrix Cz_Sun = PCtemp.getSubMatrix(temp_start + 22, temp_start + 32, 1, 1);

    temp_start += 33;
    Cx = PCtemp.getSubMatrix(temp_start, temp_start + 10, 1, 1);
    Cy = PCtemp.getSubMatrix(temp_start + 11, temp_start + 21, 1, 1);
    Cz = PCtemp.getSubMatrix(temp_start + 22, temp_start + 32, 1, 1);

    Cx_Sun = Cx_Sun.concatenate(Cx, 0);
    Cy_Sun = Cy_Sun.concatenate(Cy, 0);
    Cz_Sun = Cz_Sun.concatenate(Cz, 0);

    j = (dt > 16) ? 1 : 0;
    Mjd0 = t1 + 16 * j;

    positions.r_Sun = Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0 + 16,
                             Cx_Sun.getSubMatrix(11*j, 11*j + 10, 1, 1),
                             Cy_Sun.getSubMatrix(11*j, 11*j + 10, 1, 1),
                             Cz_Sun.getSubMatrix(11*j, 11*j + 10, 1, 1)) * 1e3;

    // Mercury calculation
    temp_start = 3 - 1;
    Matrix Cx_Mercury = PCtemp.getSubMatrix(temp_start, temp_start + 13, 1, 1);
    Matrix Cy_Mercury = PCtemp.getSubMatrix(temp_start + 14, temp_start + 27, 1, 1);
    Matrix Cz_Mercury = PCtemp.getSubMatrix(temp_start + 28, temp_start + 41, 1, 1);

    for (int k = 1; k <= 3; k++) {
        temp_start += 42;
        Cx = PCtemp.getSubMatrix(temp_start, temp_start + 13, 1, 1);
        Cy = PCtemp.getSubMatrix(temp_start + 14, temp_start + 27, 1, 1);
        Cz = PCtemp.getSubMatrix(temp_start + 28, temp_start + 41, 1, 1);

        Cx_Mercury = Cx_Mercury.concatenate(Cx, 0);
        Cy_Mercury = Cy_Mercury.concatenate(Cy, 0);
        Cz_Mercury = Cz_Mercury.concatenate(Cz, 0);
    }

    j = static_cast<int>(dt / 8);
    if (j > 3) j = 3;
    Mjd0 = t1 + 8 * j;

    positions.r_Mercury = Cheb3D(Mjd_TDB, 14, Mjd0, Mjd0 + 8,
                                 Cx_Mercury.getSubMatrix(14*j, 14*j + 13, 1, 1),
                                 Cy_Mercury.getSubMatrix(14*j, 14*j + 13, 1, 1),
                                 Cz_Mercury.getSubMatrix(14*j, 14*j + 13, 1, 1)) * 1e3;

    // Venus calculation
    temp_start = 171 - 1;
    Matrix Cx_Venus = PCtemp.getSubMatrix(temp_start, temp_start + 9, 1, 1);
    Matrix Cy_Venus = PCtemp.getSubMatrix(temp_start + 10, temp_start + 19, 1, 1);
    Matrix Cz_Venus = PCtemp.getSubMatrix(temp_start + 20, temp_start + 29, 1, 1);

    temp_start += 30;
    Cx = PCtemp.getSubMatrix(temp_start, temp_start + 9, 1, 1);
    Cy = PCtemp.getSubMatrix(temp_start + 10, temp_start + 19, 1, 1);
    Cz = PCtemp.getSubMatrix(temp_start + 20, temp_start + 29, 1, 1);

    Cx_Venus = Cx_Venus.concatenate(Cx, 0);
    Cy_Venus = Cy_Venus.concatenate(Cy, 0);
    Cz_Venus = Cz_Venus.concatenate(Cz, 0);

    j = (dt > 16) ? 1 : 0;
    Mjd0 = t1 + 16 * j;

    positions.r_Venus = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0 + 16,
                               Cx_Venus.getSubMatrix(10*j, 10*j + 9, 1, 1),
                               Cy_Venus.getSubMatrix(10*j, 10*j + 9, 1, 1),
                               Cz_Venus.getSubMatrix(10*j, 10*j + 9, 1, 1)) * 1e3;

    // Mars calculation
    temp_start = 309 - 1;
    Matrix Cx_Mars = PCtemp.getSubMatrix(temp_start, temp_start + 10, 1, 1);
    Matrix Cy_Mars = PCtemp.getSubMatrix(temp_start + 11, temp_start + 21, 1, 1);
    Matrix Cz_Mars = PCtemp.getSubMatrix(temp_start + 22, temp_start + 32, 1, 1);

    j = 0;
    Mjd0 = t1;

    positions.r_Mars = Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0 + 32,
                              Cx_Mars.getSubMatrix(11*j, 11*j + 10, 1, 1),
                              Cy_Mars.getSubMatrix(11*j, 11*j + 10, 1, 1),
                              Cz_Mars.getSubMatrix(11*j, 11*j + 10, 1, 1)) * 1e3;

    // Jupiter calculation
    temp_start = 342 - 1;
    Matrix Cx_Jupiter = PCtemp.getSubMatrix(temp_start, temp_start + 7, 1, 1);
    Matrix Cy_Jupiter = PCtemp.getSubMatrix(temp_start + 8, temp_start + 15, 1, 1);
    Matrix Cz_Jupiter = PCtemp.getSubMatrix(temp_start + 16, temp_start + 23, 1, 1);

    j = 0;
    Mjd0 = t1;

    positions.r_Jupiter = Cheb3D(Mjd_TDB, 8, Mjd0, Mjd0 + 32,
                                 Cx_Jupiter.getSubMatrix(8*j, 8*j + 7, 1, 1),
                                 Cy_Jupiter.getSubMatrix(8*j, 8*j + 7, 1, 1),
                                 Cz_Jupiter.getSubMatrix(8*j, 8*j + 7, 1, 1)) * 1e3;

    // Saturn calculation
    temp_start = 366 - 1;
    Matrix Cx_Saturn = PCtemp.getSubMatrix(temp_start, temp_start + 6, 1, 1);
    Matrix Cy_Saturn = PCtemp.getSubMatrix(temp_start + 7, temp_start + 13, 1, 1);
    Matrix Cz_Saturn = PCtemp.getSubMatrix(temp_start + 14, temp_start + 20, 1, 1);

    j = 0;
    Mjd0 = t1;

    positions.r_Saturn = Cheb3D(Mjd_TDB, 7, Mjd0, Mjd0 + 32,
                                Cx_Saturn.getSubMatrix(7*j, 7*j + 6, 1, 1),
                                Cy_Saturn.getSubMatrix(7*j, 7*j + 6, 1, 1),
                                Cz_Saturn.getSubMatrix(7*j, 7*j + 6, 1, 1)) * 1e3;

    // Uranus calculation
    temp_start = 387 - 1;
    Matrix Cx_Uranus = PCtemp.getSubMatrix(temp_start, temp_start + 5, 1, 1);
    Matrix Cy_Uranus = PCtemp.getSubMatrix(temp_start + 6, temp_start + 11, 1, 1);
    Matrix Cz_Uranus = PCtemp.getSubMatrix(temp_start + 12, temp_start + 17, 1, 1);

    j = 0;
    Mjd0 = t1;

    positions.r_Uranus = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0 + 32,
                                Cx_Uranus.getSubMatrix(6*j, 6*j + 5, 1, 1),
                                Cy_Uranus.getSubMatrix(6*j, 6*j + 5, 1, 1),
                                Cz_Uranus.getSubMatrix(6*j, 6*j + 5, 1, 1)) * 1e3;

    // Neptune calculation
    temp_start = 405 - 1;
    Matrix Cx_Neptune = PCtemp.getSubMatrix(temp_start, temp_start + 5, 1, 1);
    Matrix Cy_Neptune = PCtemp.getSubMatrix(temp_start + 6, temp_start + 11, 1, 1);
    Matrix Cz_Neptune = PCtemp.getSubMatrix(temp_start + 12, temp_start + 17, 1, 1);

    j = 0;
    Mjd0 = t1;

    positions.r_Neptune = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0 + 32,
                                 Cx_Neptune.getSubMatrix(6*j, 6*j + 5, 1, 1),
                                 Cy_Neptune.getSubMatrix(6*j, 6*j + 5, 1, 1),
                                 Cz_Neptune.getSubMatrix(6*j, 6*j + 5, 1, 1)) * 1e3;

    // Pluto calculation
    temp_start = 423 - 1;
    Matrix Cx_Pluto = PCtemp.getSubMatrix(temp_start, temp_start + 5, 1, 1);
    Matrix Cy_Pluto = PCtemp.getSubMatrix(temp_start + 6, temp_start + 11, 1, 1);
    Matrix Cz_Pluto = PCtemp.getSubMatrix(temp_start + 12, temp_start + 17, 1, 1);

    j = 0;
    Mjd0 = t1;

    positions.r_Pluto = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0 + 32,
                               Cx_Pluto.getSubMatrix(6*j, 6*j + 5, 1, 1),
                               Cy_Pluto.getSubMatrix(6*j, 6*j + 5, 1, 1),
                               Cz_Pluto.getSubMatrix(6*j, 6*j + 5, 1, 1)) * 1e3;

    // Nutations and Librations (if needed)
    // ... (similar implementation as above)

    // Final adjustments
    const double EMRAT = 81.30056907419062; // DE430
    const double EMRAT1 = 1.0 / (1.0 + EMRAT);

    positions.r_Earth = positions.r_Earth - positions.r_Moon * EMRAT1;
    positions.r_Mercury = positions.r_Mercury - positions.r_Earth;
    positions.r_Venus = positions.r_Venus - positions.r_Earth;
    positions.r_Mars = positions.r_Mars - positions.r_Earth;
    positions.r_Jupiter = positions.r_Jupiter - positions.r_Earth;
    positions.r_Saturn = positions.r_Saturn - positions.r_Earth;
    positions.r_Uranus = positions.r_Uranus - positions.r_Earth;
    positions.r_Neptune = positions.r_Neptune - positions.r_Earth;
    positions.r_Pluto = positions.r_Pluto - positions.r_Earth;
    positions.r_Sun = positions.r_Sun - positions.r_Earth;

    return positions;
}