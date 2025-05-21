#include "../include/JPL_Eph_DE430.h"
#include "../include/Cheb3D.h"
#include <cmath>
#include <stdexcept>
#include <vector>
#include <iostream>
#include "../include/global.h"

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


PlanetaryPositions JPL_Eph_DE430(double Mjd_TDB) {

    PlanetaryPositions positions;
    double JD = Mjd_TDB + 2400000.5;

    int i = -1;
    for (int r = 1; r <= PC.getFilas(); ++r) {
        if (PC(r,1) <= JD && JD <= PC(r,2)) {
            i = r;
            break;
        }
    }
    if (i < 1) {
        throw std::runtime_error("No se encontró intervalo válido para JD en PC");
    }

    Matrix PCtemp = PC.getFila(i);

    double t1 = PCtemp(1,1) - 2400000.5; // MJD al inicio del intervalo
    double dt = Mjd_TDB - t1;

    // ========== TIERRA ==========
    std::vector<int> temp = {231, 244, 257, 270}; // 231:13:270
    Matrix Cx_Earth = PCtemp.getSubMatrix(1, 1, temp[0], temp[1] - 1);
    Matrix Cy_Earth = PCtemp.getSubMatrix(1, 1, temp[1], temp[2] - 1);
    Matrix Cz_Earth = PCtemp.getSubMatrix(1, 1, temp[2], temp[3] - 1);

    // Añadir segmentos adicionales
    for (int seg = 1; seg <= 1; ++seg)  {
        for (int j = 0; j < 4; ++j) {
            temp[j] += 39;
        }
        Matrix Cx = PCtemp.getSubMatrix(1, 1, temp[0], temp[1] - 1);
        Matrix Cy = PCtemp.getSubMatrix(1, 1, temp[1], temp[2] - 1);
        Matrix Cz = PCtemp.getSubMatrix(1, 1, temp[2], temp[3] - 1);


        Cx_Earth = Matrix::concatenar(Cx_Earth, Cx);
        Cy_Earth = Matrix::concatenar(Cy_Earth, Cy);
        Cz_Earth = Matrix::concatenar(Cz_Earth, Cz);
    }

    int j;
    double Mjd0;
    if (dt >= 0 && dt <= 16) {
        j = 0;
        Mjd0 = t1;
    } else if (dt > 16 && dt <= 32) {
        j = 1;
        Mjd0 = t1 + 16;
    } else {
        throw std::runtime_error("dt fuera del rango esperado [0,32]");
    }

    int c1 = 13*j + 1;
    int c2 = 13*(j+1);
    std::cout << "[DBG] Cx_Earth columnas antes de Cheb3D: " << Cx_Earth.getColumnas() << std::endl;
    std::cout << "[DBG] c1 = " << c1 << ", c2 = " << c2 << std::endl;




    positions.r_Earth = Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0+16,
                               Cx_Earth.getSubMatrix(1, 1, c1, c2),
                               Cy_Earth.getSubMatrix(1, 1, c1, c2),
                               Cz_Earth.getSubMatrix(1, 1, c1, c2)) * 1e3;

    // ========== LUNA ==========
    temp = {441, 454, 467, 480};
    Matrix Cx_Moon = PCtemp.getSubMatrix(1, 1, temp[0], temp[1] - 1);
    Matrix Cy_Moon = PCtemp.getSubMatrix(1, 1, temp[1], temp[2] - 1);
    Matrix Cz_Moon = PCtemp.getSubMatrix(1, 1, temp[2], temp[3] - 1);


    for (int seg = 1; seg <= 7; ++seg) {
        for (int j = 0; j < 4; ++j) {
            temp[j] += 39;
        }
        Matrix Cx = PCtemp.getSubMatrix(1, 1, temp[0], temp[1] - 1);
        Matrix Cy = PCtemp.getSubMatrix(1, 1, temp[1], temp[2] - 1);
        Matrix Cz = PCtemp.getSubMatrix(1, 1, temp[2], temp[3] - 1);


        Cx_Moon = Matrix::concatenar(Cx_Moon, Cx);
        Cy_Moon = Matrix::concatenar(Cy_Moon, Cy);
        Cz_Moon = Matrix::concatenar(Cz_Moon, Cz);
    }

    if (dt <= 4) j = 0;
    else if (dt <= 8) j = 1;
    else if (dt <= 12) j = 2;
    else if (dt <= 16) j = 3;
    else if (dt <= 20) j = 4;
    else if (dt <= 24) j = 5;
    else if (dt <= 28) j = 6;
    else j = 7;

    Mjd0 = t1 + 4*j;
    positions.r_Moon = Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0+4,
                              Cx_Moon.getSubMatrix(1, 1, 13*j+1, 13*j+13),
                              Cy_Moon.getSubMatrix(1, 1, 13*j+1, 13*j+13),
                              Cz_Moon.getSubMatrix(1, 1, 13*j+1, 13*j+13)) * 1e3;


    // ========== SOL ==========
    temp = {753, 764, 775, 786};
    Matrix Cx_Sun = PCtemp.getSubMatrix(1, 1, temp[0], temp[1] - 1);
    Matrix Cy_Sun = PCtemp.getSubMatrix(1, 1, temp[1], temp[2] - 1);
    Matrix Cz_Sun = PCtemp.getSubMatrix(1, 1, temp[2], temp[3] - 1);


    for (int j = 0; j < 4; ++j) {
        temp[j] += 33;
    }
    Matrix Cx = PCtemp.getSubMatrix(1, 1, temp[0], temp[1] - 1);
    Matrix Cy = PCtemp.getSubMatrix(1, 1, temp[1], temp[2] - 1);
    Matrix Cz = PCtemp.getSubMatrix(1, 1, temp[2], temp[3] - 1);


    Cx_Sun = Matrix::concatenar(Cx_Sun, Cx);
    Cy_Sun = Matrix::concatenar(Cy_Sun, Cy);
    Cz_Sun = Matrix::concatenar(Cz_Sun, Cz);

    j = (dt <= 16) ? 0 : 1;
    Mjd0 = t1 + 16*j;
    positions.r_Sun = Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0 + 16,
                             Cx_Sun.getSubMatrix(1, 1, 11 * j + 1, 11 * j + 11),
                             Cy_Sun.getSubMatrix(1, 1, 11 * j + 1, 11 * j + 11),
                             Cz_Sun.getSubMatrix(1, 1, 11 * j + 1, 11 * j + 11)) * 1e3;


    // ========== MERCURIO ==========
    temp = {3, 17, 31, 45};
    Matrix Cx_Mercury = PCtemp.getSubMatrix(1, 1, temp[0], temp[1] - 1);
    Matrix Cy_Mercury = PCtemp.getSubMatrix(1, 1, temp[1], temp[2] - 1);
    Matrix Cz_Mercury = PCtemp.getSubMatrix(1, 1, temp[2], temp[3] - 1);


    for (int seg = 1; seg <= 3; ++seg) {
        for (int j = 0; j < 4; ++j) {
            temp[j] += 42;
        }
        Matrix Cx = PCtemp.getSubMatrix(1, 1, temp[0], temp[1] - 1);
        Matrix Cy = PCtemp.getSubMatrix(1, 1, temp[1], temp[2] - 1);
        Matrix Cz = PCtemp.getSubMatrix(1, 1, temp[2], temp[3] - 1);


        Cx_Mercury = Matrix::concatenar(Cx_Mercury, Cx);
        Cy_Mercury = Matrix::concatenar(Cy_Mercury, Cy);
        Cz_Mercury = Matrix::concatenar(Cz_Mercury, Cz);
    }

    if (dt <= 8) j = 0;
    else if (dt <= 16) j = 1;
    else if (dt <= 24) j = 2;
    else j = 3;

    Mjd0 = t1 + 8*j;
    positions.r_Mercury = Cheb3D(Mjd_TDB, 14, Mjd0, Mjd0 + 8,
                                 Cx_Mercury.getSubMatrix(1, 1, 14 * j + 1, 14 * j + 14),
                                 Cy_Mercury.getSubMatrix(1, 1, 14 * j + 1, 14 * j + 14),
                                 Cz_Mercury.getSubMatrix(1, 1, 14 * j + 1, 14 * j + 14)) * 1e3;

    // ========== VENUS ==========
    temp = {171, 181, 191, 201};
    Matrix Cx_Venus = PCtemp.getSubMatrix(1, 1, temp[0], temp[1] - 1);
    Matrix Cy_Venus = PCtemp.getSubMatrix(1, 1, temp[1], temp[2] - 1);
    Matrix Cz_Venus = PCtemp.getSubMatrix(1, 1, temp[2], temp[3] - 1);

    for (int j = 0; j < 4; ++j) temp[j] += 30;
    Cx = PCtemp.getSubMatrix(1, 1, temp[0], temp[1] - 1);
    Cy = PCtemp.getSubMatrix(1, 1, temp[1], temp[2] - 1);
    Cz = PCtemp.getSubMatrix(1, 1, temp[2], temp[3] - 1);
    Cx_Venus = Matrix::concatenar(Cx_Venus, Cx);
    Cy_Venus = Matrix::concatenar(Cy_Venus, Cy);
    Cz_Venus = Matrix::concatenar(Cz_Venus, Cz);

    j = (dt <= 16) ? 0 : 1;
    Mjd0 = t1 + 16 * j;
    positions.r_Venus = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+16,
                               Cx_Venus.getSubMatrix(1, 1, 10*j+1, 10*j+10),
                               Cy_Venus.getSubMatrix(1, 1, 10*j+1, 10*j+10),
                               Cz_Venus.getSubMatrix(1, 1, 10*j+1, 10*j+10)) * 1e3;

// ========== MARTE ==========
    temp = {309, 320, 331, 342};
    Matrix Cx_Mars = PCtemp.getSubMatrix(1, 1, temp[0], temp[1] - 1);
    Matrix Cy_Mars = PCtemp.getSubMatrix(1, 1, temp[1], temp[2] - 1);
    Matrix Cz_Mars = PCtemp.getSubMatrix(1, 1, temp[2], temp[3] - 1);

    positions.r_Mars = Cheb3D(Mjd_TDB, 11, t1, t1+32,
                              Cx_Mars.getSubMatrix(1, 1, 1, 11),
                              Cy_Mars.getSubMatrix(1, 1, 1, 11),
                              Cz_Mars.getSubMatrix(1, 1, 1, 11)) * 1e3;

// ========== JUPITER ==========
    temp = {342, 350, 358, 366};
    Matrix Cx_Jupiter = PCtemp.getSubMatrix(1, 1, temp[0], temp[1] - 1);
    Matrix Cy_Jupiter = PCtemp.getSubMatrix(1, 1, temp[1], temp[2] - 1);
    Matrix Cz_Jupiter = PCtemp.getSubMatrix(1, 1, temp[2], temp[3] - 1);

    positions.r_Jupiter = Cheb3D(Mjd_TDB, 8, t1, t1+32,
                                 Cx_Jupiter.getSubMatrix(1, 1, 1, 8),
                                 Cy_Jupiter.getSubMatrix(1, 1, 1, 8),
                                 Cz_Jupiter.getSubMatrix(1, 1, 1, 8)) * 1e3;

// ========== SATURNO ==========
    temp = {366, 373, 380, 387};
    Matrix Cx_Saturn = PCtemp.getSubMatrix(1, 1, temp[0], temp[1] - 1);
    Matrix Cy_Saturn = PCtemp.getSubMatrix(1, 1, temp[1], temp[2] - 1);
    Matrix Cz_Saturn = PCtemp.getSubMatrix(1, 1, temp[2], temp[3] - 1);

    positions.r_Saturn = Cheb3D(Mjd_TDB, 7, t1, t1+32,
                                Cx_Saturn.getSubMatrix(1, 1, 1, 7),
                                Cy_Saturn.getSubMatrix(1, 1, 1, 7),
                                Cz_Saturn.getSubMatrix(1, 1, 1, 7)) * 1e3;

// ========== URANO ==========
    temp = {387, 393, 399, 405};
    Matrix Cx_Uranus = PCtemp.getSubMatrix(1, 1, temp[0], temp[1] - 1);
    Matrix Cy_Uranus = PCtemp.getSubMatrix(1, 1, temp[1], temp[2] - 1);
    Matrix Cz_Uranus = PCtemp.getSubMatrix(1, 1, temp[2], temp[3] - 1);

    positions.r_Uranus = Cheb3D(Mjd_TDB, 6, t1, t1+32,
                                Cx_Uranus.getSubMatrix(1, 1, 1, 6),
                                Cy_Uranus.getSubMatrix(1, 1, 1, 6),
                                Cz_Uranus.getSubMatrix(1, 1, 1, 6)) * 1e3;

// ========== NEPTUNO ==========
    temp = {405, 411, 417, 423};
    Matrix Cx_Neptune = PCtemp.getSubMatrix(1, 1, temp[0], temp[1] - 1);
    Matrix Cy_Neptune = PCtemp.getSubMatrix(1, 1, temp[1], temp[2] - 1);
    Matrix Cz_Neptune = PCtemp.getSubMatrix(1, 1, temp[2], temp[3] - 1);

    positions.r_Neptune = Cheb3D(Mjd_TDB, 6, t1, t1+32,
                                 Cx_Neptune.getSubMatrix(1, 1, 1, 6),
                                 Cy_Neptune.getSubMatrix(1, 1, 1, 6),
                                 Cz_Neptune.getSubMatrix(1, 1, 1, 6)) * 1e3;

// ========== PLUTÓN ==========
    temp = {423, 429, 435, 441};
    Matrix Cx_Pluto = PCtemp.getSubMatrix(1, 1, temp[0], temp[1] - 1);
    Matrix Cy_Pluto = PCtemp.getSubMatrix(1, 1, temp[1], temp[2] - 1);
    Matrix Cz_Pluto = PCtemp.getSubMatrix(1, 1, temp[2], temp[3] - 1);

    positions.r_Pluto = Cheb3D(Mjd_TDB, 6, t1, t1+32,
                               Cx_Pluto.getSubMatrix(1, 1, 1, 6),
                               Cy_Pluto.getSubMatrix(1, 1, 1, 6),
                               Cz_Pluto.getSubMatrix(1, 1, 1, 6)) * 1e3;


    // ========== AJUSTES FINALES ==========
    const double EMRAT = 81.30056907419062; // DE430
    const double EMRAT1 = 1.0/(1.0 + EMRAT);

    positions.r_Earth = positions.r_Earth - positions.r_Moon * EMRAT1;
    positions.r_Mercury = -1*positions.r_Earth + positions.r_Mercury;
    positions.r_Venus = -1*positions.r_Earth + positions.r_Venus;
    positions.r_Mars = -1*positions.r_Earth + positions.r_Mars;
    positions.r_Jupiter = -1*positions.r_Earth + positions.r_Jupiter;
    positions.r_Saturn = -1*positions.r_Earth + positions.r_Saturn;
    positions.r_Uranus = -1*positions.r_Earth + positions.r_Uranus;
    positions.r_Neptune = -1*positions.r_Earth + positions.r_Neptune;
    positions.r_Pluto = -1*positions.r_Earth + positions.r_Pluto;
    positions.r_Sun = -1*positions.r_Earth + positions.r_Sun;

    return positions;
}