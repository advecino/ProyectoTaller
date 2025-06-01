
#include "../include/JPL_Eph_DE430.h"
#include "../include/Cheb3D.h"
#include "../include/Sat_const.h"
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <iostream>

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
%--------------------------------------------------------------------------
*/

#include "../include/Matrix.h"
#include "../include/global.h"


static std::vector<double> getCoeffs(const std::vector<double>& pc,
                                     int startIndex,
                                     int blockCount,
                                     int step,
                                     int blockLen)
{
    int maxIdx = (startIndex - 1) + (blockCount - 1) * step + (blockLen - 1);
    if (maxIdx >= (int)pc.size()) {
        std::cerr << "Error en getCoeffs: intento de acceso fuera de rango.\n";
        std::cerr << "startIndex=" << startIndex << ", blockCount=" << blockCount
                  << ", step=" << step << ", blockLen=" << blockLen << "\n";
        std::cerr << "maxIdx=" << maxIdx << ", pc.size()=" << pc.size() << "\n";
        throw std::runtime_error("getCoeffs: acceso fuera de rango");
    }

    bool hayDatos = false;
    for (int i = (startIndex-1); i <= maxIdx; ++i) {
        if (pc[i] != 0.0) {
            hayDatos = true;
            break;
        }
    }
    if (!hayDatos) {
        std::cerr << "Advertencia: getCoeffs encontró solo ceros en el rango ("
                  << startIndex << " - " << maxIdx+1 << ")\n";
    }

    std::vector<double> out;
    out.reserve(blockCount * blockLen);
    for (int b = 0; b < blockCount; ++b) {
        int base = (startIndex - 1) + b * step;
        for (int k = 0; k < blockLen; ++k) {
            out.push_back(pc[base + k]);
        }
    }
    return out;
}

PlanetaryPositions JPL_Eph_DE430(double Mjd_TDB) {

    double JD = Mjd_TDB + 2400000.5;
    int rows = PC.getFilas(), cols = PC.getColumnas();
    int idx = 1;
    for (; idx <= rows; ++idx) {
        if (PC(idx,1) <= JD && JD <= PC(idx,2)) break;
    }
    if (idx > rows)
        throw std::invalid_argument("JPL_Eph_DE430: Mjd_TDB fuera de rango en PC.");

    std::vector<double> PCrow(cols);
    for (int c = 1; c <= cols; ++c)
        PCrow[c-1] = PC(idx, c);

    double t1 = PCrow[0] - 2400000.5;
    double dt = Mjd_TDB - t1;

    PlanetaryPositions eph;

    // Earth
    {
        auto Cx = getCoeffs(PCrow, 231, 2, 39, 13);
        auto Cy = getCoeffs(PCrow, 244, 2, 39, 13);
        auto Cz = getCoeffs(PCrow, 257, 2, 39, 13);

        int blocks = 2; double interval = 16.0;
        int j = std::min(int(dt/interval), blocks-1);
        double Mjd0 = t1 + j*interval;

        Matrix Cx_mat(13,1);
        for(int k=0;k<13;++k) Cx_mat(k+1,1) = Cx[j*13+k];

        Matrix Cy_mat(13,1);
        for(int k=0;k<13;++k) Cy_mat(k+1,1) = Cy[j*13+k];

        Matrix Cz_mat(13,1);
        for(int k=0;k<13;++k) Cz_mat(k+1,1) = Cz[j*13+k];

        Matrix v = Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0+interval, Cx_mat, Cy_mat, Cz_mat);
        eph.r_Earth = v * 1e3;
    }

    // Moon
    {
        auto Cx = getCoeffs(PCrow, 441, 8, 39, 13);
        auto Cy = getCoeffs(PCrow, 454, 8, 39, 13);
        auto Cz = getCoeffs(PCrow, 467, 8, 39, 13);

        int blocks = 8;
        double interval = 4.0;
        int j = std::min(int(dt/interval), blocks-1);
        double Mjd0 = t1 + j * interval;

        // Crear Matrix desde std::vector
        Matrix Cx_mat(13,1);
        for(int k=0;k<13;++k) Cx_mat(k+1,1) = Cx[j*13+k];

        Matrix Cy_mat(13,1);
        for(int k=0;k<13;++k) Cy_mat(k+1,1) = Cy[j*13+k];

        Matrix Cz_mat(13,1);
        for(int k=0;k<13;++k) Cz_mat(k+1,1) = Cz[j*13+k];

        Matrix v = Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0 + interval, Cx_mat, Cy_mat, Cz_mat);
        eph.r_Moon = v * 1e3;
    }


    // Sun
    {
        auto Cx = getCoeffs(PCrow, 753, 2, 33, 11);
        auto Cy = getCoeffs(PCrow, 764, 2, 33, 11);
        auto Cz = getCoeffs(PCrow, 775, 2, 33, 11);

        int blocks = 2;
        double interval = 16.0;
        int j = std::min(int(dt/interval), blocks-1);
        double Mjd0 = t1 + j * interval;

        // Convertir los std::vector a Matrix
        Matrix Cx_mat(11,1);
        for(int k=0; k<11; ++k) Cx_mat(k+1,1) = Cx[j*11+k];

        Matrix Cy_mat(11,1);
        for(int k=0; k<11; ++k) Cy_mat(k+1,1) = Cy[j*11+k];

        Matrix Cz_mat(11,1);
        for(int k=0; k<11; ++k) Cz_mat(k+1,1) = Cz[j*11+k];

        Matrix v = Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0 + interval, Cx_mat, Cy_mat, Cz_mat);
        eph.r_Sun = v * 1e3;
    }


    // Mercury
    {
        auto Cx = getCoeffs(PCrow, 3, 4, 42, 14);
        auto Cy = getCoeffs(PCrow, 17, 4, 42, 14);
        auto Cz = getCoeffs(PCrow, 31, 4, 42, 14);

        int blocks = 4;
        double interval = 8.0;
        int j = std::min(int(dt/interval), blocks-1);
        double Mjd0 = t1 + j * interval;

        // Convertir los std::vector a Matrix
        Matrix Cx_mat(14,1);
        for(int k=0; k<14; ++k) Cx_mat(k+1,1) = Cx[j*14+k];

        Matrix Cy_mat(14,1);
        for(int k=0; k<14; ++k) Cy_mat(k+1,1) = Cy[j*14+k];

        Matrix Cz_mat(14,1);
        for(int k=0; k<14; ++k) Cz_mat(k+1,1) = Cz[j*14+k];

        Matrix v = Cheb3D(Mjd_TDB, 14, Mjd0, Mjd0 + interval, Cx_mat, Cy_mat, Cz_mat);
        eph.r_Mercury = v * 1e3;
    }


    // Venus
    {
        auto Cx = getCoeffs(PCrow, 171, 2, 30, 10);
        auto Cy = getCoeffs(PCrow, 181, 2, 30, 10);
        auto Cz = getCoeffs(PCrow, 191, 2, 30, 10);

        int blocks = 2;
        double interval = 16.0;
        int j = std::min(int(dt/interval), blocks-1);
        double Mjd0 = t1 + j * interval;

        // Convertir los std::vector a Matrix
        Matrix Cx_mat(10,1);
        for(int k=0; k<10; ++k) Cx_mat(k+1,1) = Cx[j*10+k];

        Matrix Cy_mat(10,1);
        for(int k=0; k<10; ++k) Cy_mat(k+1,1) = Cy[j*10+k];

        Matrix Cz_mat(10,1);
        for(int k=0; k<10; ++k) Cz_mat(k+1,1) = Cz[j*10+k];

        Matrix v = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0 + interval, Cx_mat, Cy_mat, Cz_mat);
        eph.r_Venus = v * 1e3;
    }


    // Mars
    {
        auto Cx = getCoeffs(PCrow, 309, 1, 0, 11);
        auto Cy = getCoeffs(PCrow, 320, 1, 0, 11);
        auto Cz = getCoeffs(PCrow, 331, 1, 0, 11);

        double interval = 32.0;

        // Convertir std::vector a Matrix
        Matrix Cx_mat(11,1);
        for(int k=0; k<11; ++k) Cx_mat(k+1,1) = Cx[k];

        Matrix Cy_mat(11,1);
        for(int k=0; k<11; ++k) Cy_mat(k+1,1) = Cy[k];

        Matrix Cz_mat(11,1);
        for(int k=0; k<11; ++k) Cz_mat(k+1,1) = Cz[k];

        Matrix v = Cheb3D(Mjd_TDB, 11, t1, t1 + interval, Cx_mat, Cy_mat, Cz_mat);
        eph.r_Mars = v * 1e3;
    }


    // Jupiter
    {
        auto Cx = getCoeffs(PCrow, 342, 1, 0, 8);
        auto Cy = getCoeffs(PCrow, 350, 1, 0, 8);
        auto Cz = getCoeffs(PCrow, 358, 1, 0, 8);

        double interval = 32.0;

        // Convertir std::vector a Matrix
        Matrix Cx_mat(8,1);
        for(int k=0; k<8; ++k) Cx_mat(k+1,1) = Cx[k];

        Matrix Cy_mat(8,1);
        for(int k=0; k<8; ++k) Cy_mat(k+1,1) = Cy[k];

        Matrix Cz_mat(8,1);
        for(int k=0; k<8; ++k) Cz_mat(k+1,1) = Cz[k];

        Matrix v = Cheb3D(Mjd_TDB, 8, t1, t1 + interval, Cx_mat, Cy_mat, Cz_mat);
        eph.r_Jupiter = v * 1e3;
    }


    // Saturn
    {
        auto Cx = getCoeffs(PCrow, 366, 1, 0, 7);
        auto Cy = getCoeffs(PCrow, 373, 1, 0, 7);
        auto Cz = getCoeffs(PCrow, 380, 1, 0, 7);
        double interval = 32.0;

        Matrix Cx_mat(7,1);
        for(int k=0; k<7; ++k) Cx_mat(k+1,1) = Cx[k];
        Matrix Cy_mat(7,1);
        for(int k=0; k<7; ++k) Cy_mat(k+1,1) = Cy[k];
        Matrix Cz_mat(7,1);
        for(int k=0; k<7; ++k) Cz_mat(k+1,1) = Cz[k];

        eph.r_Saturn = Cheb3D(Mjd_TDB, 7, t1, t1+interval, Cx_mat, Cy_mat, Cz_mat) * (1e3);
    }

    // Uranus
    {
        auto Cx = getCoeffs(PCrow, 387, 1, 0, 6);
        auto Cy = getCoeffs(PCrow, 393, 1, 0, 6);
        auto Cz = getCoeffs(PCrow, 399, 1, 0, 6);
        double interval = 32.0;

        Matrix Cx_mat(6,1);
        for(int k=0; k<6; ++k) Cx_mat(k+1,1) = Cx[k];
        Matrix Cy_mat(6,1);
        for(int k=0; k<6; ++k) Cy_mat(k+1,1) = Cy[k];
        Matrix Cz_mat(6,1);
        for(int k=0; k<6; ++k) Cz_mat(k+1,1) = Cz[k];

        eph.r_Uranus = Cheb3D(Mjd_TDB, 6, t1, t1+interval, Cx_mat, Cy_mat, Cz_mat) * (1e3);
    }

    // Neptune
    {
        auto Cx = getCoeffs(PCrow, 405, 1, 0, 6);
        auto Cy = getCoeffs(PCrow, 411, 1, 0, 6);
        auto Cz = getCoeffs(PCrow, 417, 1, 0, 6);
        double interval = 32.0;

        Matrix Cx_mat(6,1);
        for(int k=0; k<6; ++k) Cx_mat(k+1,1) = Cx[k];
        Matrix Cy_mat(6,1);
        for(int k=0; k<6; ++k) Cy_mat(k+1,1) = Cy[k];
        Matrix Cz_mat(6,1);
        for(int k=0; k<6; ++k) Cz_mat(k+1,1) = Cz[k];

        eph.r_Neptune = Cheb3D(Mjd_TDB, 6, t1, t1+interval, Cx_mat, Cy_mat, Cz_mat) * (1e3);
    }

    // Pluto
    {
        auto Cx = getCoeffs(PCrow, 423, 1, 0, 6);
        auto Cy = getCoeffs(PCrow, 429, 1, 0, 6);
        auto Cz = getCoeffs(PCrow, 435, 1, 0, 6);
        double interval = 32.0;

        Matrix Cx_mat(6,1);
        for(int k=0; k<6; ++k) Cx_mat(k+1,1) = Cx[k];
        Matrix Cy_mat(6,1);
        for(int k=0; k<6; ++k) Cy_mat(k+1,1) = Cy[k];
        Matrix Cz_mat(6,1);
        for(int k=0; k<6; ++k) Cz_mat(k+1,1) = Cz[k];

        eph.r_Pluto = Cheb3D(Mjd_TDB, 6, t1, t1+interval, Cx_mat, Cy_mat, Cz_mat) * (1e3);
    }


    double EMRAT  = 81.30056907419062;
    double EMRAT1 = 1.0 / (1.0 + EMRAT);

    eph.r_Earth   = eph.r_Earth   - eph.r_Moon * (EMRAT1);
    eph.r_Mercury = eph.r_Mercury - eph.r_Earth;
    eph.r_Venus   = eph.r_Venus   - eph.r_Earth;
    eph.r_Mars    = eph.r_Mars    - eph.r_Earth;
    eph.r_Jupiter = eph.r_Jupiter - eph.r_Earth;
    eph.r_Saturn  = eph.r_Saturn  - eph.r_Earth;
    eph.r_Uranus  = eph.r_Uranus  - eph.r_Earth;
    eph.r_Neptune = eph.r_Neptune - eph.r_Earth;
    eph.r_Pluto   = eph.r_Pluto   - eph.r_Earth;
    eph.r_Sun     = eph.r_Sun     - eph.r_Earth;

    return eph;
}