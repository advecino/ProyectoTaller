//
// Created by adria on 11/05/2025.
//asta squi

#include <stdexcept>
#include "IERS.h"

/*
%--------------------------------------------------------------------------
%
% IERS: Management of IERS time and polar motion data
%
% Last modified:   2018/02/01   M. Mahooti
%
%--------------------------------------------------------------------------*/



IERSResult IERS(Matrix eop, double Mjd_UTC, std::string interp) {
    IERSResult result;

    int filas = eop.getFilas();
    int cols = eop.getColumnas();

    if (cols < 2) {
        throw std::runtime_error("EOP data is empty or incomplete");
    }

    // Buscar los dos índices i e i+1 tales que eop(4, i) <= Mjd_UTC < eop(4, i+1)
    int i = -1;
    for (int k = 1; k < cols; ++k) {
        double mjd_k = eop(4, k);
        double mjd_k1 = eop(4, k + 1);

        if (mjd_k <= Mjd_UTC && Mjd_UTC < mjd_k1) {
            i = k;
            break;
        }
    }

    if (i == -1 || i + 1 > cols) {
        throw std::runtime_error("Could not find MJD interval in EOP data or insufficient data for interpolation");
    }

    if (interp == "l") {
        double mjd0 = eop(4, i);
        double mjd1 = eop(4, i + 1);
        double fixf = (Mjd_UTC - mjd0) / (mjd1 - mjd0);

        result.x_pole   = eop(5, i)   + (eop(5, i + 1)   - eop(5, i))   * fixf;
        result.y_pole   = eop(6, i)   + (eop(6, i + 1)   - eop(6, i))   * fixf;
        result.UT1_UTC  = eop(7, i)   + (eop(7, i + 1)   - eop(7, i))   * fixf;
        result.LOD      = eop(8, i)   + (eop(8, i + 1)   - eop(8, i))   * fixf;
        result.dpsi     = eop(9, i)   + (eop(9, i + 1)   - eop(9, i))   * fixf;
        result.deps     = eop(10, i)  + (eop(10, i + 1)  - eop(10, i))  * fixf;
        result.dx_pole  = eop(11, i)  + (eop(11, i + 1)  - eop(11, i))  * fixf;
        result.dy_pole  = eop(12, i)  + (eop(12, i + 1)  - eop(12, i))  * fixf;
        result.TAI_UTC  = eop(13, i); // No interpolación para TAI_UTC

    } else {
        // Sin interpolación: tomar el valor más cercano por debajo
        result.x_pole   = eop(5, i);
        result.y_pole   = eop(6, i);
        result.UT1_UTC  = eop(7, i);
        result.LOD      = eop(8, i);
        result.dpsi     = eop(9, i);
        result.deps     = eop(10, i);
        result.dx_pole  = eop(11, i);
        result.dy_pole  = eop(12, i);
        result.TAI_UTC  = eop(13, i);
    }

    // Conversión de arcosegundos a radianes
    result.x_pole  /= Arcs;
    result.y_pole  /= Arcs;
    result.dpsi    /= Arcs;
    result.deps    /= Arcs;
    result.dx_pole /= Arcs;
    result.dy_pole /= Arcs;

    return result;
}
