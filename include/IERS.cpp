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



IERSResult IERS(const Matrix& eop, double Mjd_UTC, std::string& interp) {
    IERSResult result;
    const double Arcs = 3600.0 * 180.0 / M_PI; // Define Arcs constant

    if (interp.empty()) {
        interp = "n";
    }

    double mjd = floor(Mjd_UTC);
    int i = -1;
    int cols = eop.getColumnas();

    // Find the column index where MJD matches
    for (int col = 0; col < cols; ++col) {
        if (std::abs(eop(4, col) - mjd < 1e-8)) { // Note: using row 4 (0-based) for MJD
            i = col;
            break;
        }
    }

    if (i == -1 || i >= cols-1) {
        throw std::runtime_error("Could not find MJD in EOP data or insufficient data for interpolation");
    }

    if (interp == "l") {
        // Calculate interpolation fraction
        double mfme = 1440.0 * (Mjd_UTC - mjd);
        double fixf = mfme / 1440.0;

        // Linear interpolation
        result.x_pole   = eop(5, i)   + (eop(5, i+1)   - eop(5, i))   * fixf;
        result.y_pole   = eop(6, i)   + (eop(6, i+1)   - eop(6, i))   * fixf;
        result.UT1_UTC  = eop(7, i)   + (eop(7, i+1)   - eop(7, i))   * fixf;
        result.LOD      = eop(8, i)   + (eop(8, i+1)   - eop(8, i))   * fixf;
        result.dpsi     = eop(9, i)   + (eop(9, i+1)   - eop(9, i))   * fixf;
        result.deps     = eop(10, i)  + (eop(10, i+1)  - eop(10, i))  * fixf;
        result.dx_pole  = eop(11, i)  + (eop(11, i+1)  - eop(11, i))  * fixf;
        result.dy_pole  = eop(12, i)  + (eop(12, i+1)  - eop(12, i))  * fixf;
        result.TAI_UTC  = eop(13, i); // No interpolation for TAI_UTC

        // Convert from arcseconds to radians
        result.x_pole  /= Arcs;
        result.y_pole  /= Arcs;
        result.dpsi    /= Arcs;
        result.deps    /= Arcs;
        result.dx_pole /= Arcs;
        result.dy_pole /= Arcs;
    }
    else if (interp == "n") {
        // No interpolation - direct values
        result.x_pole   = eop(5, i) / Arcs;
        result.y_pole   = eop(6, i) / Arcs;
        result.UT1_UTC  = eop(7, i);
        result.LOD      = eop(8, i);
        result.dpsi     = eop(9, i) / Arcs;
        result.deps     = eop(10, i) / Arcs;
        result.dx_pole  = eop(11, i) / Arcs;
        result.dy_pole  = eop(12, i) / Arcs;
        result.TAI_UTC  = eop(13, i);
    }

    return result;
}
