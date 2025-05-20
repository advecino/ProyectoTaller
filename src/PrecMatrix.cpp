#include "../include/PrecMatrix.h"
#include <cmath>


Matrix PrecMatrix(double Mjd_1, double Mjd_2) {
    const double MJD_J2000 = 51544.5;
    const double Arcs = 3600.0;

    double T = (Mjd_1 - MJD_J2000) / 36525.0;
    double dT = (Mjd_2 - Mjd_1) / 36525.0;

    // Precession angles
    double zeta = ((2306.2181 + (1.39656 - 0.000139 * T) * T) +
                   ((0.30188 - 0.000344 * T) + 0.017998 * dT) * dT) * dT / Arcs;

    double z = zeta + ((0.79280 + 0.000411 * T) + 0.000205 * dT) * dT * dT / Arcs;

    double theta = ((2004.3109 - (0.85330 + 0.000217 * T) * T) -
                    ((0.42665 + 0.000217 * T) + 0.041833 * dT) * dT) * dT / Arcs;

    zeta *= M_PI / 180.0;
    z *= M_PI / 180.0;
    theta *= M_PI / 180.0;

    // Precession matrix
    return R_z(-z) * R_y(theta) * R_z(-zeta);
}
