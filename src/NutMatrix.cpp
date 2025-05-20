#include "../include/NutMatrix.h"

Matrix NutMatrix(double Mjd_TT) {
    // Mean obliquity of the ecliptic
    double eps = MeanObliquity(Mjd_TT);

    // Nutation in longitude and obliquity
    double dpsi, deps;  // Cambiado de Matrix a double
    NutAngles(Mjd_TT, dpsi, deps);

    // Transformation from mean to true equator and equinox
    return R_x(-eps-deps) * R_z(-dpsi) * R_x(+eps);
}