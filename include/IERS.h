//
// Created by adria on 11/05/2025.
//

#ifndef PROYECTOTALLER_IERS_H
#define PROYECTOTALLER_IERS_H



#include "Matrix.h"
#include <cmath>
#include "Sat_const.h"


struct IERSResult {
    double x_pole;
    double y_pole;
    double UT1_UTC;
    double LOD;
    double dpsi;
    double deps;
    double dx_pole;
    double dy_pole;
    double TAI_UTC;
};

IERSResult IERS(const Matrix& eop, double Mjd_UTC, const std::string& interp = "n");

#endif //PROYECTOTALLER_IERS_H
