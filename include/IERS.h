//
// Created by adria on 11/05/2025.
//

#ifndef PROYECTOTALLER_IERS_H
#define PROYECTOTALLER_IERS_H



#include "Matrix.h"
#include <cmath>
#include "Sat_const.h"


struct IERSResult {
    double x_pole;    // Pole coordinate [rad]
    double y_pole;    // Pole coordinate [rad]
    double UT1_UTC;   // UT1-UTC time difference [s]
    double LOD;       // Length of day [s]
    double dpsi;      // Nutation angle [rad]
    double deps;      // Nutation angle [rad]
    double dx_pole;   // Pole coordinate [rad]
    double dy_pole;   // Pole coordinate [rad]
    double TAI_UTC;   // TAI-UTC time difference [s]

    // Constructor por defecto

};

// Declaración de la función IERS
IERSResult IERS(const Matrix& eop, double Mjd_UTC, std::string& interp);

#endif //PROYECTOTALLER_IERS_H
