//
// Created by adria on 11/05/2025.
//

#ifndef PROYECTOTALLER_IERS_H
#define PROYECTOTALLER_IERS_H



#include "Matrix.h"



/**
 * @brief Result container for IERS Earth orientation parameters
 */
struct IERSResult {
    double x_pole;   ///< Pole coordinate [rad]
    double y_pole;   ///< Pole coordinate [rad]
    double UT1_UTC;  ///< UT1-UTC time difference [s]
    double LOD;      ///< Length of day [s]
    double dpsi;     ///< Nutation in longitude [rad]
    double deps;     ///< Nutation in obliquity [rad]
    double dx_pole;  ///< Pole coordinate [rad]
    double dy_pole;  ///< Pole coordinate [rad]
    double TAI_UTC;  ///< TAI-UTC time difference [s]
};

/**
 * @brief Management of IERS time and polar motion data
 * @param eop     Matrix of Earth orientation parameters (at least 13 rows)
 * @param Mjd_UTC Modified Julian Date UTC
 * @param interp  'l' for linear interpolation, 'n' for no interpolation
 * @return Filled IERSResult struct
 */
IERSResult IERS(Matrix& eop, double Mjd_UTC, char interp='n');
#endif //PROYECTOTALLER_IERS_H
