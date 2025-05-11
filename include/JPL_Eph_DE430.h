#ifndef PROYECTOTALLER_JPL_EPH_DE430_H
#define PROYECTOTALLER_JPL_EPH_DE430_H


#include "Matrix.h"
#include "global.h"
#include "Cheb3D.h"


/**
 * @brief Structure containing planetary positions from JPL Ephemeris DE430
 */
struct PlanetaryPositions {
    Matrix r_Mercury;  // [m] ICRF position vector
    Matrix r_Venus;    // [m] ICRF position vector
    Matrix r_Earth;    // [m] ICRF position vector (from SSB)
    Matrix r_Mars;     // [m] ICRF position vector
    Matrix r_Jupiter;  // [m] ICRF position vector
    Matrix r_Saturn;   // [m] ICRF position vector
    Matrix r_Uranus;   // [m] ICRF position vector
    Matrix r_Neptune;  // [m] ICRF position vector
    Matrix r_Pluto;    // [m] ICRF position vector
    Matrix r_Moon;     // [m] ICRF position vector (geocentric)
    Matrix r_Sun;      // [m] ICRF position vector (geocentric)

    PlanetaryPositions()
            : r_Mercury(3, 1), r_Venus(3, 1), r_Earth(3, 1),
              r_Mars(3, 1), r_Jupiter(3, 1), r_Saturn(3, 1),
              r_Uranus(3, 1), r_Neptune(3, 1), r_Pluto(3, 1),
              r_Moon(3, 1), r_Sun(3, 1) {}
};

/**
 * @brief Computes positions of sun, moon, and planets using JPL Ephemeris DE430
 *
 * @param Mjd_TDB Modified Julian Date of TDB (Barycentric Dynamical Time)
 * @return PlanetaryPositions Struct containing all position vectors in ICRF frame [m]
 *
 * @note Light-time is already taken into account in the ephemeris
 * @last_modified 2018/01/11 M. Mahooti
 */
PlanetaryPositions JPL_Eph_DE430(double Mjd_TDB, const Matrix& PC);


#endif //PROYECTOTALLER_JPL_EPH_DE430_H
