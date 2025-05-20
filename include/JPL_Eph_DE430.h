#ifndef PROYECTOTALLER_JPL_EPH_DE430_H
#define PROYECTOTALLER_JPL_EPH_DE430_H


#include "Matrix.h"

/**
 * @brief Calcula posiciones planetarias con efemérides DE430.
 */
struct PlanetaryPositions {
    Matrix r_Mercury, r_Venus, r_Earth, r_Mars, r_Jupiter,
            r_Saturn, r_Uranus, r_Neptune, r_Pluto,
            r_Moon, r_Sun;
    Matrix Nutations, Librations;

    PlanetaryPositions()
            : r_Mercury(3,1), r_Venus(3,1), r_Earth(3,1), r_Mars(3,1),
              r_Jupiter(3,1), r_Saturn(3,1), r_Uranus(3,1), r_Neptune(3,1),
              r_Pluto(3,1), r_Moon(3,1), r_Sun(3,1),
              Nutations(10*4,1), Librations(10*4,1) {}
};


/**
 * @brief Posiciones ecuatoriales de Sol, Luna y planetas usando efemérides JPL DE430.
 * @param Mjd_TDB Fecha juliana modificada (Tiempo Baricéntrico Dinámico).
 * @return Estructura con posiciones (ICRF) y coeficientes de nutaciones/libraciones.
 */
PlanetaryPositions JPL_Eph_DE430(double Mjd_TDB);


#endif //PROYECTOTALLER_JPL_EPH_DE430_H
