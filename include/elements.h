#ifndef PROYECTOTALLER_ELEMENTS_H
#define PROYECTOTALLER_ELEMENTS_H


#include <cmath>
#include "Matrix.h"
#include "SAT_Const.h"

struct KeplerianElements {
    double p;      ///< Semilatus rectum [m]
    double a;      ///< Semieje mayor [m]
    double e;      ///< Excentricidad
    double i;      ///< Inclinación [rad]
    double Omega;  ///< Longitud del nodo ascendente [rad]
    double omega;  ///< Argumento del pericentro [rad]
    double M;      ///< Anomalía media [rad]
};

/**
 * @brief Calcula los elementos keplerianos osculantes para órbitas elípticas.
 * @throws std::invalid_argument si la órbita es circular o ecuatorial.
 */
KeplerianElements elements(const Matrix& r, const Matrix& v);

#endif //PROYECTOTALLER_ELEMENTS_H
