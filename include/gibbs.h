#ifndef PROYECTOTALLER_GIBBS_H
#define PROYECTOTALLER_GIBBS_H

#include <string>

#include "Matrix.h"
#include "unit.h"

struct GibbsResult {
    Matrix v2;       ///< Vector de velocidad en r2 [m/s]
    double theta;    ///< Ángulo entre r1 y r2 [rad]
    double theta1;   ///< Ángulo entre r2 y r3 [rad]
    double copa;     ///< Ángulo de coplanaridad [rad]
    std::string error; ///< Indicador de éxito o error
    GibbsResult() : v2(3, 1) {}
};

/**
 * @brief Determina la velocidad en r2 usando el método de Gibbs.
 * @param r1 Vector de posición 1 (3x1).
 * @param r2 Vector de posición 2 (3x1).
 * @param r3 Vector de posición 3 (3x1).
 * @return Resultados del método de Gibbs.
 */
GibbsResult gibbs(const Matrix& r1, const Matrix& r2, const Matrix& r3);

#endif //PROYECTOTALLER_GIBBS_H
