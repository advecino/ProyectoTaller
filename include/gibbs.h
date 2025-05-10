//
// Created by adria on 10/05/2025.
//

#ifndef PROYECTOTALLER_GIBBS_H
#define PROYECTOTALLER_GIBBS_H

#include <string>

#include "Matrix.h"
#include "unit.h"

struct GibbsResult {
    Matrix v2;
    double theta;
    double theta1;
    double copa;
    std::string error;
    GibbsResult() : v2(3, 1) {}
};


/**
 * @brief Implementación del método de Gibbs para determinación orbital
 * @param r1 Vector de posición 1 (3x1)
 * @param r2 Vector de posición 2 (3x1)
 * @param r3 Vector de posición 3 (3x1)
 * @return Estructura con los resultados
 */
GibbsResult gibbs(const Matrix& r1, const Matrix& r2, const Matrix& r3);

#endif //PROYECTOTALLER_GIBBS_H
