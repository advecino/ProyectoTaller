//
// Created by adria on 10/05/2025.
//

#ifndef PROYECTOTALLER_UNIT_H
#define PROYECTOTALLER_UNIT_H

#include <cmath>
#include "Matrix.h"

/**
 * @brief Devuelve el vector unitario correspondiente a un vector dado.
 * @param vec Vector de entrada (3x1).
 * @return Vector unitario (3x1) o vector cero si la magnitud es muy pequeña.
 */
Matrix unit(const Matrix& vec);


#endif //PROYECTOTALLER_UNIT_H
