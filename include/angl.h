//
// Created by adria on 20/04/2025.
//

#ifndef PROYECTOTALLER_ANGL_H
#define PROYECTOTALLER_ANGL_H


#include "Matrix.h"

/**
 * @brief Calcula el ángulo entre dos vectores (entre -π y π)
 * @param vec1 Primer vector (Matrix 3x1)
 * @param vec2 Segundo vector (Matrix 3x1)
 * @return Ángulo en radianes entre los dos vectores.
 *         Devuelve 999999.1 si no se puede calcular.
 */
double angl(const Matrix& vec1, const Matrix& vec2);

#endif //PROYECTOTALLER_ANGL_H
