#ifndef PROYECTOTALLER_ANGL_H
#define PROYECTOTALLER_ANGL_H


#include "Matrix.h"

/**
 * @brief Ángulo entre dos vectores (de -pi a pi).
 * @param vec1 Primer vector (3x1).
 * @param vec2 Segundo vector (3x1).
 * @return Ángulo en radianes. Devuelve 999999.1 si no se puede calcular.
 */
double angl(const Matrix& vec1, const Matrix& vec2);

#endif //PROYECTOTALLER_ANGL_H
