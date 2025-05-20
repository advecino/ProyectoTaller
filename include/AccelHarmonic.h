#ifndef PROYECTOTALLER_ACCELHARMONIC_H
#define PROYECTOTALLER_ACCELHARMONIC_H

#include "Matrix.h"

/**
 * @brief Calcula la aceleración debida al campo gravitatorio armónico del cuerpo central.
 * @param r Vector de posición (3x1) del satélite en el sistema inercial [m].
 * @param E Matriz de transformación del sistema inercial al sistema ligado al cuerpo (body-fixed).
 * @param n_max Grado máximo del desarrollo armónico.
 * @param m_max Orden máximo del desarrollo armónico (debe ser menor o igual que n_max).
 * @return Vector aceleración (3x1) en el sistema inercial [m/s²].
 */
Matrix AccelHarmonic(Matrix& r, Matrix& E, int n_max, int m_max);

#endif //PROYECTOTALLER_ACCELHARMONIC_H
