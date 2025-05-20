#ifndef PROYECTOTALLER_G_ACCELHARMONIC_H
#define PROYECTOTALLER_G_ACCELHARMONIC_H


#include "Matrix.h"

/**
 * @brief Gradiente del campo gravitatorio armónico de la Tierra.
 * @param r Vector de posición del satélite en el sistema true-of-date [m].
 * @param U Matriz de transformación al sistema ligado al cuerpo.
 * @param n_max Grado del modelo gravitatorio.
 * @param m_max Orden del modelo gravitatorio.
 * @return Matriz gradiente (G = da/dr) en el sistema true-of-date [1/s²].
 */
Matrix G_AccelHarmonic(Matrix& r, Matrix& U, int n_max, int m_max);


#endif //PROYECTOTALLER_G_ACCELHARMONIC_H
