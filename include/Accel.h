//
// Created by adria on 12/05/2025.
//

#ifndef PROYECTOTALLER_ACCEL_H
#define PROYECTOTALLER_ACCEL_H

#include "Matrix.h"


/**
 * @brief Calcula la aceleración total de un satélite en órbita terrestre.
 *
 * @param x  Tiempo desde epoch en segundos (offset para Mjd_UTC)
 * @param Y  Vector de estado [r; v] en ICRF/EME2000 (6×1)
 * @return   Derivada [v; a]  (6×1)
 */
Matrix Accel(double x, const Matrix& Y);

#endif //PROYECTOTALLER_ACCEL_H
