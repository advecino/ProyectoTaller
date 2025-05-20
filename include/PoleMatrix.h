#ifndef PROYECTOTALLER_POLEMATRIX_H
#define PROYECTOTALLER_POLEMATRIX_H


#include "Matrix.h"
#include "R_x.h"
#include "R_y.h"

/**
 * @brief Transformación de coordenadas pseudo-fijas a Tierra-fijas.
 * @param xp Coordenada del polo x [rad].
 * @param yp Coordenada del polo y [rad].
 * @return Matriz de corrección por movimiento polar.
 */
Matrix PoleMatrix(double xp, double yp);

#endif //PROYECTOTALLER_POLEMATRIX_H
