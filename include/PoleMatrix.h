//
// Created by adria on 11/05/2025.
//

#ifndef PROYECTOTALLER_POLEMATRIX_H
#define PROYECTOTALLER_POLEMATRIX_H


#include "Matrix.h"
#include "R_x.h"
#include "R_y.h"

/**
 * @brief Transformation from pseudo Earth-fixed to Earth-fixed coordinates
 *
 * @param xp Pole coordinate x (radians)
 * @param yp Pole coordinate y (radians)
 * @return Matrix Pole matrix
 *
 * Last modified: 2015/08/12 M. Mahooti
 */
Matrix PoleMatrix(double xp, double yp);

#endif //PROYECTOTALLER_POLEMATRIX_H
