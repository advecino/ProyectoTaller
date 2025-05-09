//
// Created by adria on 08/05/2025.
//

#ifndef PROYECTOTALLER_G_ACCELHARMONIC_H
#define PROYECTOTALLER_G_ACCELHARMONIC_H


#include "Matrix.h"

/**
 * @brief Computes the gradient of the Earth's harmonic gravity field
 *
 * @param r Satellite position vector in the true-of-date system [m]
 * @param U Transformation matrix to body-fixed system
 * @param n_max Gravity model degree
 * @param m_max Gravity model order
 * @return Matrix Gradient (G=da/dr) in the true-of-date system [1/s²]
 */
Matrix G_AccelHarmonic(const Matrix& r, const Matrix& U, int n_max, int m_max);


#endif //PROYECTOTALLER_G_ACCELHARMONIC_H
