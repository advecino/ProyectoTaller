//
// Created by adria on 11/05/2025.
//

#ifndef PROYECTOTALLER_PRECMATRIX_H
#define PROYECTOTALLER_PRECMATRIX_H


#include "Matrix.h"
#include "R_z.h"
#include "R_y.h"

/**
 * @brief Precession transformation of equatorial coordinates
 *
 * @param Mjd_1 Epoch given (Modified Julian Date TT)
 * @param Mjd_2 Epoch to precess to (Modified Julian Date TT)
 * @return Matrix Precession transformation matrix
 *
 * Last modified: 2015/08/12 M. Mahooti
 */
Matrix PrecMatrix(double Mjd_1, double Mjd_2);


#endif //PROYECTOTALLER_PRECMATRIX_H
