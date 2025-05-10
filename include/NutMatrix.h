//
// Created by adria on 11/05/2025.
//

#ifndef PROYECTOTALLER_NUTMATRIX_H
#define PROYECTOTALLER_NUTMATRIX_H


#include "Matrix.h"
#include "MeanObliquity.h"
#include "NutAngles.h"
#include "R_x.h"
#include "R_z.h"

/**
 * @brief Transformation from mean to true equator and equinox
 *
 * @param Mjd_TT Modified Julian Date (Terrestrial Time)
 * @return Matrix Nutation matrix
 *
 * Last modified: 2015/08/12 M. Mahooti
 */
Matrix NutMatrix(double Mjd_TT);


#endif //PROYECTOTALLER_NUTMATRIX_H
