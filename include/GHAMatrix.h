//
// Created by adria on 11/05/2025.
//

#ifndef PROYECTOTALLER_GHAMATRIX_H
#define PROYECTOTALLER_GHAMATRIX_H

#include "Matrix.h"
#include "gast.h"
#include "R_z.h"


/**
 * @brief Transformation from true equator and equinox to Earth equator and Greenwich meridian system
 *
 * @param Mjd_UT1 Modified Julian Date UT1
 * @return Matrix Greenwich Hour Angle matrix
 */
Matrix GHAMatrix(double Mjd_UT1);


#endif //PROYECTOTALLER_GHAMATRIX_H
