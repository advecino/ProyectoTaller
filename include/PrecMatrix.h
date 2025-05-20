#ifndef PROYECTOTALLER_PRECMATRIX_H
#define PROYECTOTALLER_PRECMATRIX_H


#include "Matrix.h"
#include "R_z.h"
#include "R_y.h"

/**
 * @brief Transformación de precesión de coordenadas ecuatoriales.
 * @param Mjd_1 Época inicial (MJD TT).
 * @param Mjd_2 Época final a la que precesar (MJD TT).
 * @return Matriz de transformación por precesión.
 */
Matrix PrecMatrix(double Mjd_1, double Mjd_2);


#endif //PROYECTOTALLER_PRECMATRIX_H
