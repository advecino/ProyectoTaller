#ifndef PROYECTOTALLER_NUTMATRIX_H
#define PROYECTOTALLER_NUTMATRIX_H


#include "Matrix.h"
#include "MeanObliquity.h"
#include "NutAngles.h"
#include "R_x.h"
#include "R_z.h"

/**
 * @brief Matriz de transformación por nutación (ecuador/equinoctio medio → verdadero).
 * @param Mjd_TT Fecha juliana modificada (Tiempo Terrestre).
 * @return Matriz de nutación.
 */
Matrix NutMatrix(double Mjd_TT);


#endif //PROYECTOTALLER_NUTMATRIX_H
