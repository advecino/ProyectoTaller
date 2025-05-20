#ifndef PROYECTOTALLER_GHAMATRIX_H
#define PROYECTOTALLER_GHAMATRIX_H

#include "Matrix.h"
#include "gast.h"
#include "R_z.h"


/**
 * @brief Transformación al sistema ecuatorial de la Tierra y meridiano de Greenwich.
 * @param Mjd_UT1 Fecha juliana modificada (UT1).
 * @return Matriz del ángulo horario de Greenwich.
 */
Matrix GHAMatrix(double Mjd_UT1);


#endif //PROYECTOTALLER_GHAMATRIX_H
