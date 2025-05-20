#ifndef PROYECTOTALLER_HGIBBS_H
#define PROYECTOTALLER_HGIBBS_H

#include <string>
#include "Matrix.h"
#include "Sat_const.h"
#include "unit.h"
#include "angl.h"
#include "gibbs.h"



/**
 * @brief Método de Herrick-Gibbs para determinación orbital.
 *
 * @param r1 Vector de posición 1 [m].
 * @param r2 Vector de posición 2 [m].
 * @param r3 Vector de posición 3 [m].
 * @param Mjd1 Fecha juliana de la primera observación.
 * @param Mjd2 Fecha juliana de la segunda observación.
 * @param Mjd3 Fecha juliana de la tercera observación.
 * @return Resultados del método (velocidad en r2, ángulos, etc.).
 */
GibbsResult hgibbs(const Matrix& r1, const Matrix& r2, const Matrix& r3,double Mjd1, double Mjd2, double Mjd3);


#endif //PROYECTOTALLER_HGIBBS_H
