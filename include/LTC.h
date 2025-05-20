#ifndef PROYECTOTALLER_LTC_H
#define PROYECTOTALLER_LTC_H

#include "Matrix.h"
#include "R_y.h"
#include "R_z.h"

/**
 * @brief Transformación al sistema de coordenadas locales tangentes.
 * @param lon Longitud geodésica este [rad].
 * @param lat Latitud geodésica [rad].
 * @return Matriz de rotación al sistema local (Este-Norte-Zenit).
 */
Matrix LTC(double lon, double lat);

#endif //PROYECTOTALLER_LTC_H
