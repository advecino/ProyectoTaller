#ifndef PROYECTOTALLER_POSITION_H
#define PROYECTOTALLER_POSITION_H


#include "Matrix.h"

/**
 * @brief Vector de posición desde coordenadas geodésicas.
 * @param lon Longitud [rad].
 * @param lat Latitud [rad].
 * @param h Altitud [m].
 * @param R_equ Radio ecuatorial de la Tierra [m].
 * @param f Achatamiento de la Tierra.
 * @return Vector de posición (3x1) en coordenadas cartesianas [m].
 */
Matrix Position(double lon, double lat, double h, double R_equ, double f);



#endif //PROYECTOTALLER_POSITION_H
