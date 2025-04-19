/**
 * @file Position.h
 * @brief Calcula el vector de posición a partir de coordenadas geodésicas
 */

#ifndef PROYECTOTALLER_POSITION_H
#define PROYECTOTALLER_POSITION_H


#include "Matrix.h"

/**
 * @brief Calcula el vector de posición desde coordenadas geodésicas
 * @param lon Longitud en radianes
 * @param lat Latitud en radianes
 * @param h Altitud en metros
 * @param R_equ Radio ecuatorial de la Tierra
 * @param f Achatamiento de la Tierra
 * @return Vector de posición (Matrix 3x1)
 */
Matrix Position(double lon, double lat, double h, double R_equ, double f);



#endif //PROYECTOTALLER_POSITION_H
