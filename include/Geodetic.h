//
// Created by adria on 11/05/2025.
//

#ifndef PROYECTOTALLER_GEODETIC_H
#define PROYECTOTALLER_GEODETIC_H

#include "Matrix.h"
#include "Sat_const.h"


struct GeodeticCoords {
    double longitude;  // Longitud en radianes
    double latitude;   // Latitud en radianes
    double altitude;   // Altitud en metros

    // Constructor por defecto
    GeodeticCoords() : longitude(0.0), latitude(0.0), altitude(0.0) {}

    // Constructor con parámetros
    GeodeticCoords(double lon, double lat, double h)
            : longitude(lon), latitude(lat), altitude(h) {}
};

/**
 * @brief Calcula coordenadas geodésicas a partir de un vector de posición
 *
 * @param r Vector de posición [m] (Matrix 3x1)
 * @param R_equ Radio ecuatorial de la Tierra [m]
 * @param f Achatamiento de la Tierra
 * @return GeodeticCoords Estructura con coordenadas geodésicas
 */
GeodeticCoords Geodetic(const Matrix& r);


#endif //PROYECTOTALLER_GEODETIC_H
