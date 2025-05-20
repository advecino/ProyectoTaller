#ifndef PROYECTOTALLER_GEODETIC_H
#define PROYECTOTALLER_GEODETIC_H

#include "Matrix.h"
#include "Sat_const.h"


struct GeodeticCoords {
    double longitude;  ///< Longitud [rad]
    double latitude;   ///< Latitud [rad]
    double altitude;   ///< Altitud [m]

    GeodeticCoords() : longitude(0.0), latitude(0.0), altitude(0.0) {}
    GeodeticCoords(double lon, double lat, double h): longitude(lon), latitude(lat), altitude(h) {}
};

/**
 * @brief Coordenadas geodésicas desde vector de posición.
 * @param r Vector de posición cartesiano (3x1) [m].
 * @return Coordenadas geodésicas: longitud, latitud y altitud.
 */
GeodeticCoords Geodetic(const Matrix& r);


#endif //PROYECTOTALLER_GEODETIC_H
