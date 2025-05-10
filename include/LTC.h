//
// Created by adria on 11/05/2025.
//

#ifndef PROYECTOTALLER_LTC_H
#define PROYECTOTALLER_LTC_H

#include "Matrix.h"
#include "R_y.h"
#include "R_z.h"

/**
 * @brief Transformation from Greenwich meridian system to local tangent coordinates
 *
 * @param lon Geodetic East longitude [rad]
 * @param lat Geodetic latitude [rad]
 * @return Matrix Rotation matrix from the Earth equator and Greenwich meridian
 *         to the local tangent (East-North-Zenith) coordinate system
 *
 * Last modified: 2015/08/12 M. Mahooti
 */
Matrix LTC(double lon, double lat);

#endif //PROYECTOTALLER_LTC_H
