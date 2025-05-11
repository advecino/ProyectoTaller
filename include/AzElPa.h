//
// Created by adria on 12/05/2025.
//

#ifndef PROYECTOTALLER_AZELPA_H
#define PROYECTOTALLER_AZELPA_H

#include "Matrix.h"


struct AzElPaResult {
    double Az;     // Azimuth [rad]
    double El;     // Elevation [rad]
    Matrix dAds;   // Partials of azimuth w.r.t. s (1x3)
    Matrix dEds;   // Partials of elevation w.r.t. s (1x3)
    AzElPaResult() : dAds(1,3), dEds(1,3) {};
};

/**
 * @brief Computes azimuth, elevation and partials from local tangent coordinates
 *
 * @param s Topocentric local tangent coordinates (East-North-Zenith frame) (3x1)
 * @return AzElPaResult Struct containing results
 *
 * @last_modified 2015/08/12 M. Mahooti
 */
AzElPaResult AzElPa(const Matrix& s);

#endif //PROYECTOTALLER_AZELPA_H
