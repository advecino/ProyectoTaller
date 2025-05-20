#ifndef PROYECTOTALLER_AZELPA_H
#define PROYECTOTALLER_AZELPA_H

#include "Matrix.h"


struct AzElPaResult {
    double Az;     ///< Azimut [rad]
    double El;     ///< Elevación [rad]
    Matrix dAds;   ///< Parciales del azimut respecto a s (1x3)
    Matrix dEds;   ///< Parciales de la elevación respecto a s (1x3)
    AzElPaResult() : dAds(1, 3), dEds(1, 3) {};
};

/**
 * @brief Calcula azimut, elevación y derivadas parciales a partir de coordenadas locales.
 * @param s Coordenadas topocéntricas en el sistema Este-Norte-Zenit (3x1).
 * @return Estructura con Az, El, dAds y dEds.
 */
AzElPaResult AzElPa(const Matrix& s);

#endif //PROYECTOTALLER_AZELPA_H
