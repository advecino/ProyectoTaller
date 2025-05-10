//
// Created by adria on 10/05/2025.
//

#ifndef PROYECTOTALLER_HGIBBS_H
#define PROYECTOTALLER_HGIBBS_H

#include <string>
#include "Matrix.h"
#include "Sat_const.h"
#include "unit.h"
#include "angl.h"


struct HGibbsResult {
    Matrix v2;
    double theta;
    double theta1;
    double copa;
    std::string error;
    HGibbsResult() : v2(3, 1) {}
};

/**
 * @brief Implementa el método de Herrick-Gibbs para determinación orbital
 *
 * @param r1 Vector posición 1 (m)
 * @param r2 Vector posición 2 (m)
 * @param r3 Vector posición 3 (m)
 * @param Mjd1 Fecha juliana de la observación 1 (días desde 4713 a.C.)
 * @param Mjd2 Fecha juliana de la observación 2
 * @param Mjd3 Fecha juliana de la observación 3
 * @return HGibbsResult Estructura con resultados
 */
HGibbsResult hgibbs(const Matrix& r1, const Matrix& r2, const Matrix& r3,
                    double Mjd1, double Mjd2, double Mjd3);


#endif //PROYECTOTALLER_HGIBBS_H
