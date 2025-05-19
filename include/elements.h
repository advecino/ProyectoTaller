//
// Created by adria on 10/05/2025.
//

#ifndef PROYECTOTALLER_ELEMENTS_H
#define PROYECTOTALLER_ELEMENTS_H


#include <cmath>
#include "Matrix.h"
#include "SAT_Const.h"

struct KeplerianElements {
    double p;
    double a;
    double e;
    double i;
    double Omega;
    double omega;
    double M;
};

/**
 * @brief Osculating Keplerian elements from r and v for elliptic, inclined orbits
 * @throws std::invalid_argument if orbit is circular or equatorial
 */
KeplerianElements elements(const Matrix& r, const Matrix& v);

#endif //PROYECTOTALLER_ELEMENTS_H
