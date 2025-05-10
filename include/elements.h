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

KeplerianElements elements(const Matrix& y);

#endif //PROYECTOTALLER_ELEMENTS_H
