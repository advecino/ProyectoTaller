//
// Created by adria on 11/05/2025.
//

#ifndef PROYECTOTALLER_DOUBLER_H
#define PROYECTOTALLER_DOUBLER_H


#include "Matrix.h"
#include "Sat_const.h"

struct DoubleRResult {
    Matrix r2;
    Matrix r3;
    double f1;
    double f2;
    double q1;
    double magr1;
    double magr2;
    double a;
    double deltae32;
    DoubleRResult() : r2(3,1), r3(3,1){};
};

/**
 * @brief Implementación del algoritmo Double-R para determinación orbital
 *
 * @param cc1
 * @param cc2
 * @param magrsite1
 * @param magrsite2
 * @param magr1in
 * @param magr2in
 * @param los1 Vector unitario dirección 1
 * @param los2 Vector unitario dirección 2
 * @param los3 Vector unitario dirección 3
 * @param rsite1 Vector sitio 1
 * @param rsite2 Vector sitio 2
 * @param rsite3 Vector sitio 3
 * @param t1 Tiempo 1
 * @param t3 Tiempo 3
 * @param direct Dirección ('y' o 'n')
 * @param GM_Earth Constante gravitacional terrestre
 * @return DoubleRResult
 */
DoubleRResult doubler(double cc1, double cc2,
                      double magrsite1, double magrsite2,
                      double magr1in, double magr2in,
                      const Matrix& los1, const Matrix& los2, const Matrix& los3,
                      const Matrix& rsite1, const Matrix& rsite2, const Matrix& rsite3,
                      double t1, double t3, char direct);


#endif //PROYECTOTALLER_DOUBLER_H
