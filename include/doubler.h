#ifndef PROYECTOTALLER_DOUBLER_H
#define PROYECTOTALLER_DOUBLER_H


#include "Matrix.h"
#include "Sat_const.h"

struct DoubleRResult {
    Matrix r2;           ///< Vector de posición en t2
    Matrix r3;           ///< Vector de posición en t3
    double f1;           ///< Función auxiliar f1
    double f2;           ///< Función auxiliar f2
    double q1;           ///< Norma de [f1 f2]
    double magr1;        ///< Módulo de r1
    double magr2;        ///< Módulo de r2
    double a;            ///< Semieje mayor de la órbita
    double deltae32;     ///< Diferencia de anomalías excéntricas entre r3 y r2
    DoubleRResult() : r2(3,1), r3(3,1) {};
};

/**
 * @brief Determinación orbital usando el método Double-R.
 *
 * @param cc1, cc2 Parámetros intermedios relacionados con ángulos y distancias.
 * @param magrsite1, magrsite2 Módulos de los vectores de sitio.
 * @param magr1in, magr2in Estimaciones iniciales de r1 y r2.
 * @param los1, los2, los3 Vectores de dirección desde los sitios de observación.
 * @param rsite1, rsite2, rsite3 Posiciones de los sitios de observación.
 * @param t1, t3 Tiempos t1 y t3 (segundos respecto a t2).
 * @param direct Dirección del movimiento ('y' directo, 'n' retrógrado).
 * @return Estructura con resultados orbitales intermedios.
 */
DoubleRResult doubler(double cc1, double cc2,
                      double magrsite1, double magrsite2,
                      double magr1in, double magr2in,
                      const Matrix& los1, const Matrix& los2, const Matrix& los3,
                      const Matrix& rsite1, const Matrix& rsite2, const Matrix& rsite3,
                      double t1, double t3, char direct);


#endif //PROYECTOTALLER_DOUBLER_H
