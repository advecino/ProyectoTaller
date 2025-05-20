#ifndef PROYECTOTALLER_CHEB3D_H
#define PROYECTOTALLER_CHEB3D_H

#include "Matrix.h"
#include <stdexcept>


/**
 * @brief Polinomio de Chebyshev de primer tipo Tn(x).
 * @param x Punto de evaluación.
 * @param n Grado del polinomio.
 * @return Valor de Tn(x).
 */
double ChebT(double x, int n);


/**
 * @brief Aproximación de Chebyshev de vectores tridimensionales.
 *
 * @param t Tiempo en el que se evalúa la aproximación.
 * @param N Número de coeficientes.
 * @param Ta Inicio del intervalo.
 * @param Tb Fin del intervalo.
 * @param Cx Coeficientes para la componente x (T0 a Tn-1).
 * @param Cy Coeficientes para la componente y (T0 a Tn-1).
 * @param Cz Coeficientes para la componente z (T0 a Tn-1).
 * @return Vector aproximado (3x1).
 */
Matrix Cheb3D(double t, int N, double Ta, double Tb,const Matrix& Cx, const Matrix& Cy, const Matrix& Cz);



#endif //PROYECTOTALLER_CHEB3D_H
