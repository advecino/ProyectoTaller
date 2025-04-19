//
// Created by adria on 20/04/2025.
//

#ifndef PROYECTOTALLER_CHEB3D_H
#define PROYECTOTALLER_CHEB3D_H

#include "Matrix.h"
#include <stdexcept>


// Declaración de la función auxiliar
double ChebT(double x, int n);


/**
 * @brief Aproximación de Chebyshev para vectores 3D
 * @param t Tiempo para evaluar
 * @param N Número de coeficientes
 * @param Ta Inicio del intervalo
 * @param Tb Fin del intervalo
 * @param Cx Coeficientes para x (ordenados de T0 a TN-1)
 * @param Cy Coeficientes para y (ordenados de T0 a TN-1)
 * @param Cz Coeficientes para z (ordenados de T0 a TN-1)
 * @return Vector 3D aproximado
 */


Matrix Cheb3D(double t, int N, double Ta, double Tb,
              const Matrix& Cx, const Matrix& Cy, const Matrix& Cz);



#endif //PROYECTOTALLER_CHEB3D_H
