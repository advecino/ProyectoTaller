//
// Created by advecino on 10/04/2025.
//

#ifndef PROYECTOTALLER_LEGENDRE_H
#define PROYECTOTALLER_LEGENDRE_H
#include "Matrix.h"

/**
 * @brief Calcula los polinomios de Legendre asociados y sus derivadas
 * @param n Grado máximo
 * @param m Orden máximo
 * @param fi Ángulo en radianes
 * @param pnm Matriz de salida para los polinomios (n+1 x m+1)
 * @param dpnm Matriz de salida para las derivadas (n+1 x m+1)
 */
void Legendre(int n, int m, double fi, Matrix& pnm, Matrix& dpnm);

#endif //PROYECTOTALLER_LEGENDRE_H
