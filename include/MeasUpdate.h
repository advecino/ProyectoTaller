//
// Created by adria on 20/04/2025.
//

#ifndef PROYECTOTALLER_MEASUPDATE_H
#define PROYECTOTALLER_MEASUPDATE_H


#include "Matrix.h"

/**
 * @brief Actualización de medición del filtro de Kalman
 * @param x Vector de estado (se actualiza in-place)
 * @param z Vector de medición
 * @param g Función de predicción de medición
 * @param s Desviación estándar de las mediciones
 * @param G Matriz de diseño de medición
 * @param P Matriz de covarianza (se actualiza in-place)
 * @param n Dimensión del estado
 * @return K Ganancia de Kalman
 */
Matrix MeasUpdate(Matrix& x,  Matrix& z,  Matrix& g, Matrix& s,  Matrix& G, Matrix& P, int n,Matrix& K_out);

#endif //PROYECTOTALLER_MEASUPDATE_H
