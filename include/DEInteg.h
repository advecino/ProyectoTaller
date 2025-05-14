//
// Created by adria on 12/05/2025.
//

#ifndef PROYECTOTALLER_DEINTEG_H
#define PROYECTOTALLER_DEINTEG_H


#include <vector>
#include <functional>
#include "Matrix.h"

/**
 * @brief Integrador de paso variable estilo Dormand–Prince (Ode45).
 *
 * @param func    Función que calcula dy/dt = f(t, y). Firma: f(t, y_in, y_dot_out).
 * @param t0      Tiempo inicial [s].
 * @param tout    Vector de tiempos de salida (monótonamente creciente).
 * @param relerr  Tolerancia relativa.
 * @param abserr  Tolerancia absoluta.
 * @param y0      Condición inicial (n_eqn×1).
 * @return        Solución en cada tout: matriz de tamaño (n_eqn × tout.size()).
 */
Matrix DEInteg(
        const std::function<void(double, const Matrix&, Matrix&)>& func,
        double t0,
        const std::vector<double>& tout,
        double relerr,
        double abserr,
        const Matrix& y0
);

#endif //PROYECTOTALLER_DEINTEG_H
