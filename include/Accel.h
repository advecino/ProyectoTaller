//
// Created by adria on 12/05/2025.
//

#ifndef PROYECTOTALLER_ACCEL_H
#define PROYECTOTALLER_ACCEL_H

#include "Matrix.h"
#include "AuxParam.h"
#include "JPL_Eph_DE430.h"

/**
 * @brief Calcula la derivada del estado [v; a] para un satélite en órbita terrestre.
 *
 * @param x      Tiempo desde la época en segundos.
 * @param Y      Vector de estado [r; v] (6×1) en ICRF/EME2000.
 * @param params Parámetros auxiliares (fechas MJD, grado/orden de gravidad, flags).
 * @param eop    Matriz de parámetros EOP para IERS (13 filas × N columnas).
 * @return dY    Vector (6×1) donde dY(1:3)=v, dY(4:6)=a.
 */
Matrix Accel(
        double x,
        const Matrix& Y,
        AuxParam& params,
        Matrix& eop
);

#endif //PROYECTOTALLER_ACCEL_H
