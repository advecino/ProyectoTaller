#ifndef PROYECTOTALLER_ACCEL_H
#define PROYECTOTALLER_ACCEL_H

#include "Matrix.h"
#include "AuxParam.h"
#include "JPL_Eph_DE430.h"

/**
 * @brief Calcula la aceleración total que actúa sobre un satélite en órbita terrestre.
 * @param x Tiempo transcurrido desde la época de referencia, en segundos.
 * @param Y Vector de estado (6x1) del satélite en el sistema ICRF/EME2000. Contiene posición [km] y velocidad [km/s].
 * @param params Estructura con parámetros auxiliares, incluyendo fecha MJD, coeficientes armónicos, y opciones activadas.
 * @param eop Matriz de parámetros de orientación terrestre (EOP) del IERS, usada para correcciones de precesión/nutación y rotación terrestre.
 * @return Vector (6x1) de derivadas del estado: dY(1:3) = velocidad, dY(4:6) = aceleración total.
 */
Matrix Accel(double x,const Matrix& Y, AuxParam& params,Matrix& eop);

#endif //PROYECTOTALLER_ACCEL_H
