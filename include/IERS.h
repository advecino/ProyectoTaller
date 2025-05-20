#ifndef PROYECTOTALLER_IERS_H
#define PROYECTOTALLER_IERS_H



#include "Matrix.h"


/**
 * @brief Consulta de parámetros EOP del IERS.
 */

struct IERSResult {
    double x_pole;   ///< Coordenada del polo [rad]
    double y_pole;   ///< Coordenada del polo [rad]
    double UT1_UTC;  ///< Diferencia UT1-UTC [s]
    double LOD;      ///< Longitud del día [s]
    double dpsi;     ///< Nutación en longitud [rad]
    double deps;     ///< Nutación en oblicuidad [rad]
    double dx_pole;  ///< Desplazamiento del polo [rad]
    double dy_pole;  ///< Desplazamiento del polo [rad]
    double TAI_UTC;  ///< Diferencia TAI-UTC [s]
};

/**
 * @brief Datos IERS de tiempo y movimiento polar.
 * @param eop Matriz de parámetros de orientación terrestre (mínimo 13 filas).
 * @param Mjd_UTC Fecha juliana modificada (UTC).
 * @param interp 'l' para interpolación lineal, 'n' sin interpolación.
 * @return Estructura con resultados: coordenadas del polo, UT1-UTC, TAI-UTC, etc.
 */
IERSResult IERS(Matrix& eop, double Mjd_UTC, char interp='n');
#endif //PROYECTOTALLER_IERS_H
