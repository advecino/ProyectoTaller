#ifndef PROYECTOTALLER_EQNEQUINOX_H
#define PROYECTOTALLER_EQNEQUINOX_H

#include "NutAngles.h"
#include "MeanObliquity.h"

/**
 * @brief Calcula la ecuación de los equinoccios.
 *
 * @param Mjd_TT Fecha juliana modificada (Tiempo Terrestre).
 * @return Ecuación de los equinoccios (dψ·cos(ε)), en radianes.
 */
double EqnEquinox(double Mjd_TT);

#endif //PROYECTOTALLER_EQNEQUINOX_H
