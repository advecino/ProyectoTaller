#ifndef PROYECTOTALLER_GMST_H
#define PROYECTOTALLER_GMST_H

#include "../include/Frac.h"

/**
 * @brief Tiempo Sidéreo Medio de Greenwich (GMST).
 * @param Mjd_UT1 Fecha juliana modificada en UT1.
 * @return GMST en radianes [0..2π].
 */
double gmst(double Mjd_UT1);


#endif //PROYECTOTALLER_GMST_H
