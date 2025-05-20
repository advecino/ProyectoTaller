#ifndef PROYECTOTALLER_GAST_H
#define PROYECTOTALLER_GAST_H

#include "gmst.h"
#include "EqnEquinox.h"

/**
 * @brief Tiempo sidéreo aparente de Greenwich (GAST).
 * @param Mjd_UT1 Fecha juliana modificada (UT1).
 * @return GAST en radianes.
 */
double gstime(double Mjd_UT1);


#endif //PROYECTOTALLER_GAST_H
