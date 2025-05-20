#ifndef PROYECTOTALLER_GLOBAL_H
#define PROYECTOTALLER_GLOBAL_H

#include <vector>
#include "Matrix.h"

/**
 * @brief Declaraciones globales usadas en el modelo dinámico y transformaciones.
 */

extern double Cnm[300][300];  ///< Coeficientes armónicos Cnm
extern double Snm[300][300];  ///< Coeficientes armónicos Snm

extern Matrix PC;             ///< Parámetros polinomiales de interpolación (precisión)
extern Matrix eopdata;        ///< Datos de orientación terrestre (EOP)


#endif //PROYECTOTALLER_GLOBAL_H
