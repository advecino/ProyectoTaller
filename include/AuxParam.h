#ifndef PROYECTOTALLER_AUXPARAM_H
#define PROYECTOTALLER_AUXPARAM_H


#include "Matrix.h"

/**
 * @brief Parámetros auxiliares para propagación y perturbaciones.
 */
struct AuxParam {
    double Mjd_UTC;  ///< Fecha juliana modificada (UTC)
    double Mjd_TT;   ///< Fecha juliana modificada (TT)
    int    n;        ///< Grado máximo para el campo armónico
    int    m;        ///< Orden máximo para el campo armónico

    bool sun;        ///< Incluir perturbación solar
    bool moon;       ///< Incluir perturbación lunar
    bool planets;    ///< Incluir perturbaciones planetarias
};


extern AuxParam AuxParamGlob;


#endif //PROYECTOTALLER_AUXPARAM_H
