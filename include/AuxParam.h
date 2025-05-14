//
// Created by adria on 14/05/2025.
//

#ifndef PROYECTOTALLER_AUXPARAM_H
#define PROYECTOTALLER_AUXPARAM_H


#include "Matrix.h"

/**
 * @brief Parámetros auxiliares pasados a VarEqn
 */
struct AuxParam {
    double Mjd_UTC;  ///< Modified Julian Date UTC
    double Mjd_TT;   ///< Modified Julian Date TT
    int    n;        ///< Grado máximo para AccelHarmonic
    int    m;        ///< Orden máximo para AccelHarmonic

    bool sun;            ///< Incluir perturbación solar
    bool moon;           ///< Incluir perturbación lunar
    bool planets;        ///< Incluir perturbaciones planetarias
};


extern AuxParam AuxParamGlob;


#endif //PROYECTOTALLER_AUXPARAM_H
