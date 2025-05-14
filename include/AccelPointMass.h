//
// Created by adria on 20/04/2025.
//

#ifndef PROYECTOTALLER_ACCELPOINTMASS_H
#define PROYECTOTALLER_ACCELPOINTMASS_H


#include "Matrix.h"

/**
 * @brief Calcula la aceleración perturbacional debido a una masa puntual
 * @param r Vector de posición del satélite (Matrix 3x1)
 * @param s Vector de posición de la masa puntual (Matrix 3x1)
 * @param GM Coeficiente gravitacional de la masa puntual
 * @return Vector de aceleración (Matrix 3x1)
 */
Matrix AccelPointMass(Matrix& r, const Matrix& s, double GM);


#endif //PROYECTOTALLER_ACCELPOINTMASS_H
