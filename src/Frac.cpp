//
// Created by advecino on 09/04/2025.
//
/*
%--------------------------------------------------------------------------
%
%  Fractional part of a number (y=x-[x])
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
 */
/**
 * @file Frac.cpp
 * @brief Implementación de la función para obtener la parte fraccionaria de un número
 */
#include "../include/Frac.h"
#include <cmath>

double Frac(double x) {
    return x - floor(x);
}