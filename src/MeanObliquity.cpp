/*%--------------------------------------------------------------------------
%
% MeanObliquity.m
%
% Purpose:
%   Computes the mean obliquity of the ecliptic
%
% Input:
%   Mjd_TT    Modified Julian Date (Terrestrial Time)
%
% Output:
%   MOblq     Mean obliquity of the ecliptic [rad]
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------*/


/**
 * @file MeanObliquity.cpp
 * @brief Implementación del cálculo de la oblicuidad media de la eclíptica
 */

#include "../include/MeanObliquity.h"
#include "../include/Sat_const.h"

double MeanObliquity(double Mjd_TT) {
    // Calcular tiempo en siglos julianos desde J2000
    double T = (Mjd_TT - MJD_J2000) / 36525.0;

    // Calcular oblicuidad media (en segundos de arco)
    double MOblq_arcsec = 84381.448 - (46.8150 + (0.00059 - 0.001813 * T) * T) * T;

    // Convertir a radianes y devolver
    return MOblq_arcsec / 3600.0 * Rad;
}
