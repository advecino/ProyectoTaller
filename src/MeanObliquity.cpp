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
#include <cmath>

double MeanObliquity(double Mjd_TT, double MJD_J2000, double Rad)
{
    // Tiempo en siglos julianos desde J2000
    double T = (Mjd_TT - MJD_J2000) / 36525.0;

    // Cálculo de la oblicuidad media (arcosegundos -> radianes)
    double MOblq_arcsec = 84381.448 - (46.8150 + (0.00059 - 0.001813 * T) * T) * T;

    // Convertir a radianes
    return Rad * (MOblq_arcsec / 3600.0);
}

