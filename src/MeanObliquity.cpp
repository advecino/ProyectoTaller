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

#include "../include/MeanObliquity.h"
#include "../include/Sat_const.h"

double MeanObliquity(double Mjd_TT) {
    double T = (Mjd_TT - MJD_J2000) / 36525.0;

    double MOblq_arcsec = 84381.448 - (46.8150 + (0.00059 - 0.001813 * T) * T) * T;

    return MOblq_arcsec / 3600.0 * Rad;
}
