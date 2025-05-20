#include <cmath>
#include "../include/gast.h"
/*
%--------------------------------------------------------------------------
%
% GAST.m
%
% Purpose:
%   Greenwich Apparent Sidereal Time
%
% Input:
%   Mjd_UT1   Modified Julian Date UT1
%
% Output:
%   gstime    GAST in [rad]
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------*/

double gstime(double Mjd_UT1) {
    double EqnEquino = EqnEquinox(Mjd_UT1);
    double gms = gmst(Mjd_UT1);
    return fmod(gms + EqnEquino, 2* M_PI);
}