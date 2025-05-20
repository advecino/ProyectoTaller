#include "../include/Mjday_TDB.h"
/*
%--------------------------------------------------------------------------
%
% Mjday_TDB: Computes the Modified Julian Date for barycentric dynamical
%            time
%
%  Inputs:
%    Mjd_TT      - Modified julian date (TT)
%
%  Output:
%    Mjd_TDB     - Modified julian date (TDB)
%
% Reference:
% Vallado D. A; Fundamentals of Astrodynamics and Applications; McGraw-Hill;
% New York; 3rd edition(2007).
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------*/
#include <cmath>

double Mjday_TDB(double Mjd_TT)
{
    // Compute Julian Centuries of TT
    double T_TT = (Mjd_TT - 51544.5) / 36525.0;

    // Compute Modified Julian Date of TDB
    double sin_terms =
            0.001658 * sin(628.3076 * T_TT + 6.2401) +
            0.000022 * sin(575.3385 * T_TT + 4.2970) +
            0.000014 * sin(1256.6152 * T_TT + 6.1969) +
            0.000005 * sin(606.9777 * T_TT + 4.0212) +
            0.000005 * sin(52.9691 * T_TT + 0.4444) +
            0.000002 * sin(21.3299 * T_TT + 5.5431) +
            0.000010 * sin(628.3076 * T_TT + 4.2490);

    double Mjd_TDB = Mjd_TT + sin_terms / 86400.0;

    return Mjd_TDB;
}

