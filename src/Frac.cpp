/*
%--------------------------------------------------------------------------
%
%  Fractional part of a number (y=x-[x])
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
 */

#include "../include/Frac.h"
#include <cmath>

double Frac(double x) {
    return x - floor(x);
}