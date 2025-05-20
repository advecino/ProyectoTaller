#include <stdexcept>
#include "../include/IERS.h"
#include "../include/global.h"
#include "../include/Sat_const.h"
#include <cmath>


/*
%--------------------------------------------------------------------------
%
% IERS: Management of IERS time and polar motion data
%
% Last modified:   2018/02/01   M. Mahooti
%
%--------------------------------------------------------------------------*/


IERSResult IERS(Matrix& eop, double Mjd_UTC, char interp) {
    if (interp!='l' && interp!='n') interp='n';

    double mjd_floor = std::floor(Mjd_UTC);
    int cols = eop.getColumnas();
    int idx = -1;
    for (int j = 1; j <= cols; ++j) {
        if (std::floor(eop(4, j)) == mjd_floor) { idx = j; break; }
    }
    if (idx < 1 || idx > cols)
        throw std::out_of_range("Mjd_UTC fuera de rango en eop");

    if (interp=='l' && idx==cols) interp='n';

    Matrix pre = eop.getSubMatrix(1,13,idx,idx);

    Matrix nxt(0,0);
    double fixf = 0.0;
    if (interp=='l') {
        nxt   = eop.getSubMatrix(1,13,idx+1,idx+1);
        fixf  = Mjd_UTC - mjd_floor;
    }

    auto lerp = [&](int row)->double {
        double v0 = pre(row,1);
        if (interp=='l') {
            double v1 = nxt(row,1);
            return v0 + (v1 - v0)*fixf;
        } else {
            return v0;
        }
    };

    IERSResult res;
    res.x_pole  = lerp(5)/Arcs;
    res.y_pole  = lerp(6)/Arcs;
    res.UT1_UTC = lerp(7);
    res.LOD     = lerp(8);
    res.dpsi    = lerp(9)/Arcs;
    res.deps    = lerp(10)/Arcs;
    res.dx_pole = lerp(11)/Arcs;
    res.dy_pole = lerp(12)/Arcs;
    res.TAI_UTC = lerp(13);

    return res;
}


