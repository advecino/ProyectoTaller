#include "VarEqn.h"
#include "global.h"
#include "IERS.h"
#include "TimeDiff.h"
#include "PrecMatrix.h"
#include "NutMatrix.h"
#include "PoleMatrix.h"
#include "GHAMatrix.h"
#include "AccelHarmonic.h"
#include "G_AccelHarmonic.h"
#include "Sat_const.h"


/*
%------------------------------------------------------------------------------
%
% VarEqn.m
%
% Purpose:
%   Computes the variational equations, i.e. the derivative of the state vector
%   and the state transition matrix
%
% Input:
%   x           Time since epoch in [s]
%   yPhi        (6+36)-dim vector comprising the state vector (y) and the
%               state transition matrix (Phi) in column wise storage order
%
% Output:
%   yPhip       Derivative of yPhi
%
% Last modified:   2015/08/12   M. Mahooti
%
%------------------------------------------------------------------------------
*/



Matrix VarEqn(
        double t,
        Matrix& yPhi,
        AuxParam& params,
        Matrix& eop
) {
    // 1) Recuperar EOP & UT1
    double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;
    IERSResult ires = IERS(eop, params.Mjd_UTC, 'l');
    x_pole   = ires.x_pole;
    y_pole   = ires.y_pole;
    UT1_UTC  = ires.UT1_UTC;
    LOD      = ires.LOD;
    dpsi     = ires.dpsi;
    deps     = ires.deps;
    dx_pole  = ires.dx_pole;
    dy_pole  = ires.dy_pole;
    TAI_UTC  = ires.TAI_UTC;

    // 2) Diferencias UT1/TAI/GPS/TT/UTC
    TimeDiffs tres = timediff(UT1_UTC, TAI_UTC);
    double UT1_TAI  = tres.UT1_TAI;
    double UTC_GPS = tres.UTC_GPS;
    double UT1_GPS = tres.UT1_GPS;
    double TT_UTC  = tres.TT_UTC;
    double GPS_UTC = tres.GPS_UTC;

    // 3) Mjd_UT1
    double Mjd_UT1 = params.Mjd_TT + (UT1_UTC - TT_UTC)/86400.0;

    // 4) Matriz de transformación E = Pole * GHA * Nut * Precession
    Matrix P = PrecMatrix(MJD_J2000, params.Mjd_TT + t/86400.0);
    Matrix N = NutMatrix(params.Mjd_TT + t/86400.0);
    Matrix Tmat = N * P;
    Matrix E = PoleMatrix(x_pole, y_pole)
               * GHAMatrix(Mjd_UT1)
               * Tmat;

    // 5) Extraer r(3×1), v(3×1) y Φ(6×6)
    Matrix r = yPhi.getSubMatrix(1,3,1,1);
    Matrix v = yPhi.getSubMatrix(4,6,1,1);
    Matrix Phi(6,6);
    for(int j=1; j<=6; ++j) {
        // columnas de Phi almacenadas columna-mayor en yPhi
        Matrix col = yPhi.getSubMatrix(6*j+1, 6*j+6, 1, 1);
        Phi.setColumn(j, col);
    }

    // 6) Aceleración y su jacobiano G = ∂a/∂r
    Matrix a = AccelHarmonic(r, E, params.n, params.m);
    Matrix Gmat = G_AccelHarmonic(r, E, params.n, params.m);

    // 7) Construir df/dy (6×6)
    Matrix dfdy(6,6);
    // [ 0    I ]
    // [ Gmat 0 ]
    for(int i=1;i<=3;++i) {
        for(int j=1;j<=3;++j) {
            dfdy(i, j)     = 0.0;            // dv/dr
            dfdy(i, j+3)   = (i==j ? 1.0 : 0.0); // dv/dv
            dfdy(i+3, j)   = Gmat(i,j);     // da/dr
            dfdy(i+3, j+3) = 0.0;           // da/dv
        }
    }

    // 8) dΦ/dt = dfdy * Φ
    Matrix Phip = dfdy * Phi;

    // 9) Construir salida yPhip (42×1)
    Matrix yPhip(42,1);
    // dr/dt = v, dv/dt = a
    for(int i=1; i<=3; ++i) {
        yPhip(i,1)   = v(i,1);
        yPhip(i+3,1) = a(i,1);
    }
    // dΦ/dt columna-mayor
    for(int j=1; j<=6; ++j) {
        for(int i=1; i<=6; ++i) {
            yPhip(6*j + i, 1) = Phip(i,j);
        }
    }

    return yPhip;
}
