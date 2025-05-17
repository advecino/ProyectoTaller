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
        const Matrix& yPhi,
        const AuxParam& params,
        Matrix& eop
) {
    // 1) Recuperar EOP & UT1 (sin interpolación lineal para evitar out-of-range)
    IERSResult ires = IERS(eop, params.Mjd_UTC + t/86400.0);
    double x_pole  = ires.x_pole;
    double y_pole  = ires.y_pole;
    double UT1_UTC = ires.UT1_UTC;
    double LOD     = ires.LOD;
    double dpsi    = ires.dpsi;
    double deps    = ires.deps;
    double dx_pole = ires.dx_pole;
    double dy_pole = ires.dy_pole;
    double TAI_UTC = ires.TAI_UTC;

    // 2) UT1/TAI/GPS/TT/UTC
    TimeDiffs td = timediff(UT1_UTC, TAI_UTC);
    double TT_UTC = td.TT_UTC;

    // 3) Fecha UT1 en días (MJD)
    double Mjd_UT1 = params.Mjd_TT + (UT1_UTC - TT_UTC)/86400.0;

    // 4) Construir matriz de transformación E
    double Mjd_TT_now = params.Mjd_TT + t/86400.0;
    Matrix P = PrecMatrix(MJD_J2000, Mjd_TT_now);
    Matrix N = NutMatrix(Mjd_TT_now);
    Matrix Tmat = N * P;
    Matrix E = PoleMatrix(x_pole, y_pole)
               * GHAMatrix(Mjd_UT1)
               * Tmat;

    // 5) Extraer r(3×1), v(3×1), Φ(6×6) de yPhi (42×1)
    Matrix r = yPhi.getSubMatrix(1,3, 1,1);
    Matrix v = yPhi.getSubMatrix(4,6, 1,1);
    Matrix Phi(6,6);
    for(int j=1; j<=6; ++j) {
        // columnas de Φ en orden columna-mayor
        Matrix col = yPhi.getSubMatrix(6*j+1, 6*j+6, 1,1);
        Phi.setColumn(j, col);
    }

    // 6) Aceleración y jacobiano
    Matrix a    = AccelHarmonic(    r, E, params.n, params.m);
    Matrix Gmat = G_AccelHarmonic(r, E, params.n, params.m);

    // 7) Construir df/dy = [ 0 I ; Gmat 0 ]
    Matrix dfdy(6,6);
    // llenamos a 0
    for(int i=1;i<=6;++i) for(int j=1;j<=6;++j) dfdy(i,j)=0.0;
    // identidad en el tope
    for(int i=1;i<=3;++i) dfdy(i,i+3)=1.0;
    // Gmat en la parte inferior izquierda
    for(int i=1;i<=3;++i)
        for(int j=1;j<=3;++j)
            dfdy(i+3,j)=Gmat(i,j);

    // 8) dΦ/dt = dfdy * Φ
    Matrix Phip = dfdy * Phi;

    // 9) Construir yPhip (42×1)
    Matrix yPhip(42,1);
    // dr/dt = v, dv/dt = a
    for(int i=1;i<=3;++i) {
        yPhip(i,1)   = v(i,1);
        yPhip(i+3,1) = a(i,1);
    }
    // luego dΦ/dt en orden columna-mayor
    for(int j=1;j<=6;++j)
        for(int i=1;i<=6;++i)
            yPhip(6*j + i,1) = Phip(i,j);

    return yPhip;
}

