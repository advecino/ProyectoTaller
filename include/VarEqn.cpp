//
// Created by adria on 12/05/2025.
//

#include <stdexcept>
#include <iostream>
#include "VarEqn.h"

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


Matrix VarEqn(double x, const Matrix& yPhi, const Matrix& eopdata) {

    // Validate input dimensions
    if (yPhi.getFilas() != 42 || yPhi.getColumnas() != 1) {
        throw std::invalid_argument("yPhi must be a 42x1 vector");
    }

    double Mjd_UTC=49746.116354;
    std::string s = "l";
    double Mjd_TT= 49746.117062;
    // IERS Earth orientation parameters
    IERSResult iers = IERS(eopdata, Mjd_UTC, s);

    // Time differences
    TimeDiffs td = timediff(iers.UT1_UTC, iers.TAI_UTC);

    // Modified Julian Date UT1
    double Mjd_UT1 = Mjd_TT + (iers.UT1_UTC - td.TT_UTC)/86400.0;

    // Transformation matrices
    Matrix P = PrecMatrix(MJD_J2000, Mjd_TT + x/86400.0);
    Matrix N = NutMatrix(Mjd_TT + x/86400.0);
    Matrix T = N * P;
    Matrix E = PoleMatrix(iers.x_pole, iers.y_pole) * GHAMatrix(Mjd_UT1) * T;

    // State vector components
    Matrix r = yPhi.getSubMatrix(0+1, 2+1, 0+1, 0+1);
    Matrix v = yPhi.getSubMatrix(3+1, 5+1, 0+1, 0+1);

    // State transition matrix (reshape from vector to 6x6 matrix)
    Matrix Phi(6, 6);
    for (int j = 1; j < 6; j++) {
        for (int i = 1; i < 6; i++) {
            Phi(i, j) = yPhi(6*(j+1) + i, 1);
        }
    }

    int n=20;
    int m =20;
    // Acceleration and gradient
    Matrix a = AccelHarmonic(r, E, n, m);
    Matrix G = G_AccelHarmonic(r, E, n, m);

    // Time derivative of state transition matrix
    Matrix dfdy(6, 6);
    for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 3; j++) {
            dfdy(i, j) = 0.0;                 // dv/dr(i,j)
            dfdy(i+3, j) = G(i, j);           // da/dr(i,j)
            dfdy(i, j+3) = (i == j) ? 1 : 0;  // dv/dv(i,j)
            dfdy(i+3, j+3) = 0.0;             // da/dv(i,j)
        }
    }

    Matrix Phip = dfdy * Phi;

    // Derivative of combined state vector and state transition matrix
    Matrix yPhip(42, 1);
    for (int i = 1; i <= 3; i++) {
        yPhip(i, 1) = v(i, 1);               // dr/dt(i)
        yPhip(i+3, 1) = a(i, 1);             // dv/dt(i)
    }

    for (int j = 1; j <= 6; j++) {
        for (int i = 1; i <= 6; i++) {
            yPhip(6*(j-1) + i, 1) = Phip(i, j); // dPhi/dt(i,j)
        }
    }

    return yPhip;
}