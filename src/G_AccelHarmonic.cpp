#include "../include/G_AccelHarmonic.h"
#include "../include/AccelHarmonic.h"

/*
%--------------------------------------------------------------------------
%
% G_AccelHarmonic.m
%
% Purpose:
%   Computes the gradient of the Earth's harmonic gravity field
%
% Inputs:
%   r           Satellite position vector in the true-of-date system
%   U           Transformation matrix to body-fixed system
%   n           Gravity model degree
%   m 			Gravity model order
%
% Output:
%   G    		Gradient (G=da/dr) in the true-of-date system
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------*/

Matrix G_AccelHarmonic(const Matrix& r, const Matrix& U, int n_max, int m_max) {
    double d = 1.0;   // Position increment [m]

    Matrix G = Matrix(3,3);
    Matrix dr = Matrix(3,1);

    // Gradient
    for(int i = 1; i <= 3; i++) {  // Cambiado a índices basados en 1
        // Set offset in i-th component of the position vector
        dr(1,1) = 0.0;
        dr(2,1) = 0.0;
        dr(3,1) = 0.0;
        dr(i,1) = d;

        // Acceleration difference
        Matrix r_plus = dr * 0.5;
        r_plus = r_plus + r;
        Matrix r_minus = dr * 0.5;
        r_minus = -1 * r_minus;
        r_minus = r_minus + r;

        Matrix da = AccelHarmonic(r_plus, U, n_max, m_max) -
                    AccelHarmonic(r_minus, U, n_max, m_max);

        // Derivative with respect to i-th axis
        G(1,i) = da(1,1)/d;
        G(2,i) = da(2,1)/d;
        G(3,i) = da(3,1)/d;
    }
    return G;
}



