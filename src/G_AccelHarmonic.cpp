#include <iostream>
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

Matrix G_AccelHarmonic(Matrix& r, Matrix& U, int n_max, int m_max) {
    const double delta = 1.0;

    // Sanity checks
    if (r.getFilas() != 3 || r.getColumnas() != 1)
        throw std::invalid_argument("r must be 3×1");
    if (U.getFilas() != 3 || U.getColumnas() != 3)
        throw std::invalid_argument("U must be 3×3");
    if (n_max < 0 || m_max < 0 || m_max > n_max)
        throw std::invalid_argument("Invalid degree/order in G_AccelHarmonic");

    Matrix G(3, 3);
    // for each coordinate direction i
    for (int i = 1; i <= 3; ++i) {
        // build offset vector dr = δ * e_i
        Matrix dr(3,1);
        dr(i,1) = delta;

        // central points
        Matrix r_plus  = r + dr * 0.5;
        Matrix r_minus = r - dr * 0.5;

        // compute accelerations
        Matrix a_plus  = AccelHarmonic(r_plus,  U, n_max, m_max);
        Matrix a_minus = AccelHarmonic(r_minus, U, n_max, m_max);

        // finite‐difference derivative for column i
        Matrix da = a_plus - a_minus;
        for (int k = 1; k <= 3; ++k) {
            G(k, i) = da(k,1) / delta;
        }
    }

    return G;
}



