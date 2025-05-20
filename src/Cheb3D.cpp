/*%--------------------------------------------------------------------------
%
% Chebyshev approximation of 3-dimensional vectors
%
% Inputs:
%     N       Number of coefficients
%     Ta      Begin interval
%     Tb      End interval
%     Cx      Coefficients of Chebyshev polyomial (x-coordinate)
%     Cy      Coefficients of Chebyshev polyomial (y-coordinate)
%     Cz      Coefficients of Chebyshev polyomial (z-coordinate)
%
% Last modified:   2018/01/27   M. Mahooti
%
%--------------------------------------------------------------------------*/

#include "../include/Cheb3D.h"
#include <cmath>

Matrix Cheb3D(double t, int N, double Ta, double Tb,
              const Matrix& Cx, const Matrix& Cy, const Matrix& Cz)
{
    if (t < Ta || t > Tb) {
        throw std::runtime_error("Time out of range in Cheb3D");
    }

    double tau = (2.0*t - Ta - Tb)/(Tb - Ta);

    Matrix result(3, 1);

    result(1,1) = 0.0;
    for (int i = 1; i <= N; ++i) {
        double T;
        if (i == 1) T = 1.0;
        else if (i == 2) T = tau;
        else T = 2.0 * tau * ChebT(tau, i-2) - ChebT(tau, i-3);

        result(1,1) += Cx(i,1) * T;
    }

    result(2,1) = 0.0;
    for (int i = 1; i <= N; ++i) {
        double T;
        if (i == 1) T = 1.0;
        else if (i == 2) T = tau;
        else T = 2.0 * tau * ChebT(tau, i-2) - ChebT(tau, i-3);

        result(2,1) += Cy(i,1) * T;
    }

    result(3,1) = 0.0;
    for (int i = 1; i <= N; ++i) {
        double T;
        if (i == 1) T = 1.0;
        else if (i == 2) T = tau;
        else T = 2.0 * tau * ChebT(tau, i-2) - ChebT(tau, i-3);

        result(3,1) += Cz(i,1) * T;
    }

    return result;
}

// Función auxiliar para calcular Tn(x)
double ChebT(double x, int n) {
    if (n == 0) return 1.0;
    if (n == 1) return x;
    return 2.0 * x * ChebT(x, n-1) - ChebT(x, n-2);
}