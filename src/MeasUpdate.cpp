#include <string>
#include <stdexcept>
#include "../include/Matrix.h"

void MeasUpdate( Matrix& x, Matrix& P, Matrix& K, Matrix& z, const Matrix& g, const Matrix& s,const Matrix& G,int n) {
    int m = z.getFilas();
    if (z.getColumnas() != 1 || g.getFilas() != m || g.getColumnas() != 1
        || s.getFilas() != m || s.getColumnas() != 1
        || G.getFilas() != m || G.getColumnas() != n)
    {
        throw std::invalid_argument("Dimension mismatch in MeasUpdate inputs");
    }

    // Build inverse weight (measurement covariance)
    Matrix Inv_W(m, m);
    for (int i = 1; i <= m; ++i) {
        double si = s(i,1);
        Inv_W(i, i) = si * si;
    }

    Matrix S = Inv_W + G * P * G.transpuesta();

    Matrix S_inv = S.inversa();

    // Kalman gain
    K = P * G.transpuesta() * S_inv;

    // State update
    x = x + K * (z - g);

    // Covariance update
    Matrix I(n, n);
    for (int i = 1; i <= n; ++i) {
        I(i, i) = 1.0;
    }
    P = (I - K * G) * P;
}