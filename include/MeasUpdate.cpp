#include "Matrix.h"
void MeasUpdate(Matrix& x,  Matrix& z,  Matrix& g,  Matrix& s,  Matrix& G, Matrix& P, int n, Matrix& K_out) {
    int m = z.getFilas();


    Matrix Inv_W(m, m);
    for (int i = 1; i <= m; ++i) {
        double val = s(i, 1);
        Inv_W(i, i) = val * val;
    }

    // Kalman gain
    Matrix Gt = G.transpuesta();
    Matrix GP = G * P;
    Matrix GPGt = GP * Gt;
    Matrix S = Inv_W + GPGt;
    Matrix S_inv = S.inversa();
    Matrix K = P * Gt * S_inv;

    // State update
    Matrix innovation = z - g;
    x = x + K * innovation;

    // Covariance update
    Matrix I(n, n);
    for (int i = 1; i <= n; ++i)
        I(i, i) = 1.0;

    Matrix KG = K * G;
    Matrix I_KG = I - KG;
    P = I_KG * P;

    K_out = K;
}