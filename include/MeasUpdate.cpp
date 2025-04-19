#include "Matrix.h"
void MeasUpdate(Matrix& x,  Matrix& z,  Matrix& g,  Matrix& s,  Matrix& G, Matrix& P, int n, Matrix& K_out) {
    int m = z.getFilas();  // Longitud del vector de mediciones

    // Crear matriz diagonal Inv_W (covarianza de las medidas invertida)
    Matrix Inv_W(m, m);
    for (int i = 1; i <= m; ++i) {
        double val = s(i, 1);  // <- CORREGIDO para base-1
        Inv_W(i, i) = val * val;
    }

    // Kalman Gain: K = P * G^T * inv(Inv_W + G * P * G^T)
    Matrix Gt = G.transpuesta();
    Matrix GP = G * P;
    Matrix GPGt = GP * Gt;
    Matrix S = Inv_W + GPGt;
    Matrix S_inv = S.inversa();
    Matrix K = P * Gt * S_inv;

    // State update: x = x + K * (z - g)
    Matrix innovation = z - g;
    x = x + K * innovation;

    // Covariance update: P = (I - K * G) * P
    Matrix I(n, n);
    for (int i = 1; i <= n; ++i)  // <- CORREGIDO para base-1
        I(i, i) = 1.0;

    Matrix KG = K * G;
    Matrix I_KG = I - KG;
    P = I_KG * P;

    // Devolver K si se necesita fuera
    K_out = K;
}

