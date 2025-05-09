#include <string>
#include <stdexcept>
#include "Matrix.h"
void MeasUpdate(Matrix& x,  Matrix z,  Matrix g,  Matrix s,  Matrix G, Matrix& P, int n, Matrix& K) {
    int m = z.getFilas();

    // Verificación crítica de dimensiones
    if (s.getFilas() != m || s.getColumnas() != 1) {
        throw std::invalid_argument(
                "Dimensiones de s (" + std::to_string(s.getFilas()) + "x" + std::to_string(s.getColumnas()) +
                ") no coinciden con z (" + std::to_string(m) + "x1)");
    }

    Matrix Inv_W(m, m);
    for (int i = 1; i <= m; ++i) {
        double val = s(i, 1);
        Inv_W(i, i) = val * val; //Inverse weight (measurement covariance)
    }

    // Kalman gain
    Matrix Gt = G.transpuesta();
    Matrix S = Inv_W + (G * P * Gt);
    Matrix S_inv = S.inversa();
    K = P * Gt * S_inv;

    // State update
    Matrix innovation = z - g;
    x = x + (K * innovation);


    // Covariance update
    Matrix I(n, n);
    for (int i = 1; i <= n; ++i) {
        I(i, i) = 1.0;
    }
    Matrix KG = K * G;
    Matrix I_KG = I - KG;
    P = I_KG * P * I_KG.transpuesta() + (K * Inv_W * K.transpuesta());
}