//
// Created by advecino on 10/04/2025.
//

#include "../include/Legendre.h"
#include <cmath>
// fi [rad]
void Legendre(int n, int m, double fi, Matrix& pnm, Matrix& dpnm) {

    pnm = Matrix(n+1, m+1);
    dpnm = Matrix(n+1, m+1);


    pnm(0, 0) = 1.0;
    dpnm(0, 0) = 0.0;
    //diagonal coefficients
    if (n >= 1 && m >= 1) {
        pnm(1, 1) = sqrt(3.0) * cos(fi);
        dpnm(1, 1) = -sqrt(3.0) * sin(fi);
    }

    for (int i = 2; i <= n; i++) {
        double factor = sqrt((2.0*i + 1.0)/(2.0*i));
        pnm(i, i) = factor * cos(fi) * pnm(i-1, i-1);
    }

    for (int i = 2; i <= n; i++) {
        double factor = sqrt((2.0*i + 1.0)/(2.0*i));
        dpnm(i, i) = factor * (cos(fi) * dpnm(i-1, i-1) - sin(fi) * pnm(i-1, i-1));
    }


    for (int i = 1; i <= n; i++) {
        pnm(i, i-1) = sqrt(2.0*i + 1.0) * sin(fi) * pnm(i-1, i-1);
    }

    for (int i = 1; i <= n; i++) {
        dpnm(i, i-1) = sqrt(2.0*i + 1.0) * (cos(fi) * pnm(i-1, i-1) + sin(fi) * dpnm(i-1, i-1));
    }


    int j = 0;
    int k = 2;
    while (true) {
        for (int i = k; i <= n; i++) {
            double term1 = sqrt((2.0*i + 1.0)/((i-j)*(i+j)));
            double term2 = sqrt(2.0*i - 1.0) * sin(fi) * pnm(i-1, j);
            double term3 = sqrt(((i+j-1.0)*(i-j-1.0))/(2.0*i-3.0)) * pnm(i-2, j);
            pnm(i, j) = term1 * (term2 - term3);
        }
        j++;
        k++;
        if (j > m) break;
    }

    j = 0;
    k = 2;
    while (true) {
        for (int i = k; i <= n; i++) {
            double term1 = sqrt((2.0*i + 1.0)/((i-j)*(i+j)));
            double term2 = sqrt(2.0*i - 1.0) * sin(fi) * dpnm(i-1, j);
            double term3 = sqrt(2.0*i - 1.0) * cos(fi) * pnm(i-1, j);
            double term4 = sqrt(((i+j-1.0)*(i-j-1.0))/(2.0*i-3.0)) * dpnm(i-2, j);
            dpnm(i, j) = term1 * (term2 + term3 - term4);
        }
        j++;
        k++;
        if (j > m) break;
    }
}