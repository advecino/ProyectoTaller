#include "../include/Legendre.h"
#include <cmath>


void Legendre(int n, int m, double fi, Matrix& pnm, Matrix& dpnm) {
    pnm = Matrix(n+2, m+2);
    dpnm = Matrix(n+2, m+2);


    pnm(1,1) = 1.0;
    dpnm(1,1) = 0.0;
    if (n >= 1 && m >= 1) {
        pnm(2,2) = std::sqrt(3.0)*std::cos(fi);
        dpnm(2,2) = -std::sqrt(3.0)*std::sin(fi);
    }


    for (int i = 2; i <= n; ++i) {
        double factor = std::sqrt((2.0*i+1)/(2.0*i));
        pnm(i+1,i+1) = factor * std::cos(fi) * pnm(i,i);
        dpnm(i+1,i+1) = factor * (std::cos(fi)*dpnm(i,i) - std::sin(fi)*pnm(i,i));
    }


    for (int i = 1; i <= n; ++i) {
        double f = std::sqrt(2.0*i+1);
        pnm(i+1,i)   = f * std::sin(fi) * pnm(i,i);
        dpnm(i+1,i) = f * (std::cos(fi)*pnm(i,i) + std::sin(fi)*dpnm(i,i));
    }


    for (int j = 0; j <= m; ++j) {
        for (int i = j+2; i <= n; ++i) {
            double denom = (i-j)*(i+j);
            if (denom <= 0) continue;
            double a = std::sqrt((2.0*i+1)/denom);
            double b = std::sqrt(2.0*i-1);
            double c = ((i+j-1)*(i-j-1) > 0)
                       ? std::sqrt(((i+j-1)*(i-j-1))/(2.0*i-3))
                       : 0.0;
            pnm(i+1,j+1) = a*(b*std::sin(fi)*pnm(i,j+1) - c*pnm(i-1,j+1));
            dpnm(i+1,j+1)= a*(b*std::sin(fi)*dpnm(i,j+1)
                              + b*std::cos(fi)*pnm(i,j+1)
                              - c*dpnm(i-1,j+1));
        }
    }
}