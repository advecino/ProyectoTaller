#include "../include/TimeUpdate.h"

void TimeUpdate(Matrix& P, const Matrix& Phi, const Matrix* Qdt)
{
    const int n = 3;

    Matrix temp(n, n);

    for(int i = 1; i <= n; ++i) {
        for(int j = 1; j <= n; ++j) {
            double sum = 0.0;
            for(int k = 1; k <= n; ++k) {
                sum += Phi(i,k) * P(k,j);
            }
            temp(i,j) = sum;
        }
    }


    Matrix result(n, n);
    for(int i = 1; i <= n; ++i) {
        for(int j = 1; j <= n; ++j) {
            double sum = 0.0;
            for(int k = 1; k <= n; ++k) {
                sum += temp(i,k) * Phi(j,k);
            }
            result(i,j) = sum;
        }
    }


    if(Qdt != nullptr) {
        for(int i = 1; i <= n; ++i) {
            for(int j = 1; j <= n; ++j) {
                result(i,j) += (*Qdt)(i,j);
            }
        }
    }

    P = result;
}