#include "../include/unit.h"
/*
%--------------------------------------------------------------------------
%
%  unit.m
%
%  this function calculates a unit vector given the original vector. if a
%  zero vector is input, the vector is set to zero.
%
%  input:
%    vec         - vector
%
%  output:
%    outvec      - unit vector
%
%--------------------------------------------------------------------------*/

Matrix unit(const Matrix& vec) {
    const double small = 0.000001;
    Matrix outvec(3, 1);


    double magv = 0.0;
    for (int i = 1; i <= 3; ++i) {
        magv += vec(i, 1) * vec(i, 1);
    }
    magv = sqrt(magv);

    if (magv > small) {
        for (int i = 1; i <= 3; ++i) {
            outvec(i, 1) = vec(i, 1) / magv;
        }
    } else {
        for (int i = 1; i <= 3; ++i) {
            outvec(i, 1) = 0.0;
        }
    }

    return outvec;
}