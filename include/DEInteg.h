//
// Created by adria on 12/05/2025.
//

#ifndef PROYECTOTALLER_DEINTEG_H
#define PROYECTOTALLER_DEINTEG_H


#include "Matrix.h"
#include "Sat_const.h"
#include <vector>
#include <functional>

enum DEState {
    DE_INIT    = 1,
    DE_DONE    = 2,
    DE_BADACC  = 3,
    DE_NUMSTEPS= 4,
    DE_STIFF   = 5,
    DE_INVPARAM= 6
};

Matrix DEInteg(
        void (*func)(double t, const Matrix& y, Matrix& dydt),
        double t,
        double tout,
        double relerr,
        double abserr,
        int n_eqn,
        const Matrix& y
);

#endif //PROYECTOTALLER_DEINTEG_H
