//
// Created by adria on 12/05/2025.
//

#ifndef PROYECTOTALLER_DEINTEG_H
#define PROYECTOTALLER_DEINTEG_H


#include "Matrix.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <functional>

class DEInteg {
private:
    // Constantes de estado
    enum DE_STATE {
        DE_INIT = 1,    // Restart integration
        DE_DONE = 2,    // Successful step
        DE_BADACC = 3,  // Accuracy requirement could not be achieved
        DE_NUMSTEPS = 4,// Permitted number of steps exceeded
        DE_STIFF = 5,   // Stiff problem suspected
        DE_INVPARAM = 6 // Invalid input parameters
    };

    // Variables miembro
    DE_STATE State_;
    bool PermitTOUT;
    double told;
    bool OldPermit;
    double delsgn;
    bool start;
    double x;
    Matrix yy;
    double h;
    double hold;
    double hnew;
    int k;
    int kold;
    bool phase1;
    bool nornd;
    bool stiff;
    int nostep;
    int kle4;
    int ifail;
    bool crash;

    // Vectores de trabajo
    Matrix wt;
    Matrix p;
    Matrix yp;
    Matrix phi;
    std::vector<double> g;
    std::vector<double> sig;
    std::vector<double> rho;
    std::vector<double> w;
    std::vector<double> alpha;
    std::vector<double> beta;
    std::vector<double> v;
    std::vector<double> psi_;

    // Constantes
    static const std::vector<double> two;
    static const std::vector<double> gstr;

public:
    DEInteg();
    Matrix integrate(std::function<Matrix(double, const Matrix&)> func,
                     double t, double tout, double relerr, double abserr,
                     Matrix& y);
};
#endif //PROYECTOTALLER_DEINTEG_H
