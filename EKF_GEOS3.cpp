#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <cstring>
#include "include/Sat_const.h"
#include "include/Mjday.h"
#include "include/Position.h"
#include "include/anglesg.h"
#include "include/Matrix.h"
#include "include/DEInteg.h"
#include "include/TimeDiff.h"
#include "include/IERS.h"
#include "include/R_z.h"
#include "include/LTC.h"
#include "include/AuxParam.h"
#include "include/AzElPa.h"
#include "include/TimeUpdate.h"
#include "include/MeasUpdate.h"
#include "include/Accel.h"
#include "include/global.h"

using namespace std;

int main() {
    // Constantes y estructuras
    Matrix Cnm(181,181), Snm(181,181);

    // Leer GGM03S.txt
    ifstream fin("./data/GGM03S.txt");
    if (!fin) {
        cerr << "No se pudo abrir GGM03S.txt" << endl;
        return 1;
    }
    for (int n = 0; n <= 180; ++n) {
        for (int m = 0; m <= n; ++m) {
            int dummy1, dummy2;
            double cnm, snm, dummy3, dummy4;
            fin >> dummy1 >> dummy2 >> cnm >> snm >> dummy3 >> dummy4;
            Cnm(n+1, m+1) = cnm;
            Snm(n+1, m+1) = snm;
        }
    }
    fin.close();

    // Leer eop19620101.txt
    ifstream eopfile("./data/eop19620101.txt");
    if (!eopfile) {
        cerr << "No se pudo abrir eop19620101.txt" << endl;
        return 1;
    }
    vector<vector<double>> eop_vec(13);
    double value;
    while (eopfile >> value) {
        static int col = 0;
        eop_vec[col % 13].push_back(value);
        ++col;
    }
    int rows = eop_vec[0].size();
    eopdata = Matrix(13, rows);
    for (int i = 1; i <= 13; ++i)
        for (int j = 1; j <= rows; ++j)
            eopdata(i,j) = eop_vec[i-1][j-1];

    // Leer DE430
    ifstream pcfile("./data/DE430Coeff.txt");
    if (!pcfile) {
        cerr << "No se pudo abrir DE430Coeff.txt" << endl;
        return 1;
    }
    Matrix PC(2285, 1020);
    for (int i = 1; i <= 2285; ++i)
        for (int j = 1; j <= 1020; ++j)
            pcfile >> PC(i,j);

    // Leer observaciones
    ifstream obsfile("./data/GEOS3.txt");
    if (!obsfile) {
        cerr << "No se pudo abrir GEOS3.txt" << endl;
        return 1;
    }
    int nobs = 46;
    Matrix obs(nobs, 4);
    for (int i = 1; i <= nobs; ++i) {
        string line;
        getline(obsfile, line);
        int Y = stoi(line.substr(0,4));
        int M = stoi(line.substr(5,2));
        int D = stoi(line.substr(8,2));
        int h = stoi(line.substr(12,2));
        int m = stoi(line.substr(15,2));
        double s = stod(line.substr(18,6));
        double az = stod(line.substr(25,8));
        double el = stod(line.substr(35,8));
        double dist = stod(line.substr(44,9));

        obs(i,1) = Mjday(Y,M,D,h,m,s);
        obs(i,2) = az * Rad;
        obs(i,3) = el * Rad;
        obs(i,4) = dist * 1e3;
    }

    // Parámetros
    double sigma_range = 92.5;
    double sigma_az = 0.0224 * Rad;
    double sigma_el = 0.0139 * Rad;

    Matrix Rs = Position(-158.2706 * Rad, 21.5748 * Rad, 300.20, R_Earth, f_Earth);

    double Mjd1 = obs(1,1);
    double Mjd2 = obs(9,1);
    double Mjd3 = obs(18,1);

    AuxParam aux;
    aux.Mjd_UTC = Mjd2;
    aux.n = 20;
    aux.m = 20;
    aux.sun = 1;
    aux.moon = 1;
    aux.planets = 1;

    AnglesGResult ag = anglesg(obs(1,2), obs(9,2), obs(18,2),
                               obs(1,3), obs(9,3), obs(18,3),
                               Mjd1, Mjd2, Mjd3,
                               Rs, Rs, Rs,
                               aux, eopdata);

    Matrix Y0_apr = Matrix::concatenate(ag.r2,ag.v2, 0);

    double Mjd0 = Mjday(1995, 1, 29, 2, 38, 0);
    aux.Mjd_UTC = obs(9,1);
    aux.Mjd_TT = aux.Mjd_UTC; // temporal

    DEInteg integrator;
    Matrix Y = integrator.integrate(Accel, 0.0, -(obs(9,1)-Mjd0)*86400.0, 1e-13, 1e-6, Y0_apr);

    Matrix P = Matrix::identity(6);
    for (int i = 1; i <= 3; ++i) P(i,i) = 1e8;
    for (int i = 4; i <= 6; ++i) P(i,i) = 1e3;

    Matrix LT = LTC(-158.2706*Rad, 21.5748*Rad);

    Matrix yPhi(42,1);
    Matrix Phi(6,6);

    double t = 0.0;
    for (int i = 1; i <= nobs; ++i) {
        double t_old = t;
        Matrix Y_old = Y;
        aux.Mjd_UTC = obs(i,1);
        t = (obs(i,1)-Mjd0)*86400.0;

        IERSResult eop = IERS(eopdata, aux.Mjd_UTC, 'l');
        TimeDiffs dT = timediff(eop.UT1_UTC, eop.TAI_UTC);
        double Mjd_TT = aux.Mjd_UTC + dT.TT_UTC/86400.0;
        double Mjd_UT1 = Mjd_TT + (eop.UT1_UTC - dT.TT_UTC)/86400.0;
        aux.Mjd_TT = Mjd_TT;

        for (int ii = 1; ii <= 6; ++ii) {
            yPhi(ii,1) = Y_old(ii,1);
            for (int jj = 1; jj <= 6; ++jj)
                yPhi(6*jj + ii,1) = (ii==jj ? 1.0 : 0.0);
        }

        yPhi = integrator.integrate(VarEqn, 0.0, t - t_old, 1e-13, 1e-6, yPhi);

        for (int j = 1; j <= 6; ++j)
            for (int k = 1; k <= 6; ++k)
                Phi(k,j) = yPhi(6*j + k,1);

        Y = integrator.integrate(Accel, 0.0, t - t_old, 1e-13, 1e-6, Y_old);

        double theta = gmst(Mjd_UT1);
        Matrix U = R_z(theta);
        Matrix r = Y.getSubMatrix(1,3,1,1);
        Matrix s = LT * (U * r - Rs);

        P = TimeUpdate(P, Phi);

        AzElPaResult azel = AzElPa(s);
        Matrix dAdY = azel.dAds * LT * U;
        dAdY = dAdY.concatenate(Matrix(1,3), 2);

        MeasUpdate(Y, P, Matrix(), Matrix(1,1,obs(i,2)), Matrix(1,1,azel.Az), Matrix(1,1,sigma_az), dAdY, 6);

        r = Y.getSubMatrix(1,3,1,1);
        s = LT * (U * r - Rs);
        azel = AzElPa(s);
        Matrix dEdY = azel.dEds * LT * U;
        dEdY = dEdY.concatenate(Matrix(1,3), 2);
        MeasUpdate(Y, P, Matrix(), Matrix(1,1,obs(i,3)), Matrix(1,1,azel.El), Matrix(1,1,sigma_el), dEdY, 6);

        r = Y.getSubMatrix(1,3,1,1);
        s = LT * (U * r - Rs);
        double Dist = s.norm();
        Matrix dDdY = (s / Dist).transpuesta() * LT * U;
        dDdY = dDdY.concatenate(Matrix(1,3), 2);
        MeasUpdate(Y, P, Matrix(), Matrix(1,1,obs(i,4)), Matrix(1,1,Dist), Matrix(1,1,sigma_range), dDdY, 6);
    }

    Matrix Y0 = integrator.integrate(Accel, 0.0, -(obs(46,1) - obs(1,1))*86400.0, 1e-13, 1e-6, Y);
    Matrix Y_true(6,1);
    Y_true(1,1) = 5753.173e3;
    Y_true(2,1) = 2673.361e3;
    Y_true(3,1) = 3440.304e3;
    Y_true(4,1) = 4.324207e3;
    Y_true(5,1) = -1.924299e3;
    Y_true(6,1) = -5.728216e3;

    std::cout << "\nError de estimación de posición [m]:\n";
    std::cout << "dX = " << Y0(1,1) - Y_true(1,1) << endl;
    std::cout << "dY = " << Y0(2,1) - Y_true(2,1) << endl;
    std::cout << "dZ = " << Y0(3,1) - Y_true(3,1) << endl;

    std::cout << "\nError de estimación de velocidad [m/s]:\n";
    std::cout << "dVx = " << Y0(4,1) - Y_true(4,1) << endl;
    std::cout << "dVy = " << Y0(5,1) - Y_true(5,1) << endl;
    std::cout << "dVz = " << Y0(6,1) - Y_true(6,1) << endl;

    return 0;
}
