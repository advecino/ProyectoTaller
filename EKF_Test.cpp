/**
 * @brief Fichero para probar todos los métodos con test unitarios
 * @return Devolverá el numero de test ejecutados y si estan correctos mostrará PASSED
 */

#include <iostream>
#include <cstdio>
#include <cmath>
#include "./include/Mjday.h"
#include "./include/Matrix.h"
#include "./include/R_x.h"
#include "./include/TimeUpdate.h"
#include "./include/Position.h"
#include "./include/sign_.h"
#include "include/MeasUpdate.h"
#include "include/AccelPointMass.h"
#include "include/Mjday_TDB.h"
#include "include/angl.h"
#include "include/Cheb3D.h"
#include "include/MeanObliquity.h"
#include "include/NutAngles.h"
#include "include/AccelHarmonic.h"
#include "include/G_AccelHarmonic.h"
#include "include/EqnEquinox.h"
#include "include/gast.h"
#include "include/unit.h"
#include "include/gibbs.h"
#include "include/hgibbs.h"
#include "include/elements.h"
#include "include/LTC.h"
#include "include/GHAMatrix.h"
#include "include/PoleMatrix.h"
#include "include/NutMatrix.h"
#include "include/PrecMatrix.h"
#include "include/timediff.h"
#include "include/Geodetic.h"
#include "include/doubler.h"
#include "include/IERS.h"
#include "include/JPL_Eph_DE430.h"
#include "include/AzElPa.h"
#include "include/VarEqn.h"
#include "include/anglesdr.h"
#include "include/anglesg.h"
#include "include/Legendre.h"
#include "include/global.h"
#include "include/Accel.h"
#include "include/DEInteg.h"
#include "src/cargarPCdesdeFichero.cpp"

#define TOL_ 10e-14

int tests_run = 0;

#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)

using namespace std;


bool MatrixEqual(const Matrix& a, const Matrix& b, double tol = TOL_) {
    for(int i = 1; i <= 3; ++i) {
        for(int j = 1; j <= 3; ++j) {
            if(fabs(a(i,j) - b(i,j)) > tol) {
                return false;
            }
        }
    }
    return true;
}

int Matrix_Basico() {
    try {

        Matrix m(2, 3);

        _assert(m.getFilas() == 2);
        _assert(m.getColumnas() == 3);


        for (int i = 1; i <= 2; ++i) {
            for (int j = 1; j <= 3; ++j) {
                _assert(m(i,j) == 0.0);
            }
        }


        m(1,1) = 1.5;
        m(2,3) = 2.0;
        _assert(m(1,1) == 1.5);
        _assert(m(2,3) == 2.0);

        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error en Test_Matrix_Basico: " << e.what() << std::endl;
        return 1;
    }
}

int Position_01() {

    double lon = 0.0;
    double lat = 0.0;
    double h   = 0.0;
    Matrix r = Position(lon, lat, h, R_Earth, f_Earth);


    double TOL = 1e-6;
    _assert(std::fabs(r(1,1) - R_Earth) < TOL);
    _assert(std::fabs(r(2,1))          < TOL);
    _assert(std::fabs(r(3,1))          < TOL);

    return 0;
}

int Position_02() {

    double lon = 1.234;
    double lat = M_PI/2.0;                   // 90°
    double h   = 1000.0;


    const double e2 = f_Earth * (2.0 - f_Earth);

    double N = R_Earth / std::sqrt(1.0 - e2);

    double expected_z = (1.0 - e2) * N + h;

    Matrix r = Position(lon, lat, h, R_Earth, f_Earth);


    constexpr double TOL_POS = 1e-6;
    constexpr double TOL_Z   = 1e-3; // 1 mm


    _assert(std::fabs(r(1,1))             < TOL_POS);
    _assert(std::fabs(r(2,1))             < TOL_POS);
    _assert(std::fabs(r(3,1) - expected_z) < TOL_Z);

    return 0;
}

int TimeUpdate_01() {
    double phi_data[] = {1,0,0,0,1,0,0,0,1};
    double p_data[] = {2,0,0,0,3,0,0,0,4};

    Matrix Phi(3, 3, phi_data, 9);
    Matrix P(3, 3, p_data, 9);
    Matrix expected(3, 3, p_data, 9);

    TimeUpdate(P, Phi);

    if(!MatrixEqual(P, expected)) {
        std::cout << "TimeUpdate_01 failed:\nExpected:\n";
        expected.print();
        std::cout << "Got:\n";
        P.print();
        FAIL();
    }
    return 0;
}

int Mjday_01()
{
    _assert(fabs(Mjday(2025,4,3,15,37,5)-60768.6507523148) < pow(10,-10));

    return 0;
}

int Mjday_02()
{
    _assert(fabs(Mjday(2025,4,3,0,0,0.0)-Mjday(2025,4,3)) < pow(10,-10));

    return 0;
}

int R_x_01()
{
    double alpha = 1.0;
    Matrix sol(3, 3);

    sol = R_x(alpha);

    _assert(fabs(sol(1,1)) - 1 < TOL_ && fabs(sol(1,2)) < TOL_ && fabs(sol(1,3)) < TOL_);
    _assert(sol(2,1) < TOL_ && fabs(sol(2,2) - 0.54030230586814 ) < TOL_ && fabs(sol(2,3) -0.841470984807897) < TOL_);


    return 0;
}

int R_y_01()
{
    double alpha = 1.0;
    Matrix sol(3, 3);

    sol = R_y(alpha);

    _assert(fabs(sol(1,1) - 0.54030230586814) < TOL_);
    _assert(fabs(sol(1,2)) < TOL_);
    _assert(fabs(sol(1,3) + 0.841470984807897) < TOL_);

    _assert(fabs(sol(2,1)) < TOL_);
    _assert(fabs(sol(2,2) - 1.0) < TOL_);
    _assert(fabs(sol(2,3)) < TOL_);

    _assert(fabs(sol(3,1) - 0.841470984807897) < TOL_);
    _assert(fabs(sol(3,2)) < TOL_);
    _assert(fabs(sol(3,3) - 0.54030230586814) < TOL_);

    return 0;
}

int R_z_01()
{
    double alpha = 1.0;
    Matrix sol(3, 3);

    sol = R_z(alpha);

    _assert(fabs(sol(1,1) - 0.54030230586814) < TOL_);
    _assert(fabs(sol(1,2) - 0.841470984807897) < TOL_);
    _assert(fabs(sol(1,3)) < TOL_);

    _assert(fabs(sol(2,1) + 0.841470984807897) < TOL_);
    _assert(fabs(sol(2,2) - 0.54030230586814) < TOL_);
    _assert(fabs(sol(2,3)) < TOL_);

    _assert(fabs(sol(3,1)) < TOL_);
    _assert(fabs(sol(3,2)) < TOL_);
    _assert(fabs(sol(3,3) - 1.0) < TOL_);

    return 0;
}

int Legendre_01()
{
    int n = 2, m = 2;
    double fi = M_PI/6; // 30°
    Matrix pnm(0,0), dpnm(0,0);
    Legendre(n, m, fi, pnm, dpnm);

    _assert(fabs(pnm(1,1) - 1.0) < TOL_);
    _assert(fabs(pnm(2,2) - std::sqrt(3.0)*std::cos(fi)) < TOL_);
    _assert(fabs(dpnm(2,2) + std::sqrt(3.0)*std::sin(fi)) < TOL_);

    return 0;
}

int sign_() {
    // Casos positivos
    _assert(sign_(5.0, 3.0) == 5.0);
    _assert(sign_(-5.0, 3.0) == 5.0);

    // Casos negativos
    _assert(sign_(5.0, -3.0) == -5.0);
    _assert(sign_(-5.0, -3.0) == -5.0);

    // Casos con cero
    _assert(sign_(5.0, 0.0) == 5.0);
    _assert(sign_(-5.0, 0.0) == 5.0);

    return 0;
}

int AccelPointMass_01() {
    Matrix r(3, 1, new double[3]{7000.0, 0.0, 0.0}, 3);
    Matrix s(3, 1, new double[3]{0.0, 0.0, 0.0}, 3);
    double GM = 3.986004418e14; // GM de la Tierra

    Matrix a = AccelPointMass(r, s, GM);

    double expected = -GM / (7000.0*7000.0*7000.0) * 7000.0;

    if(fabs(a(1,1) - expected) > TOL_ ||
       fabs(a(2,1)) > TOL_ ||
       fabs(a(3,1)) > TOL_) {
        return 1;
    }

    return 0;
}

int Mjday_TDB_01() {

    double Mjd_TT = 54930.5;

    double Mjd_TDB = Mjday_TDB(Mjd_TT);


    double expected = 54930.5000007235;

    if(fabs(Mjd_TDB - expected) > 1e-5) {
        return 1;
    }

    return 0;
}

int Mjday_TDB_02() {

    double Mjd_TT = 51544.5;
    double Mjd_TDB = Mjday_TDB(Mjd_TT);
    double delta = Mjd_TDB - Mjd_TT;

    _assert(std::fabs(delta) < 1e-4);

    return 0;
}

int Mjday_TDB_03() {

    double Mjd0 = 58000.0;
    double Mjd1 = Mjd0 + 10.5;
    double TDB0 = Mjday_TDB(Mjd0);
    double TDB1 = Mjday_TDB(Mjd1);
    double dt_in  = Mjd1 - Mjd0;
    double dt_out = TDB1 - TDB0;

    _assert(std::fabs(dt_out - dt_in) < 1e-4);

    return 0;
}

int angl_01() {

    Matrix vec1(3, 1, new double[3]{1.0, 0.0, 0.0}, 3);
    Matrix vec2(3, 1, new double[3]{2.0, 0.0, 0.0}, 3);

    double angle = angl(vec1, vec2);

    if(fabs(angle - 0.0) > TOL_) {
        std::cout << "angl_01 failed! Angle: " << angle << std::endl;
        return 1;
    }

    return 0;
}

int angl_02() {

    Matrix vec1(3, 1, new double[3]{1.0, 0.0, 0.0}, 3);
    Matrix vec2(3, 1, new double[3]{0.0, 1.0, 0.0}, 3);

    double angle = angl(vec1, vec2);
    double expected = M_PI/2.0;

    if(fabs(angle - expected) > TOL_) {
        std::cout << "angl_02 failed! Angle: " << angle << " Expected: " << expected << std::endl;
        return 1;
    }

    return 0;
}

int Cheb3D_01() {
    try {

        const int N = 2;
        double Ta = 0.0;
        double Tb = 1.0;



        Matrix Cx(N,1, new double[N]{0.5, 0.5}, N);

        Matrix Cy(3,1, new double[3]{0.5, 0.5, 0.25}, 3);

        double t = 0.5;
        Matrix result = Cheb3D(t, 2, Ta, Tb, Cx, Cx, Cx);


        double expected = 0.5;
        const double TOL = 1e-10;

        if (fabs(result(1,1) - expected) > TOL ||
            fabs(result(2,1) - expected) > TOL ||
            fabs(result(3,1) - expected) > TOL) {
            std::cerr << "Cheb3D_01 failed!\n";
            return 1;
        }

        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Exception: " << e.what() << "\n";
        return 1;
    }
}

int Cheb3D_02() {
    try {

        Matrix Cx(2,1, new double[2]{0.5, 0.5}, 2);

        Matrix result = Cheb3D(1.0, 2, 0.0, 1.0, Cx, Cx, Cx);

        const double TOL = 1e-10;
        double expected = 1.0;

        if(fabs(result(1,1)-expected) > TOL ||
           fabs(result(2,1)-expected) > TOL ||
           fabs(result(3,1)-expected) > TOL) {
            std::cerr << "Cheb3D_04 failed at t=1\n";
            std::cerr << "Expected: " << expected << "\n";
            std::cerr << "Got: [" << result(1,1) << ", " << result(2,1) << ", " << result(3,1) << "]\n";
            return 1;
        }
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Exception in Cheb3D_04: " << e.what() << "\n";
        return 1;
    }
}

int MeanObliquity_01() {


    double Mjd_TT = MJD_J2000;


    double MOblq = MeanObliquity(Mjd_TT);

    double expected = 23.4392911 * Rad;

    if(fabs(MOblq - expected) > 1e-4) {
        std::cout << "MeanObliquity_01 failed!\n";
        std::cout << "Expected: " << expected << "\n";
        std::cout << "Got:      " << MOblq << "\n";
        return 1;
    }

    return 0;
}

int EqnEquinox_01() {
    try {

        double Mjd_TT = 51544.5;

        double EqE = EqnEquinox(Mjd_TT);

        _assert(fabs(EqE) < 1e-3);
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Exception in EqnEquinox_01: " << e.what() << std::endl;
        return 1;
    }
}

int EqnEquinox_02() {
    try {
        double Mjd_TT = 58849.0;

        double EqE = EqnEquinox(Mjd_TT);

        double expected = -0.0033534 * M_PI/180.0/3600.0;

        _assert(fabs(EqE - expected) < 1e-3);
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Exception in EqnEquinox_02: " << e.what() << std::endl;
        return 1;
    }
}

int NutAngles_01() {
    try {
        double Mjd_TT = 51544.5;

        double dpsi, deps;
        NutAngles(Mjd_TT, dpsi, deps);

        if(fabs(dpsi) > 1e-3 || fabs(deps) > 1e-3) {
            return 1;
        }
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Exception in NutAngles_01: " << e.what() << "\n";
        return 1;
    }
}

int NutAngles_02() {
    double dpsi, deps;
    NutAngles(51544.5, dpsi, deps);




    _assert(std::fabs(dpsi) < 1e-2);
    _assert(std::fabs(deps) < 1e-2);

    return 0;
}

int NutAngles_03() {
    double dpsi0, deps0, dpsi1, deps1;
    NutAngles(51544.5,     dpsi0, deps0);
    NutAngles(51544.5+18262.5, dpsi1, deps1);


    _assert(std::fabs(dpsi1-dpsi0) > 1e-6 ||
            std::fabs(deps1-deps0) > 1e-6);


    return 0;
}

int AccelHarmonic_01() {

    for(int i=0;i<300;i++) for(int j=0;j<300;j++){
            Cnm[i][j]=0; Snm[i][j]=0;
        }
    Cnm[1][1]=1.0;

    double r_val=7000e3;
    Matrix r(3,1); r(1,1)=r_val; r(2,1)=0; r(3,1)=0;
    Matrix E(3,3);
    for(int i=1;i<=3;i++) E(i,i)=1.0;

    Matrix a=AccelHarmonic(r,E,0,0);
    double gm=398600.4415e9;
    double expected_ax=-gm/(r_val*r_val);
    _assert(fabs(a(1,1)-expected_ax)<TOL_);
    _assert(fabs(a(2,1))<TOL_);
    _assert(fabs(a(3,1))<TOL_);


    return 0;
}

int AccelHarmonic_02() {

    for(int i=0;i<300;i++) for(int j=0;j<300;j++){
            Cnm[i][j]=0; Snm[i][j]=0;
        }
    Cnm[1][1]=1.0;

    double r_val=8000e3;
    Matrix r(3,1); r(1,1)=0; r(2,1)=r_val; r(3,1)=0;
    Matrix E(3,3);
    for(int i=1;i<=3;i++) E(i,i)=1.0;

    Matrix a=AccelHarmonic(r,E,0,0);
    double gm=398600.4415e9;
    double expected_ay=-gm/(r_val*r_val);
    _assert(fabs(a(1,1))<TOL_);
    _assert(fabs(a(2,1)-expected_ay)<TOL_);
    _assert(fabs(a(3,1))<TOL_);

    return 0;
}

int AccelHarmonic_03() {
    Matrix I(3,3);
    for(int i=1;i<=3;++i) I(i,i)=1.0;
    Matrix r(3,1);
    const double R = 7000e3;
    r(1,1)=R; r(2,1)=0; r(3,1)=0;
    Matrix a = AccelHarmonic(r, I, 0, 0);
    double expect = -398600.4415e9/(R*R);
    _assert(fabs(a(1,1)-expect) < 1e-8);
    _assert(fabs(a(2,1)) < 1e-12);
    _assert(fabs(a(3,1)) < 1e-12);

    return 0;
}

int AccelHarmonic_04() {
    Matrix I(3,3);
    for(int i=1;i<=3;++i) I(i,i)=1.0;
    Matrix r(3,1);
    Matrix a = AccelHarmonic(r, I, 5, 5);
    _assert(a.norm() == 0.0);

    return 0;
}

int G_AccelHarmonic_01() {

    for(int i=0;i<300;i++) for(int j=0;j<300;j++){
            Cnm[i][j]=0.0; Snm[i][j]=0.0;
        }
    Cnm[1][1] = 1.0;

    double R = 7000e3;
    Matrix r(3,1); r(1,1)=R; r(2,1)=0; r(3,1)=0;
    Matrix E = Matrix(3,3);
    for(int i=1;i<=3;i++) E(i,i)=1.0;

    Matrix G = G_AccelHarmonic(r, E, 0, 0);

    double gm = 398600.4415e9;
    double r3 = R*R*R;
    double Gxx = 2*gm/r3;
    double Gyy = -gm/r3;
    double Gzz = Gyy;

    _assert(fabs(G(1,1) - Gxx) < 1e-6);
    _assert(fabs(G(2,2) - Gyy) < 1e-6);
    _assert(fabs(G(3,3) - Gzz) < 1e-6);

    for(int i=1;i<=3;i++) for(int j=1;j<=3;j++){
            if(i!=j) _assert(fabs(G(i,j)) < TOL_);
        }


    return 0;
}

int G_AccelHarmonic_02() {

    for(int i=0;i<300;i++) for(int j=0;j<300;j++){
            Cnm[i][j]=0.0; Snm[i][j]=0.0;
        }
    Cnm[1][1] = 1.0;

    double R = 8000e3;
    Matrix r(3,1); r(1,1)=0; r(2,1)=R; r(3,1)=0;
    Matrix E = Matrix(3,3);
    for(int i=1;i<=3;i++) E(i,i)=1.0;

    Matrix G = G_AccelHarmonic(r, E, 0, 0);

    double gm = 398600.4415e9;
    double r3 = R*R*R;
    double Gyy = 2*gm/r3;
    double Gxx = -gm/r3;
    double Gzz = Gxx;

    _assert(fabs(G(2,2) - Gyy) < 1e-6);
    _assert(fabs(G(1,1) - Gxx) < 1e-6);
    _assert(fabs(G(3,3) - Gzz) < 1e-6);

    for(int i=1;i<=3;i++) for(int j=1;j<=3;j++){
            if(i!=j) _assert(fabs(G(i,j)) < TOL_);
        }


    return 0;
}

int MeasUpdate_01() {
    int n = 2, m = 2;
    Matrix x(n,1), P(n,n), K(0,0);
    Matrix z(m,1), g(m,1), s(m,1), G(m,n);


    x(1,1) = 1.0; x(2,1) = -2.0;
    P(1,1) = 2.0; P(1,2) = 0.5;
    P(2,1) = 0.5; P(2,2) = 1.0;


    z(1,1) = 10.0; g(1,1) = 10.0;
    z(2,1) = -5.0; g(2,1) = -5.0;


    s(1,1) = 1.0;
    s(2,1) = 2.0;


    G(1,1) = 1.0; G(1,2) = 0.0;
    G(2,1) = 0.0; G(2,2) = 1.0;


    MeasUpdate(x, P, K, z, g, s, G, n);


    const double P11 = 0.66101694915254;
    const double P12 = 0.13559322033898;
    const double P21 = 0.13559322033898;
    const double P22 = 0.74576271186441;

    _assert(fabs(P(1,1) - P11) < TOL_);
    _assert(fabs(P(1,2) - P12) < TOL_);
    _assert(fabs(P(2,1) - P21) < TOL_);
    _assert(fabs(P(2,2) - P22) < TOL_);

    _assert(fabs(x(1,1) - 1.0) < TOL_);
    _assert(fabs(x(2,1) + 2.0) < TOL_);


    return 0;
}

int MeasUpdate_02() {

    int n = 1, m = 1;
    Matrix x(n,1), P(n,n), K(0,0);
    Matrix z(m,1), g(m,1), s(m,1), G(m,n);


    x(1,1) = 0.0;
    P(1,1) = 1.0;

    z(1,1) = 2.0;
    g(1,1) = 1.0;

    s(1,1) = 1.0;

    G(1,1) = 1.0;

    MeasUpdate(x, P, K, z, g, s, G, n);

    _assert(fabs(K(1,1) - 0.5) < TOL_);

    _assert(fabs(x(1,1) - 0.5) < TOL_);

    _assert(fabs(P(1,1) - 0.5) < TOL_);

    return 0;
}

int gstime_01() {
    try {
        double Mjd_UT1 = 51544.5;
        double expected = 4.894961;
        double result = gstime(Mjd_UT1);

        _assert(fabs(result - expected) < 1e-3);
        return 0;
    } catch(...) {
        return 1;
    }
}

int unit_01() {
    try {



        double vec_data[] = {1.0, 2.0, 3.0};
        Matrix vec(3, 1, vec_data, 3);
        Matrix outvec = unit(vec);

        double mag = sqrt(outvec(1,1)*outvec(1,1) + outvec(2,1)*outvec(2,1) + outvec(3,1)*outvec(3,1));


        _assert(fabs(mag - 1.0) < 1e-6);

        double expected[] = {0.267261, 0.534522, 0.801784};
        for (int i = 1; i <= 3; ++i) {
            _assert(fabs(outvec(i,1) - expected[i-1]) < 1e-6);
        }

        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error en unit_01 (Caso 1): " << e.what() << std::endl;
        return 1;
    }
}

int unit_02() {
    try {


        double vec_data[] = {1e-7, 1e-7, 1e-7};
        Matrix vec(3, 1, vec_data, 3);
        Matrix outvec = unit(vec);


        for (int i = 1; i <= 3; ++i) {
            _assert(outvec(i,1) == 0.0);
        }


        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error en unit_02 (Caso 2): " << e.what() << std::endl;
        return 1;
    }
}

int unit_03() {
    try {



        double vec_data[] = {0.577350, 0.577350, 0.577350}; // ~1/√3
        Matrix vec(3, 1, vec_data, 3);
        Matrix outvec = unit(vec);


        for (int i = 1; i <= 3; ++i) {
            _assert(fabs(outvec(i,1) - vec(i,1)) < 1e-6);
        }

        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error en unit_03 (Caso 3): " << e.what() << std::endl;
        return 1;
    }
}

int Gibbs_01() {


    Matrix r1(3, 1);
    r1(1,1) = 7000e3;
    r1(2,1) = 0;
    r1(3,1) = 0;

    Matrix r2(3, 1);
    r2(1,1) = 0;
    r2(2,1) = 7000e3;
    r2(3,1) = 0;

    Matrix r3(3, 1);
    r3(1,1) = -7000e3;
    r3(2,1) = 0;
    r3(3,1) = 0;


    GibbsResult result = gibbs(r1, r2, r3);
    double TOL = 1e-6;

    _assert(result.error == "          ok");
    _assert(fabs(result.theta - 1.570796) < TOL);
    _assert(fabs(result.theta1 - 1.570796) < TOL);
    _assert(fabs(result.copa) < TOL);


    _assert((std::abs(result.v2(1,1)) - (7546.05329011)) < TOL);
    _assert(fabs(result.v2(2,1)) < TOL);
    _assert(fabs(result.v2(3,1)) < TOL);

    return 0;
}

int HGibbs_01() {
    try {
        Matrix r1(3, 1), r2(3, 1), r3(3, 1);


        r1(1,1) = 7000e3; r1(2,1) = 0; r1(3,1) = 0;


        double angle = 0.5 * M_PI/180.0;
        r2(1,1) = 7000e3 * cos(angle);
        r2(2,1) = 7000e3 * sin(angle);
        r2(3,1) = 0;


        r3(1,1) = 7000e3 * cos(2*angle);
        r3(2,1) = 7000e3 * sin(2*angle);
        r3(3,1) = 0;


        double Mjd1 = 60000.0;
        double Mjd2 = Mjd1 + 1.0/(60.0*24.0);
        double Mjd3 = Mjd2 + 1.0/(60.0*24.0);


        GibbsResult result = hgibbs(r1, r2, r3, Mjd1, Mjd2, Mjd3);

        double TOL = 1e-6;


        _assert(result.error == "          ok");
        _assert(result.theta * 180.0/M_PI <= 1.0 + TOL);
        _assert(result.theta1 * 180.0/M_PI <= 1.0 + TOL);

        _assert(fabs(result.v2(1,1)) > 1.0);
        _assert(fabs(result.v2(2,1)) > 1.0);
        _assert(fabs(result.v2(3,1)) < TOL);

        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error en HGibbs_01: " << e.what() << std::endl;
        return 1;
    }
}

int HGibbs_02() {
    try {


        Matrix r1(3, 1), r2(3, 1), r3(3, 1);
        r1(1,1) = 7000e3; r1(2,1) = 0;    r1(3,1) = 0;
        r2(1,1) = 0;      r2(2,1) = 7000e3; r2(3,1) = 0;
        r3(1,1) = 0;      r3(2,1) = 0;     r3(3,1) = 7000e3;

        double Mjd1 = 60000.0;
        double Mjd2 = Mjd1 + 5.0/(60.0*24.0);
        double Mjd3 = Mjd2 + 5.0/(60.0*24.0);

        GibbsResult result = hgibbs(r1, r2, r3, Mjd1, Mjd2, Mjd3);

        _assert(result.error == "not coplanar");

        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error en test_HGibbs_02: " << e.what() << std::endl;
        return 1;
    }
}

int HGibbs_03() {
    try {

        Matrix r1(3, 1), r2(3, 1), r3(3, 1);
        r1(1,1) = 7000e3; r1(2,1) = 0;      r1(3,1) = 0;
        r2(1,1) = 5000e3; r2(2,1) = 5000e3; r2(3,1) = 0;
        r3(1,1) = 0;      r3(2,1) = 7000e3; r3(3,1) = 0;

        double Mjd1 = 60000.0;
        double Mjd2 = Mjd1 + 5.0/(60.0*24.0);
        double Mjd3 = Mjd2 + 5.0/(60.0*24.0);

        GibbsResult result = hgibbs(r1, r2, r3, Mjd1, Mjd2, Mjd3);

        _assert(result.error == "   angl > 1ø");

        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error en test_HGibbs_03: " << e.what() << std::endl;
        return 1;
    }
}

int HGibbs_04() {
    try {

        Matrix r1(3, 1), r2(3, 1), r3(3, 1);
        r1(1,1) = 7000e3; r1(2,1) = 0;      r1(3,1) = 0;
        r2(1,1) = 0;      r2(2,1) = 8000e3; r2(3,1) = 0;
        r3(1,1) = -9000e3; r3(2,1) = 0;     r3(3,1) = 0;


        double Mjd1 = 60000.0;
        double Mjd2 = Mjd1 + 30.0/(60.0*24.0);
        double Mjd3 = Mjd2 + 30.0/(60.0*24.0);

        GibbsResult result = hgibbs(r1, r2, r3, Mjd1, Mjd2, Mjd3);

        double TOL = 1e-6;
        _assert(result.error == "   angl > 1ø");
        _assert(fabs(result.theta - 1.570796) < TOL);
        _assert(fabs(result.theta1 - 1.570796) < TOL);
        _assert(fabs(result.copa - 0.0) < TOL);

        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error en test_HGibbs_04: " << e.what() << std::endl;
        return 1;
    }
}

int Elements_01() {
    // Ejemplo: órbita inclinada, a=7000km, e=0.1, i=45°
    double a = 7000e3;
    double e = 0.1;
    double i = M_PI/4;

    double p = a*(1-e*e);
    double r0 = p/(1+e);
    Matrix r(3,1), v(3,1);
    r(1,1)=r0; r(2,1)=0; r(3,1)=0;
    double h = std::sqrt(GM_Earth*p);
    v(1,1)=0; v(2,1)=h/r0; v(3,1)=0;

    Matrix Rinc(3,3);
    Rinc(1,1)=1; Rinc(1,2)=0;             Rinc(1,3)=0;
    Rinc(2,1)=0; Rinc(2,2)=cos(i); Rinc(2,3)=-sin(i);
    Rinc(3,1)=0; Rinc(3,2)=sin(i); Rinc(3,3)= cos(i);
    r = Rinc * r;
    v = Rinc * v;

    auto E = elements(r, v);
    _assert(std::fabs(E.a - a) < 1e-3);
    _assert(std::fabs(E.e - e) < 1e-6);
    _assert(std::fabs(E.i - i) < 1e-6);

    return 0;
}

int Elements_02() {

    Matrix r(3,1), v(3,1);
    r(1,1)=7000e3; r(2,1)=0; r(3,1)=0;
    double speed = std::sqrt(GM_Earth/r.norm());
    v(1,1)=0; v(2,1)=speed; v(3,1)=0;
    try {
        auto E = elements(r,v);
    } catch(const std::invalid_argument&) {
        return 0;
    }
    std::cout<<"Elements_02 failed\n";
    return 1;
}

int LTC_01() {
    try {

        double lon = 0.0;
        double lat = 0.0;

        Matrix result = LTC(lon, lat);


        Matrix expected(3, 3);
        expected(1,1) = 0.0; expected(1,2) = 1.0; expected(1,3) = 0.0;
        expected(2,1) = 0.0; expected(2,2) = 0.0; expected(2,3) = 1.0;
        expected(3,1) = 1.0; expected(3,2) = 0.0; expected(3,3) = 0.0;


        _assert(MatrixEqual(result, expected));

        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error en LTC_test_01: " << e.what() << std::endl;
        return 1;
    }
}

int LTC_02() {
    try {

        double lon = M_PI/6;      // 30°
        double lat = M_PI/4;      // 45°

        Matrix result = LTC(lon, lat);


        Matrix expected(3, 3);
        expected(1,1) = -0.500000; expected(1,2) = 0.866025; expected(1,3) = 0.000000;
        expected(2,1) = -0.612372; expected(2,2) = -0.353553; expected(2,3) = 0.707107;
        expected(3,1) = 0.612372;  expected(3,2) = 0.353553;  expected(3,3) = 0.707107;

        _assert(MatrixEqual(result, expected, 1e-6));

        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error en LTC_test_02: " << e.what() << std::endl;
        return 1;
    }
}

int GHAMatrix_01() {
    try {

        double Mjd_UT1 = 51544.5;  // J2000.0
        Matrix result = GHAMatrix(Mjd_UT1);


        double expected_GAST = gstime(Mjd_UT1);

        Matrix expected = R_z(expected_GAST);

        _assert(MatrixEqual(result, expected));

        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error en GHAMatrix_test_01: " << e.what() << std::endl;
        return 1;
    }
}

int GHAMatrix_02() {
    try {


        double Mjd_UT1 = 59754.0;
        Matrix result = GHAMatrix(Mjd_UT1);


        double expected_GAST = gstime(Mjd_UT1);


        Matrix expected = R_z(expected_GAST);

        _assert(MatrixEqual(result, expected));

        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error en GHAMatrix_test_02: " << e.what() << std::endl;
        return 1;
    }
}

int PoleMatrix_01() {
    try {



        double arcsec_to_rad = M_PI/(180.0*3600.0);
        double xp = 0.1 * arcsec_to_rad;
        double yp = 0.2 * arcsec_to_rad;

        Matrix result = PoleMatrix(xp, yp);

        Matrix expected = R_y(-xp) * R_x(-yp);

        _assert(MatrixEqual(result, expected, 1e-12));

        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error en PoleMatrix_test_02: " << e.what() << std::endl;
        return 1;
    }
}

int NutMatrix_01() {
    try {


        double Mjd_TT = 51544.5;  // J2000.0
        Matrix result = NutMatrix(Mjd_TT);

        double eps = MeanObliquity(Mjd_TT);
        double dpsi, deps;
        NutAngles(Mjd_TT, dpsi, deps);
        Matrix expected = R_x(-eps-deps) * R_z(-dpsi) * R_x(+eps);

        _assert(MatrixEqual(result, expected, 1e-12));

        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}

int PrecMatrix_01() {
    try {


        double Mjd_1 = 51544.5;  // J2000.0
        double Mjd_2 = 58952.0;  // 2020-03-21 12:00:00 TT

        Matrix result = PrecMatrix(Mjd_1, Mjd_2);


        Matrix expected(3, 3);
        expected(1,1) = 0.9999999999999929; expected(1,2) = -0.0000271361741621; expected(1,3) = -0.0000117760568353;
        expected(2,1) = 0.0000271361741621; expected(2,2) = 0.9999999999996318; expected(2,3) = -0.0000000000000000;
        expected(3,1) = 0.0000117760568353; expected(3,2) = -0.0000000000000000; expected(3,3) = 0.9999999999993063;



        _assert(MatrixEqual(result, expected, 1e-2));

        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error en PrecMatrix_01: " << e.what() << std::endl;
        return 1;
    }
}

int timediff_01() {
    try {

        double UT1_UTC = 0.0;
        double TAI_UTC = 0.0;

        TimeDiffs result = timediff(UT1_UTC, TAI_UTC);


        _assert(fabs(result.UT1_TAI - 0.0) < TOL_);
        _assert(fabs(result.UTC_GPS - 19.0) < TOL_);
        _assert(fabs(result.UT1_GPS - 19.0) < TOL_);
        _assert(fabs(result.TT_UTC - 32.184) < TOL_);
        _assert(fabs(result.GPS_UTC - (-19.0)) < TOL_);
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error en timediff_01: " << e.what() << std::endl;
        return 1;
    }
}

int Geodetic_01() {
    try {

        Matrix r(3, 1);
        r(1,1) = R_Earth; r(2,1) = 0.0; r(3,1) = 0.0;

        GeodeticCoords coords = Geodetic(r);

        _assert(fabs(coords.longitude - 0.0) < TOL_);
        _assert(fabs(coords.latitude - 0.0) < TOL_);
        _assert(fabs(coords.altitude - 0.0) < 1.0);

        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error en Geodetic_test_01: " << e.what() << std::endl;
        return 1;
    }
}

int Geodetic_02() {

    Matrix r(3,1);
    r(1,1) = R_Earth;
    r(2,1) = 0;
    r(3,1) = 0;
    auto G = Geodetic(r);
    if (std::abs(G.longitude - 0.0) > TOL_ ||
        std::abs(G.latitude  - 0.0) > TOL_ ||
        std::abs(G.altitude  - 0.0) > 1e-6) {
        std::cout<<"lon="<<G.longitude
                 <<" lat="<<G.latitude
                 <<" h="<<G.altitude<<"\n";
        return 1;
    }
    return 0;
}

int Geodetic_03() {

    Matrix r(3,1);
    r(1,1) = 0;
    r(2,1) = 0;
    r(3,1) = R_Earth * (1 - f_Earth);
    auto G = Geodetic(r);
    if (std::abs(G.longitude) > TOL_||
        std::abs(G.latitude  - M_PI_2) > 1e-9 ||
        std::abs(G.altitude) > 1e-6) {
        std::cout<<"lon="<<G.longitude
                 <<" lat="<<G.latitude
                 <<" h="<<G.altitude<<"\n";
        return 1;
    }
    return 0;
}

int doubler_01() {
    try {
        double los1_data[] = {0.2673, 0.5345, 0.8018};
        double los2_data[] = {0.4082, 0.8165, 0.4082};
        double los3_data[] = {0.7071, 0.0,    0.7071};

        double rsite1_data[] = {6378.137, 0.0, 0.0};
        double rsite2_data[] = {0.0, 6378.137, 0.0};
        double rsite3_data[] = {0.0, 0.0, 6378.137};

        Matrix los1(3,1, los1_data, 3);
        Matrix los2(3,1, los2_data, 3);
        Matrix los3(3,1, los3_data, 3);

        Matrix rsite1(3,1, rsite1_data, 3);
        Matrix rsite2(3,1, rsite2_data, 3);
        Matrix rsite3(3,1, rsite3_data, 3);

        double cc1 = 1.0;
        double cc2 = 1.0;
        double magrsite1 = rsite1.norm();
        double magrsite2 = rsite2.norm();
        double magr1in = 7000.0;
        double magr2in = 7100.0;
        double t1 = 10.0;
        double t3 = 40.0;
        char direct = 'y';

        DoubleRResult res = doubler(cc1, cc2, magrsite1, magrsite2, magr1in, magr2in,
                                    los1, los2, los3, rsite1, rsite2, rsite3,
                                    t1, t3, direct);


        _assert(res.r2.getFilas() == 3 && res.r2.getColumnas() == 1);
        _assert(res.r3.getFilas() == 3 && res.r3.getColumnas() == 1);
        _assert(!std::isnan(res.magr1) && !std::isinf(res.magr1));
        _assert(!std::isnan(res.magr2) && !std::isinf(res.magr2));
        _assert(!std::isnan(res.f1) && !std::isinf(res.f1));
        _assert(!std::isnan(res.f2) && !std::isinf(res.f2));
        _assert(!std::isnan(res.q1) && !std::isinf(res.q1));
        _assert(!std::isnan(res.a) && !std::isinf(res.a));
        _assert(!std::isnan(res.deltae32) && !std::isinf(res.deltae32));

        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error en doubler_01: " << e.what() << std::endl;
        return 1;
    }
}

int Doubler_01() {
    bool threw=false;
    try {
        doubler(
                0.0, 0.0,
                1.0, 3.0,
                1.0, 1.0,
                Matrix(3,1),Matrix(3,1),Matrix(3,1),
                Matrix(3,1),Matrix(3,1),Matrix(3,1),
                0.0, 0.0, 'y'
        );
    } catch(...) { threw=true; }
    _assert(threw && "Esperaba excepción por delta2<0");

    return 0;
}

int Doubler_02() {
    bool threw=false;

    Matrix los1(3,1),los2(3,1),los3(3,1);
    Matrix rs1(3,1), rs2(3,1), rs3(3,1);
    los1(1,1)=1; los1(2,1)=0; los1(3,1)=0;
    los2(1,1)=0; los2(2,1)=1; los2(3,1)=0;
    los3(1,1)=1; los3(2,1)=0; los3(3,1)=0;
    rs1(1,1)=1; rs1(2,1)=0; rs1(3,1)=0;
    rs2(1,1)=0; rs2(2,1)=1; rs2(3,1)=0;
    rs3 = Matrix(3,1);
    try {
        doubler(
                2.0, 2.0,
                1.0, 1.0,
                1.0, 1.0,
                los1, los2, los3,
                rs1, rs2, rs3,
                0.0, 0.0, 'y'
        );
    } catch(...) {
        threw=true;
    }
    _assert(threw && "Esperaba excepción por división por cero en rho3");
    return 0;
}

int Doubler_03() {

    double cc1 = 0.0, cc2 = 0.0;
    double magrsite1 = 0.0, magrsite2 = 0.0;
    double magr1in = 10000.0, magr2in = 10000.0;
    Matrix los1(3,1), los2(3,1), los3(3,1);
    los1(1,1)=1; los1(2,1)=0; los1(3,1)=0;
    los2(1,1)=0; los2(2,1)=1; los2(3,1)=0;
    los3(1,1)=0; los3(2,1)=0; los3(3,1)=1;
    Matrix rsite1(3,1), rsite2(3,1), rsite3(3,1);

    DoubleRResult R = doubler(
            cc1, cc2,
            magrsite1, magrsite2,
            magr1in, magr2in,
            los1, los2, los3,
            rsite1, rsite2, rsite3,
            1.0, 1.0, 'y'
    );
    double expected_rho2 = (-cc2 + sqrt(cc2*cc2 - 4*(magrsite2*magrsite2 - magr2in*magr2in))) / 2.0;
    if (std::fabs(R.r2(1,1) - 0.0) > TOL_ ||
        std::fabs(R.r2(2,1) - expected_rho2) > TOL_ ||
        std::fabs(R.r2(3,1) - 0.0) > TOL_) {
        std::cout<<"Doubler_03 FAILED\n";
        return 1;
    }
    return 0;
}

int IERS_01() {


    IERSResult r = IERS(eopdata, 58000.0, 'l');

    _assert(fabs(r.x_pole   - 0.000001151607026) < TOL_);
    _assert(fabs(r.y_pole   - 0.000001704629143) < TOL_);
    _assert(fabs(r.UT1_UTC  - 0.3357933) < TOL_);
    _assert(fabs(r.LOD      - 0.0003011) < TOL_);
    _assert(fabs(r.dpsi     + 0.000000538947977) < TOL_); // signo negativo
    _assert(fabs(r.deps     + 0.000000060611406) < TOL_); // signo negativo
    _assert(fabs(r.dx_pole  - 0.000000000799943) < TOL_);
    _assert(fabs(r.dy_pole  - 0.000000000911450) < TOL_);
    _assert(fabs(r.TAI_UTC  - 37.0) < TOL_);
    return 0;
}

int JPL_Eph_01() {
    double Mjd_TDB = 106351.0;

    PlanetaryPositions p = JPL_Eph_DE430(Mjd_TDB);

    double tol = 10e-6;
    _assert(fabs(p.r_Earth(1,1)    - (-31598008156.192417)) < tol);
    _assert(fabs(p.r_Earth(2,1)    - (132460433151.290695)) < tol);
    _assert(fabs(p.r_Earth(3,1)    - (57373501576.184982))  < tol);

    _assert(fabs(p.r_Sun(1,1)      - (31713109510.304195))  < tol);
    _assert(fabs(p.r_Sun(2,1)      - (-131813766788.852325))< tol);
    _assert(fabs(p.r_Sun(3,1)      - (-57095479016.666679)) < tol);

    _assert(fabs(p.r_Moon(1,1)     - (-364222823.780679))   < tol);
    _assert(fabs(p.r_Moon(2,1)     - (-167003158.299063))   < tol);
    _assert(fabs(p.r_Moon(3,1)     - (-33007954.970030))    < tol);

    _assert(fabs(p.r_Mercury(1,1)  - (-26806176952.609943)) < tol);
    _assert(fabs(p.r_Mercury(2,1)  - (-138002470934.937256))< tol);
    _assert(fabs(p.r_Mercury(3,1)  - (-54357817587.275719)) < tol);

    _assert(fabs(p.r_Venus(1,1)    - (-27325902768.835911)) < tol);
    _assert(fabs(p.r_Venus(2,1)    - (-51364561590.153427)) < tol);
    _assert(fabs(p.r_Venus(3,1)    - (-17133255020.696655)) < tol);

    _assert(fabs(p.r_Mars(1,1)     - (16285514723.482100))  < tol);
    _assert(fabs(p.r_Mars(2,1)     - (-330699687867.098938))< tol);
    _assert(fabs(p.r_Mars(3,1)     - (-147917974938.416229))< tol);

    _assert(fabs(p.r_Jupiter(1,1)  - (-66644823641.162491)) < tol);
    _assert(fabs(p.r_Jupiter(2,1)  - (-853591820524.296997))< tol);
    _assert(fabs(p.r_Jupiter(3,1)  - (-364009407514.953491))< tol);

    _assert(fabs(p.r_Saturn(1,1)   - (184661249713.726105)) < tol);
    _assert(fabs(p.r_Saturn(2,1)   - (1109215158423.234863))< tol);
    _assert(fabs(p.r_Saturn(3,1)   - (449187708416.220093)) < tol);

    _assert(fabs(p.r_Uranus(1,1)   - (-1320196279348.428955)) < tol);
    _assert(fabs(p.r_Uranus(2,1)   - (-2405446527895.510742))< tol);
    _assert(fabs(p.r_Uranus(3,1)   - (-1033761148857.640991))< tol);

    _assert(fabs(p.r_Neptune(1,1)  - (183888341140.650604)) < tol);
    _assert(fabs(p.r_Neptune(2,1)  - (-4316591998386.310547))< tol);
    _assert(fabs(p.r_Neptune(3,1)  - (-1773748114193.736328))< tol);

    _assert(fabs(p.r_Pluto(1,1)    - (1327605134864.141602)) < tol);
    _assert(fabs(p.r_Pluto(2,1)    - (6519854043373.274414)) < tol);
    _assert(fabs(p.r_Pluto(3,1)    - (1628319706193.649658)) < tol);

    return 0;
}

int JPL_Eph_02() {
    PC = Matrix(1, 10);
    PC(1,1) = 1000.0;
    PC(1,2) = 2000.0;

    bool threw = false;
    try {
        JPL_Eph_DE430(3000.0);
    } catch(const std::runtime_error& e) {
        threw = true;
        _assert(std::string(e.what()).find("JD fuera de rango")!=std::string::npos);
    }
    _assert(threw);
    return 0;
}

int AzElPa_Test_01() {
    try {

        Matrix s(3, 1);
        s(1,1) = 0.0;
        s(2,1) = 1.0;
        s(3,1) = 0.0;

        AzElPaResult result = AzElPa(s);


        double expected_Az = 0.0;
        double expected_El = 0.0;
        Matrix expected_dAds(1, 3);
        expected_dAds(1,1) = 1.0;
        expected_dAds(1,2) = 0.0;
        expected_dAds(1,3) = 0.0;
        Matrix expected_dEds(1, 3);
        expected_dEds(1,1) = 0.0;
        expected_dEds(1,2) = 0.0;
        expected_dEds(1,3) = 1.0;

        _assert(fabs(result.Az - expected_Az) < 1e-9);
        _assert(fabs(result.El - expected_El) < 1e-9);
        /*
        result.dAds.print();
        std::cout<<"==================="<<std::endl;
        expected_dAds.print();
        std::cout<<"==================="<<std::endl;
        result.dEds.print();
        std::cout<<"==================="<<std::endl;
        expected_dEds.print();
        std::cout<<"==================="<<std::endl;*/



        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error in AzElPa_Test_01: " << e.what() << std::endl;
        return 1;
    }
}

int AzElPa_Test_02() {
    try {

        Matrix s(3, 1);
        s(1,1) = 1.0;
        s(2,1) = 1.0;
        s(3,1) = sqrt(2.0);  // Up

        AzElPaResult result = AzElPa(s);

        // Expected results
        double expected_Az = M_PI/4.0; // 45°
        double expected_El = M_PI/4.0; // 45°


        _assert(fabs(result.Az - expected_Az) < 1e-9);
        _assert(fabs(result.El - expected_El) < 1e-9);

        _assert(fabs(result.dAds(1,1) - 0.5) < 1e-9);
        _assert(fabs(result.dAds(1,2) - (-0.5)) < 1e-9);
        _assert(fabs(result.dAds(1,3) - 0.0) < 1e-9);

        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error in AzElPa_Test_02: " << e.what() << std::endl;
        return 1;
    }
}

int AzElPa_Test_03() {
    try {

        Matrix s(2, 2);

        try {
            AzElPaResult result = AzElPa(s);
            std::cerr << "Error: Expected exception not thrown\n";
            return 1;
        } catch (const std::invalid_argument& e) {
            _assert(std::string(e.what()).find("3x1") != std::string::npos);
        }

        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error in AzElPa_Test_03: " << e.what() << std::endl;
        return 1;
    }
}

int VarEqn_Test_01() {
    double Mjd0 = 58000.0;

    AuxParam params;
    params.Mjd_UTC = Mjd0;
    params.Mjd_TT  = Mjd0;
    params.n = 0;
    params.m = 0;

    Matrix yPhi(42,1);
    for(int i=1;i<=6;++i) yPhi(i,1)=0.0;

    for(int j=1;j<=6;++j) {
        yPhi(6*j+ j, 1) = 1.0;
    }


    Matrix yPhip = VarEqn(0.0, yPhi, params, eopdata);

    _assert(std::fabs(yPhip(1,1)) < 1e-12);
    _assert(std::fabs(yPhip(6,1)) < 1e-2);

    return 0;
}

int VarEqn_Test_02() {
    Matrix yPhi_bad(10,1);
    AuxParam params{59000.0, 59000.0, 0, 0};

    double JD = params.Mjd_UTC + 2400000.5;

    bool threw = false;
    try {
        Matrix out = VarEqn(0.0, yPhi_bad, params, eopdata);
    } catch (const std::exception& e) {
        threw = true;
    }
    _assert(threw && "Expected exception for wrong yPhi size");

    return 0;
}

int VarEqn_Test_03() {
    double Mjd0 = 51544.5;

    AuxParam params;
    params.Mjd_UTC = Mjd0;
    params.Mjd_TT  = Mjd0;
    params.n = 0; params.m = 0;


    const double R = 7000e3;
    Matrix yPhi(42,1);
    // r
    yPhi(1,1)= R; yPhi(2,1)=0; yPhi(3,1)=0;
    // v
    yPhi(4,1)=0;  yPhi(5,1)=0; yPhi(6,1)=0;
    // Phi = identidad
    for(int j=1;j<=6;++j) yPhi(6*j+j,1)=1.0;

    Matrix yPhip = VarEqn(0.0, yPhi, params, eopdata);

    for(int i=1;i<=3;++i) _assert(fabs(yPhip(i,1))<1e-12);
    double expect_ax = -GM_Earth/(R*R);
    _assert(fabs(yPhip(4,1)-expect_ax)<1e-6);
    _assert(fabs(yPhip(5,1))<1e-12 && fabs(yPhip(6,1))<1e-12);
    return 0;
}

int Accel_01() {
    AuxParam p{58000.0, 58000.0, 0, 0, false,false,false};

    double R = 7000e3;
    Matrix Y(6,1);
    Y(1,1)=R;
    Matrix dY = Accel(0.0, Y, p, eopdata);
    double expect = -GM_Earth/(R*R);
    _assert(fabs(dY(4,1)-expect)<1e-8);
    _assert(fabs(dY(5,1))<1e-12);
    _assert(fabs(dY(6,1))<1e-12);
    return 0;
}

int Accel_02(){
    AuxParam params{58000.0,58000.0,0,0};
    Matrix Y_bad(5,1);
    bool threw=false;
    try {
        Accel(0.0, Y_bad, params, eopdata);
    } catch(const std::invalid_argument&) {
        threw=true;
    }
    _assert(threw);
    return 0;
}

int Accel_03(){
    AuxParam params{58000.0,58000.0,2,3};

    Matrix Y(6,1);
    for(int i=1;i<=6;++i) Y(i,1)=1.0;
    bool threw=false;
    try {
        Accel(0.0, Y, params, eopdata);
    } catch(const std::invalid_argument&) {
        threw=true;
    }
    _assert(threw);
    return 0;
}

int Accel_04(){
    AuxParam params{58000.0,58000.0,0,0};
    params.sun = params.moon = params.planets = false;

    Matrix Y(6,1);
    for(int i=1;i<=6;++i) Y(i,1)=0.0;
    Matrix dY = Accel(0.0, Y, params, eopdata);

    for(int i=1;i<=6;++i){
        _assert(std::fabs(dY(i,1))<1e-12);
    }

    return 0;
}

int Accel_05(){
    AuxParam params{58000.0,58000.0,0,0};
    params.sun = params.moon = params.planets = false;

    const double R = 7000e3;
    const double mu = GM_Earth;
    const double v_circ = std::sqrt(mu/(R));

    Matrix Y(6,1);
    Y(1,1)= R;  Y(2,1)=0;      Y(3,1)=0;
    Y(4,1)=0;   Y(5,1)=v_circ; Y(6,1)=0;

    Matrix dY = Accel(0.0, Y, params, eopdata);

    _assert(std::fabs(dY(1,1) - 0.0)     < 1e-12);
    _assert(std::fabs(dY(2,1) - v_circ)  < 1e-8);
    _assert(std::fabs(dY(3,1) - 0.0)     < 1e-12);

    double expect_ax = -mu/(R*R);
    _assert(std::fabs(dY(4,1) - expect_ax) < 1e-6);
    _assert(std::fabs(dY(5,1) - 0.0)       < 1e-12);
    _assert(std::fabs(dY(6,1) - 0.0)       < 1e-12);

    return 0;
}

int anglesdr_01(){
    const double R = R_Earth + 1000e3;
    const double GM = GM_Earth;
    double w = std::sqrt(GM/(R*R*R));
    double dt = 60.0;
    double M2 = 58000.0;
    double M1 = M2 - dt/86400.0, M3 = M2 + dt/86400.0;

    Matrix site(3,1); site(1,1)=0; site(2,1)=0; site(3,1)=0;
    auto sat = [&](double t){
        Matrix r(3,1);
        r(1,1)=R*std::cos(w*t);
        r(2,1)=R*std::sin(w*t);
        r(3,1)=0;
        return r;
    };
    Matrix r1 = sat(-dt), r2 = sat(0), r3 = sat(+dt);
    Matrix los1 = (r1)*(1.0/r1.norm());
    Matrix los2 = (r2)*(1.0/r2.norm());
    Matrix los3 = (r3)*(1.0/r3.norm());
    double az1 = std::atan2(los1(1,1),los1(2,1)), el1 = std::asin(los1(3,1));
    double az2 = std::atan2(los2(1,1),los2(2,1)), el2 = std::asin(los2(3,1));
    double az3 = std::atan2(los3(1,1),los3(2,1)), el3 = std::asin(los3(3,1));

    AuxParam params;
    params.Mjd_UTC = M2; params.Mjd_TT = M2;
    params.n=0; params.m=0;
    params.sun = params.moon = params.planets = false;


    auto out = anglesdr(
            az1,az2,az3,
            el1,el2,el3,
            M1,M2,M3,
            site,site,site,
            params,
            eopdata
    );

    Matrix expect_r2 = r2;
    Matrix expect_v2(3,1);
    expect_v2(1,1) = -R*w*std::sin(0.0);
    expect_v2(2,1) =  R*w*std::cos(0.0);
    expect_v2(3,1) =  0.0;

    bool ok = true;
    for(int i=1;i<=3;++i){
        if(std::fabs(out.r2(i,1)-expect_r2(i,1))>1e-2) ok=false;
        if(std::fabs(out.v2(i,1)-expect_v2(i,1))>1e-4) ok=false;
    }
    if(!ok){
        std::cout<<"anglesdr_01 FAILED\n";
        return 1;
    }
    return 0;
}

static bool is_finite(double x) {
    return std::isfinite(x);
}

int anglesg_01() {
    double Mjd0 = 58000.0;


    double az1=0.1, az2=0.2, az3=0.3;
    double el1=0.15, el2=0.25, el3=0.35;
    double Mjd1=Mjd0, Mjd2=Mjd0+0.01, Mjd3=Mjd0+0.02;

    Matrix Rs1(3,1), Rs2(3,1), Rs3(3,1);
    Rs1(1,1)=1000; Rs1(2,1)=2000; Rs1(3,1)=3000;
    Rs2(1,1)=1100; Rs2(2,1)=2100; Rs2(3,1)=3100;
    Rs3(1,1)=1200; Rs3(2,1)=2200; Rs3(3,1)=3200;

    AuxParam params;
    params.Mjd_UTC = Mjd2;
    params.Mjd_TT  = Mjd2;

    try {
        AnglesGResult out = anglesg(
                az1, az2, az3,
                el1, el2, el3,
                Mjd1, Mjd2, Mjd3,
                Rs1, Rs2, Rs3,
                params, eopdata
        );

        if (out.r2.getFilas()!=3 || out.r2.getColumnas()!=1) {
            std::cerr<<"r2 dimension wrong\n";
            return 1;
        }
        if (out.v2.getFilas()!=3 || out.v2.getColumnas()!=1) {
            std::cerr<<"v2 dimension wrong\n";
            return 1;
        }

        for(int i=1;i<=3;++i){
            if (!is_finite(out.r2(i,1))) {
                std::cerr<<"r2("<<i<<") is not finite\n";
                return 1;
            }
            if (!is_finite(out.v2(i,1))) {
                std::cerr<<"v2("<<i<<") is not finite\n";
                return 1;
            }
        }

        return 0;
    }
    catch(const std::exception &ex) {
        std::cerr<<"anglesg_BasicFinite_Test threw: "<<ex.what()<<"\n";
        return 1;
    }
}

int anglesg_02() {

    double M2 = 58000.0;
    Matrix Rs(3,1); Rs(1,1)=0; Rs(2,1)=0; Rs(3,1)=0;
    AuxParam params;
    params.Mjd_UTC = M2;
    params.Mjd_TT  = M2;

    Matrix eop_bad(12,3);
    bool threw = false;
    try {
        anglesg(0,0,0, 0,0,0, M2,M2,M2, Rs,Rs,Rs, params, eop_bad);
    } catch(const std::exception&) {
        threw = true;
    }
    _assert(threw && "Expected exception for bad eopdata size");
    return 0;
}

int anglesg_03() {
    double M2 = 58000.0;
    Matrix Rs(3,1); Rs(1,1)=0; Rs(2,1)=0; Rs(3,1)=0;
    AuxParam params;
    params.Mjd_UTC = M2;
    params.Mjd_TT  = M2;
    Matrix eop_bad(12,3);
    bool threw = false;
    try {
        anglesg(0,0,0, 0,0,0, M2,M2,M2, Rs,Rs,Rs, params, eop_bad);
    } catch(const std::exception&) {
        threw = true;
    }
    _assert(threw && "Expected exception for bad eopdata size");
    return 0;
}

int DEInteg_01() {
    DEInteg integ;

    auto f = [&](double t, const Matrix& y) {
        Matrix yp(1,1);
        yp(1,1) = 123.0;  // no importa
        return yp;
    };

    Matrix y0(1,1);
    y0(1,1) = 42.0;
    double t0   = 5.0;
    double tOut = 5.0;
    Matrix y = y0;

    Matrix y1 = integ.integrate(f, t0, tOut, 1e-8, 1e-8, y);
    if (std::fabs(y1(1,1) - 42.0) > TOL_) {
        std::cout << " got " << y1(1,1) << ", expected 42\n";
        return 1;
    }
    return 0;
}

int DEInteg_02() {
    DEInteg solver;
    double t0 = 0.0, t1 = M_PI/2;
    Matrix y(2,1);
    y(1,1) = 1.0; // y1(0)
    y(2,1) = 0.0; // y2(0)
    auto f = [&](double t, const Matrix& ys)->Matrix {
        Matrix dy(2,1);
        dy(1,1) = ys(2,1);
        dy(2,1) = -ys(1,1);
        return dy;
    };
    Matrix yout = solver.integrate(f, t0, t1, 1e-8, 1e-10, y);
    double y1_exp = 1.0, y2_exp = 0.0;
    if (std::fabs(yout(1,1) - y1_exp) > 1e-2 ||
        std::fabs(yout(2,1) - y2_exp) > 1e-2) {
        std::cout << "DEInteg_SHO_Test failed: got ("
                  << yout(1,1)<<","<<yout(2,1)<<"), expected ("
                  << y1_exp<<","<<y2_exp<<")\n";
        return 1;
    }
    return 0;
}









int all_tests()
{

    _verify(Matrix_Basico);
    _verify(Mjday_01);
    _verify(Mjday_02);
    _verify(R_x_01);
    _verify(R_y_01);
    _verify(R_z_01);
    _verify(TimeUpdate_01);
    _verify(Position_01);
    _verify(Position_02);
    _verify(Legendre_01);
    _verify(sign_);
    _verify(AccelPointMass_01);
    _verify(Mjday_TDB_01);
    _verify(Mjday_TDB_02);
    _verify(Mjday_TDB_03);
    _verify(angl_01);
    _verify(angl_02);
    _verify(Cheb3D_01);
    _verify(Cheb3D_02);
    _verify(MeanObliquity_01);
    _verify(NutAngles_01);
    _verify(NutAngles_02);
    _verify(NutAngles_03);
    _verify(AccelHarmonic_01);
    _verify(AccelHarmonic_02);
    _verify(AccelHarmonic_03);
    _verify(AccelHarmonic_04);
    _verify(G_AccelHarmonic_01);
    _verify(G_AccelHarmonic_02);
    _verify(EqnEquinox_01);
    _verify(EqnEquinox_02);
    _verify(MeasUpdate_01);
    _verify(MeasUpdate_02);
    _verify(gstime_01);
    _verify(unit_01);
    _verify(unit_02);
    _verify(unit_03);
    _verify(Doubler_01);
    _verify(Doubler_02);
    _verify(Doubler_03);
    _verify(Gibbs_01);
    _verify(HGibbs_01);
    _verify(HGibbs_02);
    _verify(HGibbs_03);
    _verify(HGibbs_04);
    _verify(Elements_01);
    _verify(Elements_02);
    _verify(LTC_01);
    _verify(LTC_02);
    _verify(GHAMatrix_01);
    _verify(GHAMatrix_02);
    _verify(PoleMatrix_01);
    _verify(NutMatrix_01);
    _verify(PrecMatrix_01);
    _verify(timediff_01);
    _verify(Geodetic_01);
    _verify(Geodetic_02);
    _verify(Geodetic_03);
    _verify(doubler_01);
    _verify(IERS_01);
    _verify(JPL_Eph_01);
    _verify(JPL_Eph_02);
    _verify(AzElPa_Test_01);
    _verify(AzElPa_Test_02);
    _verify(AzElPa_Test_03);
    _verify(VarEqn_Test_01);
    _verify(VarEqn_Test_02);
    _verify(VarEqn_Test_03);
    _verify(Accel_01);
    _verify(Accel_02);
    _verify(Accel_03);
    _verify(Accel_04);
    _verify(Accel_05);
    _verify(anglesdr_01);
    _verify(anglesg_01);
    _verify(anglesg_02);
    _verify(anglesg_03);
    _verify(DEInteg_01);
    _verify(DEInteg_02);


    return 0;
}




int main()
{
    try {
        cargarPC("../data/DE430Coeff.txt");
    } catch (const std::exception& e) {
        std::cerr << "Error al cargar PC: " << e.what() << std::endl;
    }
    try{
        cargarEOP("../data/eop19620101.txt");
    }catch(const std::exception& e){
        std::cerr << "Error al cargar eopdata: " << e.what() << std::endl;
    }
    int result = all_tests();

    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);
/*
    double v[] = {1.0, 2.0, 3.0, 4.0};
    double v2[] = {5.0, 6.0, 7.0, 8.0};

    Matrix m1(2, 2);
    m1.print();

    Matrix m2(2, 2, v, 4);
    m2.print();

    Matrix m3(2, 2, v2, 4);
    m3.print();

    m1 = m2 * m3 * m2;
    m1.print();
 */
    return result != 0;
}
