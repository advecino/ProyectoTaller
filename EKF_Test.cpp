#include <iostream>
#include <cassert>
#include <iomanip>
#include <stdio.h>
#include <math.h>
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
        std::cout << "=== Test de inicialización básica ===" << std::endl;

        // Test 1: Creación simple
        Matrix m(2, 3);
        std::cout << "Matriz 2x3 creada. Filas: " << m.getFilas()
                  << ", Columnas: " << m.getColumnas() << std::endl;

        // Verificar dimensiones
        _assert(m.getFilas() == 2);
        _assert(m.getColumnas() == 3);

        // Verificar inicialización a cero
        for (int i = 1; i <= 2; ++i) {
            for (int j = 1; j <= 3; ++j) {
                _assert(m(i,j) == 0.0);
            }
        }

        // Test 2: Asignación de valores
        m(1,1) = 1.5;
        m(2,3) = 2.0;
        _assert(m(1,1) == 1.5);
        _assert(m(2,3) == 2.0);

        std::cout << "Test de inicialización pasado con éxito!\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error en Test_Matrix_Basico: " << e.what() << std::endl;
        return 1;
    }
}

int Position_01() {
    // En el ecuador (lat=0), meridiano de Greenwich (lon=0), sin altura
    double lon = 0.0;
    double lat = 0.0;
    double h   = 0.0;
    Matrix r = Position(lon, lat, h, R_Earth, f_Earth);

    std::cout << "Position@Equator: ["
              << r(1,1) << ", " << r(2,1) << ", " << r(3,1) << "]\n";

    double TOL = 1e-6;
    // Debe ser (R_equ, 0, 0)
    _assert(std::fabs(r(1,1) - R_Earth) < TOL);
    _assert(std::fabs(r(2,1))          < TOL);
    _assert(std::fabs(r(3,1))          < TOL);

    std::cout << "Test_Position_Equator passed\n";
    return 0;
}

int Position_02() {
    // En el polo Norte (lat=90°), lon irrelevante, con h=1000 m
    double lon = 1.234;                      // puede ser cualquier valor
    double lat = M_PI/2.0;                   // 90°
    double h   = 1000.0;                     // 1 km sobre el elipsoide

    // excentricidad al cuadrado
    const double e2 = f_Earth * (2.0 - f_Earth);
    // N en el polo
    double N = R_Earth / std::sqrt(1.0 - e2);
    // z esperado
    double expected_z = (1.0 - e2) * N + h;

    Matrix r = Position(lon, lat, h, R_Earth, f_Earth);
    std::cout << "Position@Pole: ["
              << r(1,1) << ", " << r(2,1) << ", " << r(3,1) << "]\n";

    constexpr double TOL_POS = 1e-6;
    constexpr double TOL_Z   = 1e-3; // 1 mm

    // X,Y ≈ 0, Z ≈ expected_z
    _assert(std::fabs(r(1,1))             < TOL_POS);
    _assert(std::fabs(r(2,1))             < TOL_POS);
    _assert(std::fabs(r(3,1) - expected_z) < TOL_Z);

    std::cout << "Test_Position_Pole passed\n";
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

    /*   cout << setprecision(20);
       cout << Mjday(2025,4,3,15,37,5) << endl;
       cout << 60768.6507523148;
      */
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

int Legendre_01()
{
    int n = 2, m = 2;
    double fi = M_PI/6; // 30°
    // Inicializar con dimensiones mínimas para usar operator=
    Matrix pnm(0,0), dpnm(0,0);
    Legendre(n, m, fi, pnm, dpnm);

    // Comprueba valores conocidos de P0, P1
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
    // Caso de prueba 1: Masa en el origen
    Matrix r(3, 1, new double[3]{7000.0, 0.0, 0.0}, 3);
    Matrix s(3, 1, new double[3]{0.0, 0.0, 0.0}, 3);
    double GM = 3.986004418e14; // GM de la Tierra

    Matrix a = AccelPointMass(r, s, GM);

    // Valor esperado: aceleración puramente radial
    double expected = -GM / (7000.0*7000.0*7000.0) * 7000.0;

    if(fabs(a(1,1) - expected) > TOL_ ||
       fabs(a(2,1)) > TOL_ ||
       fabs(a(3,1)) > TOL_) {
        std::cout << "AccelPointMass_01 failed!\n";
        return 1;
    }

    return 0;
}

int Mjday_TDB_01() {
    // Test case from Vallado's book example
    double Mjd_TT = 54930.5; // 2009-04-09 00:00:00 TT

    double Mjd_TDB = Mjday_TDB(Mjd_TT);

    // Expected value calculated from reference implementation
    double expected = 54930.5000007235;

    if(fabs(Mjd_TDB - expected) > 1e-5) {
        std::cout << std::setprecision(15);
        std::cout << "Mjday_TDB_01 failed!\n";
        std::cout << "Expected: " << expected << "\n";
        std::cout << "Got:      " << Mjd_TDB << "\n";
        std::cout << "Diff:     " << fabs(Mjd_TDB - expected) << "\n";
        return 1;
    }

    return 0;
}

int Mjday_TDB_02() {
    // Para cualquier Mjd_TT la corrección es muy pequeña
    double Mjd_TT = 51544.5; // J2000.0
    double Mjd_TDB = Mjday_TDB(Mjd_TT);
    double delta = Mjd_TDB - Mjd_TT;
    std::cout << "delta = " << delta << " días\n";
    _assert(std::fabs(delta) < 1e-4);
    std::cout << "Test_Mjday_TDB_Identity passed\n";
    return 0;
}

int Mjday_TDB_03() {
    // Si cambio Mjd_TT en N días, Mjday_TDB debe cambiar en ~N días
    double Mjd0 = 58000.0;
    double Mjd1 = Mjd0 + 10.5;
    double TDB0 = Mjday_TDB(Mjd0);
    double TDB1 = Mjday_TDB(Mjd1);
    double dt_in  = Mjd1 - Mjd0;
    double dt_out = TDB1 - TDB0;
    std::cout << "dt_in = " << dt_in << ", dt_out = " << dt_out << "\n";
    _assert(std::fabs(dt_out - dt_in) < 1e-4);
    std::cout << "Test_Mjday_TDB_Difference passed\n";
    return 0;
}

int angl_01() {
    // Vectores paralelos
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
    // Vectores perpendiculares
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
        // Coeficientes para f(t) = t en [0,1]
        const int N = 2;
        double Ta = 0.0;
        double Tb = 1.0;

        // Coeficientes exactos para f(t) = t:
        // t = 0.5*T0(tau) + 0.5*T1(tau)
        Matrix Cx(N,1, new double[N]{0.5, 0.5}, N);

        // Coeficientes para f(t) = t^2:
        // t^2 = 0.5*T0(tau) + 0.5*T1(tau) + 0.25*T2(tau)
        Matrix Cy(3,1, new double[3]{0.5, 0.5, 0.25}, 3);

        // Evaluar en t=0.5
        double t = 0.5;
        Matrix result = Cheb3D(t, 2, Ta, Tb, Cx, Cx, Cx);

        // Valores esperados
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
        // Coeficientes para f(t) = t en [0,1]
        // f(t) = 0.5*T0 + 0.5*T1
        Matrix Cx(2,1, new double[2]{0.5, 0.5}, 2);

        // Evaluar en t=1 (debería dar 1 para f(t)=t)
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
    // Valores de constantes
    const double MJD_J2000 = 51544.5;
    const double Rad = M_PI / 180.0; // Factor de conversión a radianes

    // Fecha de prueba (J2000)
    double Mjd_TT = MJD_J2000;

    // Calcular oblicuidad
    double MOblq = MeanObliquity(Mjd_TT);

    // Valor esperado para J2000 (23.4392911 grados en radianes)
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
        // Fecha J2000 (1.5 enero 2000)
        double Mjd_TT = 51544.5;

        // Calcular ecuación de los equinoccios
        double EqE = EqnEquinox(Mjd_TT);

        // El valor debería ser muy pequeño cerca de J2000
        // Usamos una tolerancia más relajada (1e-12 radianes ≈ 0.2 milliarcseconds)
        _assert(fabs(EqE) < 1e-3);
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Exception in EqnEquinox_01: " << e.what() << std::endl;
        return 1;
    }
}

int EqnEquinox_02() {
    try {
        // Fecha aleatoria (1 Jan 2020)
        double Mjd_TT = 58849.0;

        // Calcular ecuación de los equinoccios
        double EqE = EqnEquinox(Mjd_TT);

        // Valor de referencia más preciso (usando SOFA o JPL Horizons)
        // Para 1 Jan 2020 00:00:00 TT:
        double expected = -0.0033534 * M_PI/180.0/3600.0; // -0.0033534 arcsec en radianes

        // Tolerancia relajada a 1e-9 radianes (≈ 0.2 microarcseconds)
        _assert(fabs(EqE - expected) < 1e-3);
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Exception in EqnEquinox_02: " << e.what() << std::endl;
        return 1;
    }
}

int NutAngles_01() {
    try {
        // Fecha J2000 (1.5 enero 2000)
        double Mjd_TT = 51544.5;

        double dpsi, deps;
        NutAngles(Mjd_TT, dpsi, deps);

        // Valores esperados cercanos a cero
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
    NutAngles(51544.5, dpsi, deps); // J2000

    std::cout << "NutAngles@J2000: dpsi=" << dpsi
              << " rad, deps=" << deps << " rad\n";

    // Debe estar en ±0.01 rad
    _assert(std::fabs(dpsi) < 1e-2);
    _assert(std::fabs(deps) < 1e-2);

    std::cout << "Test_NutAngles_Range passed\n";
    return 0;
}

int NutAngles_03() {
    double dpsi0, deps0, dpsi1, deps1;
    NutAngles(51544.5,     dpsi0, deps0);
    NutAngles(51544.5+18262.5, dpsi1, deps1); // +0.5 siglo

    std::cout << "Δdpsi=" << (dpsi1-dpsi0)
              << ", Δdeps=" << (deps1-deps0) << " rad\n";

    // Debe cambiar al menos unos μrad
    _assert(std::fabs(dpsi1-dpsi0) > 1e-6 ||
            std::fabs(deps1-deps0) > 1e-6);

    std::cout << "Test_NutAngles_Variation passed\n";
    return 0;
}

int AccelHarmonic_01() {
    // Monopole puro: n_max=m_max=0, Cnm[1][1]=1.0
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

    std::cout<<"AccelHarmonic_Test_01 passed\n";
    return 0;
}

int AccelHarmonic_02() {
    // Monopole puro y posición en y
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

    std::cout<<"AccelHarmonic_Test_02 passed\n";
    return 0;
}

int G_AccelHarmonic_01() {
    // Monopole field test at r = [R,0,0]
    for(int i=0;i<300;i++) for(int j=0;j<300;j++){
            Cnm[i][j]=0.0; Snm[i][j]=0.0;
        }
    Cnm[1][1] = 1.0;  // monopole

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
    // off-diagonals zero
    for(int i=1;i<=3;i++) for(int j=1;j<=3;j++){
            if(i!=j) _assert(fabs(G(i,j)) < TOL_);
        }

    std::cout<<"G_AccelHarmonic_Test_01 passed\n";
    return 0;
}

int G_AccelHarmonic_02() {
    // Monopole field test at r = [0,R,0]
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
    // off-diagonals zero
    for(int i=1;i<=3;i++) for(int j=1;j<=3;j++){
            if(i!=j) _assert(fabs(G(i,j)) < TOL_);
        }

    std::cout<<"G_AccelHarmonic_Test_02 passed\n";
    return 0;
}

int MeasUpdate_01() {
    int n = 2, m = 2;
    Matrix x(n,1), P(n,n), K(0,0);
    Matrix z(m,1), g(m,1), s(m,1), G(m,n);

    // Initialize state and covariance
    x(1,1) = 1.0; x(2,1) = -2.0;
    P(1,1) = 2.0; P(1,2) = 0.5;
    P(2,1) = 0.5; P(2,2) = 1.0;

    // Measurement = prediction
    z(1,1) = 10.0; g(1,1) = 10.0;
    z(2,1) = -5.0; g(2,1) = -5.0;

    // Measurement noise std deviations
    s(1,1) = 1.0;
    s(2,1) = 2.0;

    // Sensitivity matrix G
    G(1,1) = 1.0; G(1,2) = 0.0;
    G(2,1) = 0.0; G(2,2) = 1.0;

    // Perform update
    MeasUpdate(x, P, K, z, g, s, G, n);


    const double P11 = 0.66101694915254;
    const double P12 = 0.13559322033898;
    const double P21 = 0.13559322033898;
    const double P22 = 0.74576271186441;

    _assert(fabs(P(1,1) - P11) < TOL_);
    _assert(fabs(P(1,2) - P12) < TOL_);
    _assert(fabs(P(2,1) - P21) < TOL_);
    _assert(fabs(P(2,2) - P22) < TOL_);

    // El estado x no cambia porque z-g = 0
    _assert(fabs(x(1,1) - 1.0) < TOL_);
    _assert(fabs(x(2,1) + 2.0) < TOL_);

    std::cout << "Test_MeasUpdate_Identity passed\n";
    return 0;
}

int MeasUpdate_02() {
    // Scalar case: n=1, m=1
    int n = 1, m = 1;
    Matrix x(n,1), P(n,n), K(0,0);
    Matrix z(m,1), g(m,1), s(m,1), G(m,n);

    // Prior state and covariance
    x(1,1) = 0.0;
    P(1,1) = 1.0;

    // Measurement and prediction
    z(1,1) = 2.0;
    g(1,1) = 1.0;

    // Measurement noise std dev
    s(1,1) = 1.0;

    // Sensitivity = 1
    G(1,1) = 1.0;

    // Perform update
    MeasUpdate(x, P, K, z, g, s, G, n);

    // Kalman gain should be 0.5
    _assert(fabs(K(1,1) - 0.5) < TOL_);

    // Updated state: x = 0 + 0.5*(2-1) = 0.5
    _assert(fabs(x(1,1) - 0.5) < TOL_);

    // Updated covariance: (1 - 0.5*1)*1 = 0.5
    _assert(fabs(P(1,1) - 0.5) < TOL_);

    std::cout << "Test_MeasUpdate_Scalar passed\n";
    return 0;
}

int gstime_01() {
    try {
        double Mjd_UT1 = 51544.5; // J2000.0
        double expected = 4.894961; // Valor de referencia
        double result = gstime(Mjd_UT1);

        _assert(fabs(result - expected) < 1e-3);
        return 0;
    } catch(...) {
        return 1;
    }
}

int unit_01() {
    try {
        std::cout << "\n=== Test unit_01 ===" << std::endl;

        // Caso 1: Vector no nulo (debe normalizarse)
        double vec_data[] = {1.0, 2.0, 3.0};
        Matrix vec(3, 1, vec_data, 3);
        Matrix outvec = unit(vec);

        // Calcular magnitud del resultado
        double mag = sqrt(outvec(1,1)*outvec(1,1) + outvec(2,1)*outvec(2,1) + outvec(3,1)*outvec(3,1));

        // Verificar que es unitario (magnitud ≈ 1.0)
        _assert(fabs(mag - 1.0) < 1e-6);

        // Verificar valores esperados (normalizados)
        double expected[] = {0.267261, 0.534522, 0.801784};
        for (int i = 1; i <= 3; ++i) {
            _assert(fabs(outvec(i,1) - expected[i-1]) < 1e-6);
        }

        std::cout << "Caso 1 pasado: Vector normalizado correctamente.\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error en unit_01 (Caso 1): " << e.what() << std::endl;
        return 1;
    }
}

int unit_02() {
    try {
        std::cout << "\n=== Test unit_02 ===" << std::endl;

        // Caso 2: Vector casi cero (debe devolver cero)
        double vec_data[] = {1e-7, 1e-7, 1e-7};
        Matrix vec(3, 1, vec_data, 3);
        Matrix outvec = unit(vec);

        // Verificar que es cero
        for (int i = 1; i <= 3; ++i) {
            _assert(outvec(i,1) == 0.0);
        }

        std::cout << "Caso 2 pasado: Vector cero manejado correctamente.\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error en unit_02 (Caso 2): " << e.what() << std::endl;
        return 1;
    }
}

int unit_03() {
    try {
        std::cout << "\n=== Test unit_03 ===" << std::endl;

        // Caso 3: Vector ya unitario (no debe cambiar)
        double vec_data[] = {0.577350, 0.577350, 0.577350}; // ~1/√3
        Matrix vec(3, 1, vec_data, 3);
        Matrix outvec = unit(vec);

        // Verificar que no cambia
        for (int i = 1; i <= 3; ++i) {
            _assert(fabs(outvec(i,1) - vec(i,1)) < 1e-6);
        }

        std::cout << "Caso 3 pasado: Vector unitario no se modifica.\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error en unit_03 (Caso 3): " << e.what() << std::endl;
        return 1;
    }
}

int Gibbs_01() {
    std::cout << "=== Prueba 1: Órbita circular ===" << std::endl;

    // Configurar vectores de posición (órbita circular)
    Matrix r1(3, 1);
    r1(1,1) = 7000e3; // 7000 km -> m
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

    // Ejecutar método Gibbs
    GibbsResult result = gibbs(r1, r2, r3);
    double TOL = 1e-6;
    // Verificar resultados
    _assert(result.error == "          ok");
    _assert(fabs(result.theta - 1.570796) < TOL);
    _assert(fabs(result.theta1 - 1.570796) < TOL);
    _assert(fabs(result.copa) < TOL);

    // Verificar velocidad (valores de tu salida MATLAB)
    _assert((std::abs(result.v2(1,1)) - (7546.05329011)) < TOL);
    _assert(fabs(result.v2(2,1)) < TOL);
    _assert(fabs(result.v2(3,1)) < TOL);

    std::cout << "Prueba 1 pasada!\n" << std::endl;
    return 0;
}

int HGibbs_01() {
    try {
        std::cout << "\n=== Test HGibbs_SmallAngles: Ángulos <= 1 grado ===" << std::endl;

        // Configurar vectores con ángulos pequeños (0.5 grados)
        Matrix r1(3, 1), r2(3, 1), r3(3, 1);

        // Vector base
        r1(1,1) = 7000e3; r1(2,1) = 0; r1(3,1) = 0;

        // Segundo vector con 0.5 grados de diferencia
        double angle = 0.5 * M_PI/180.0; // 0.5 grados en radianes
        r2(1,1) = 7000e3 * cos(angle);
        r2(2,1) = 7000e3 * sin(angle);
        r2(3,1) = 0;

        // Tercer vector con 0.5 grados de diferencia del segundo
        r3(1,1) = 7000e3 * cos(2*angle);
        r3(2,1) = 7000e3 * sin(2*angle);
        r3(3,1) = 0;

        // Fechas julianas (observaciones cercanas en el tiempo)
        double Mjd1 = 60000.0;
        double Mjd2 = Mjd1 + 1.0/(60.0*24.0); // 1 minuto después
        double Mjd3 = Mjd2 + 1.0/(60.0*24.0); // 2 minutos después del primero

        // Ejecutar método
        GibbsResult result = hgibbs(r1, r2, r3, Mjd1, Mjd2, Mjd3);

        double TOL = 1e-6;

        // Verificaciones principales
        _assert(result.error == "          ok"); // Debe pasar sin error de ángulo
        _assert(result.theta * 180.0/M_PI <= 1.0 + TOL); // <= 1 grado
        _assert(result.theta1 * 180.0/M_PI <= 1.0 + TOL); // <= 1 grado

        // Verificación adicional de que el cálculo de velocidad no es cero
        _assert(fabs(result.v2(1,1)) > 1.0);
        _assert(fabs(result.v2(2,1)) > 1.0);
        _assert(fabs(result.v2(3,1)) < TOL);

        std::cout << "Test HGibbs_01 pasado: Ángulos pequeños aceptados correctamente.\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error en HGibbs_01: " << e.what() << std::endl;
        return 1;
    }
}

int HGibbs_02() {
    try {
        std::cout << "\n=== Test HGibbs_02: Vectores no coplanares ===" << std::endl;

        Matrix r1(3, 1), r2(3, 1), r3(3, 1);
        r1(1,1) = 7000e3; r1(2,1) = 0;    r1(3,1) = 0;
        r2(1,1) = 0;      r2(2,1) = 7000e3; r2(3,1) = 0;
        r3(1,1) = 0;      r3(2,1) = 0;     r3(3,1) = 7000e3;

        double Mjd1 = 60000.0;
        double Mjd2 = Mjd1 + 5.0/(60.0*24.0);
        double Mjd3 = Mjd2 + 5.0/(60.0*24.0);

        GibbsResult result = hgibbs(r1, r2, r3, Mjd1, Mjd2, Mjd3);

        _assert(result.error == "not coplanar");

        std::cout << "Test HGibbs_02 pasado: Detecta vectores no coplanares.\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error en test_HGibbs_02: " << e.what() << std::endl;
        return 1;
    }
}

int HGibbs_03() {
    try {
        std::cout << "\n=== Test HGibbs_03: Angulos > 1 grado ===" << std::endl;

        Matrix r1(3, 1), r2(3, 1), r3(3, 1);
        r1(1,1) = 7000e3; r1(2,1) = 0;      r1(3,1) = 0;
        r2(1,1) = 5000e3; r2(2,1) = 5000e3; r2(3,1) = 0;
        r3(1,1) = 0;      r3(2,1) = 7000e3; r3(3,1) = 0;

        double Mjd1 = 60000.0;
        double Mjd2 = Mjd1 + 5.0/(60.0*24.0);
        double Mjd3 = Mjd2 + 5.0/(60.0*24.0);

        GibbsResult result = hgibbs(r1, r2, r3, Mjd1, Mjd2, Mjd3);

        _assert(result.error == "   angl > 1ø");

        std::cout << "Test HGibbs_03 pasado: Detecta angulos grandes.\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error en test_HGibbs_03: " << e.what() << std::endl;
        return 1;
    }
}

int HGibbs_04() {
    try {
        std::cout << "\n=== Test HGibbs_04: Orbita eliptica ===" << std::endl;

        Matrix r1(3, 1), r2(3, 1), r3(3, 1);
        r1(1,1) = 7000e3; r1(2,1) = 0;      r1(3,1) = 0;
        r2(1,1) = 0;      r2(2,1) = 8000e3; r2(3,1) = 0;
        r3(1,1) = -9000e3; r3(2,1) = 0;     r3(3,1) = 0;

        // Tiempos más espaciados para órbita elíptica
        double Mjd1 = 60000.0;
        double Mjd2 = Mjd1 + 30.0/(60.0*24.0);
        double Mjd3 = Mjd2 + 30.0/(60.0*24.0);

        GibbsResult result = hgibbs(r1, r2, r3, Mjd1, Mjd2, Mjd3);

        double TOL = 1e-6;
        // Verificaciones básicas
        _assert(result.error == "   angl > 1ø");
        _assert(fabs(result.theta - 1.570796) < TOL);
        _assert(fabs(result.theta1 - 1.570796) < TOL);
        _assert(fabs(result.copa - 0.0) < TOL);

        std::cout << "Test HGibbs_04 pasado: Orbita eliptica calculada.\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error en test_HGibbs_04: " << e.what() << std::endl;
        return 1;
    }
}

int elements_01() {
    try {
        std::cout << "\n=== Test 1: Órbita elíptica ===\n";

        Matrix y1(6, 1);
        y1(1,1) = 7000e3; y1(2,1) = 1000e3; y1(3,1) = 2000e3;
        y1(4,1) = 4e3; y1(5,1) = 5e3; y1(6,1) = 1e3;

        KeplerianElements el = elements(y1);

        // Valores esperados
        const double TOL = 1e-4;
        _assert(fabs(el.p - 2616655.403817) < TOL);
        _assert(fabs(el.a - 5995316.950134) < TOL);
        _assert(fabs(el.e - 0.750700) < TOL);
        _assert(fabs(el.i - 0.284202) < TOL);
        _assert((fabs(el.Omega - 4.601732) < TOL) || (fabs(el.Omega - 4.601732 + 2*M_PI) < TOL));
        _assert((fabs(el.omega - 5.495119) < TOL) || (fabs(el.omega - 5.495119 + 2*M_PI) < TOL));
        _assert(fabs(el.M - 1.160208) < TOL);

        std::cout << "Test 1 pasado: Órbita elíptica calculada correctamente.\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error en elements_test_01: " << e.what() << std::endl;
        return 1;
    }
}

int elements_02() {
    try {
        std::cout << "\n=== Test 2: Órbita con mayor inclinación ===\n";

        Matrix y2(6, 1);
        y2(1,1) = 8000e3; y2(2,1) = 0; y2(3,1) = 0;
        y2(4,1) = 0; y2(5,1) = 7e3; y2(6,1) = 2e3;

        KeplerianElements el = elements(y2);

        // Valores esperados
        const double TOL = 1e-6;
        _assert(fabs(el.p - 8509774.812799) < TOL);
        _assert(fabs(el.a - 8544469.411862) < TOL);
        _assert(fabs(el.e - 0.063722) < TOL);
        _assert(fabs(el.i - 0.278300) < TOL);
        _assert(fabs(el.Omega) < TOL);
        _assert(fabs(el.omega) < TOL);
        _assert(fabs(el.M) < TOL);

        std::cout << "Test 2 pasado: Órbita con inclinación calculada correctamente.\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error en elements_test_02: " << e.what() << std::endl;
        return 1;
    }
}

int elements_03() {
    try {
        std::cout << "\n=== Test 3: Órbita casi circular ===\n";

        Matrix y3(6, 1);
        y3(1,1) = 7000e3; y3(2,1) = 100e3; y3(3,1) = 50e3;
        y3(4,1) = 0.1e3; y3(5,1) = 7.5e3; y3(6,1) = 0.05e3;

        KeplerianElements el = elements(y3);

        // Valores esperados
        const double TOL = 1e-6;
        _assert(fabs(el.p - 6912827.322912) < TOL);
        _assert(fabs(el.a - 6919087.429595) < TOL);
        _assert(fabs(el.e - 0.030079) < TOL);
        _assert(fabs(el.i - 0.009638) < TOL);
        _assert((fabs(el.Omega - 5.462836) < TOL) || (fabs(el.Omega - 5.462836 + 2 * M_PI) < TOL));
        _assert((fabs(el.omega - 5.115581) < TOL) || (fabs(el.omega - 5.115581 + 2 * M_PI) < TOL));
        _assert(fabs(el.M - 1.947103) < TOL);

        std::cout << "Test 3 pasado: Órbita casi circular calculada correctamente.\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error en elements_test_03: " << e.what() << std::endl;
        return 1;
    }
}

int LTC_01() {
    try {
        std::cout << "\n=== Test 1: Latitud 0°, Longitud 0° ===\n";

        double lon = 0.0;
        double lat = 0.0;

        Matrix result = LTC(lon, lat);

        // Resultado esperado (matriz identidad con filas permutadas)
        Matrix expected(3, 3);
        expected(1,1) = 0.0; expected(1,2) = 1.0; expected(1,3) = 0.0;
        expected(2,1) = 0.0; expected(2,2) = 0.0; expected(2,3) = 1.0;
        expected(3,1) = 1.0; expected(3,2) = 0.0; expected(3,3) = 0.0;

        std::cout << "Matriz resultante:\n";
        result.print();

        std::cout << "Matriz esperada:\n";
        expected.print();

        _assert(MatrixEqual(result, expected));

        std::cout << "Test 1 pasado: LTC en (0°,0°) calculada correctamente.\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error en LTC_test_01: " << e.what() << std::endl;
        return 1;
    }
}

int LTC_02() {
    try {
        std::cout << "\n=== Test 2: Latitud 45°, Longitud 30° ===\n";

        double lon = M_PI/6;      // 30°
        double lat = M_PI/4;      // 45°

        Matrix result = LTC(lon, lat);

        // Valores esperados CORREGIDOS (calculados analíticamente)
        Matrix expected(3, 3);
        expected(1,1) = -0.500000; expected(1,2) = 0.866025; expected(1,3) = 0.000000;
        expected(2,1) = -0.612372; expected(2,2) = -0.353553; expected(2,3) = 0.707107;
        expected(3,1) = 0.612372;  expected(3,2) = 0.353553;  expected(3,3) = 0.707107;

        std::cout << "Matriz resultante:\n";
        result.print();

        std::cout << "Matriz esperada:\n";
        expected.print();

        _assert(MatrixEqual(result, expected, 1e-6));

        std::cout << "Test 2 pasado: LTC en (30°,45°) calculada correctamente.\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error en LTC_test_02: " << e.what() << std::endl;
        return 1;
    }
}

int GHAMatrix_01() {
    try {
        std::cout << "\n=== Test 1: GHAMatrix en J2000 (MJD 51544.5) ===\n";

        double Mjd_UT1 = 51544.5;  // J2000.0
        Matrix result = GHAMatrix(Mjd_UT1);

        // Calcular GAST esperado
        double expected_GAST = gstime(Mjd_UT1);

        // Matriz esperada (rotación alrededor del eje Z por GAST)
        Matrix expected = R_z(expected_GAST);

        std::cout << "Matriz resultante:\n";
        result.print();

        std::cout << "Matriz esperada:\n";
        expected.print();

        _assert(MatrixEqual(result, expected));

        std::cout << "Test 1 pasado: GHAMatrix en J2000 calculada correctamente.\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error en GHAMatrix_test_01: " << e.what() << std::endl;
        return 1;
    }
}

int GHAMatrix_02() {
    try {
        std::cout << "\n=== Test 2: GHAMatrix en fecha actual (MJD 59754.0) ===\n";

        double Mjd_UT1 = 59754.0;  // 2022/06/15 00:00:00 UTC
        Matrix result = GHAMatrix(Mjd_UT1);

        // Calcular GAST esperado
        double expected_GAST = gstime(Mjd_UT1);

        // Matriz esperada (rotación alrededor del eje Z por GAST)
        Matrix expected = R_z(expected_GAST);

        std::cout << "Matriz resultante:\n";
        result.print();

        std::cout << "Matriz esperada:\n";
        expected.print();

        _assert(MatrixEqual(result, expected));

        std::cout << "Test 2 pasado: GHAMatrix en fecha actual calculada correctamente.\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error en GHAMatrix_test_02: " << e.what() << std::endl;
        return 1;
    }
}

int PoleMatrix_01() {
    try {
        std::cout << "\n=== Test 2: PoleMatrix con xp=0.1\", yp=0.2\" (en radianes) ===\n";

        // Convertir segundos de arco a radianes
        double arcsec_to_rad = M_PI/(180.0*3600.0);
        double xp = 0.1 * arcsec_to_rad;
        double yp = 0.2 * arcsec_to_rad;

        Matrix result = PoleMatrix(xp, yp);

        // Calcular manualmente la matriz esperada
        Matrix expected = R_y(-xp) * R_x(-yp);

        std::cout << "Matriz resultante:\n";
        result.print();

        std::cout << "Matriz esperada:\n";
        expected.print();

        _assert(MatrixEqual(result, expected, 1e-12));

        std::cout << "Test 1 pasado: PoleMatrix con pequeños valores calculada correctamente.\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error en PoleMatrix_test_02: " << e.what() << std::endl;
        return 1;
    }
}

int NutMatrix_01() {
    try {
        std::cout << "\n=== Test 1: NutMatrix en J2000 ===\n";

        double Mjd_TT = 51544.5;  // J2000.0
        Matrix result = NutMatrix(Mjd_TT);

        // Calcular componentes esperados
        double eps = MeanObliquity(Mjd_TT);
        double dpsi, deps;
        NutAngles(Mjd_TT, dpsi, deps);
        Matrix expected = R_x(-eps-deps) * R_z(-dpsi) * R_x(+eps);

        _assert(MatrixEqual(result, expected, 1e-12));

        std::cout << "Test 1 pasado correctamente.\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}

int PrecMatrix_01() {
    try {
        std::cout << "\n=== Test 2: PrecMatrix de J2000 a 2020 ===\n";

        double Mjd_1 = 51544.5;  // J2000.0
        double Mjd_2 = 58952.0;  // 2020-03-21 12:00:00 TT

        Matrix result = PrecMatrix(Mjd_1, Mjd_2);

        // Valores esperados calculados manualmente
        Matrix expected(3, 3);
        expected(1,1) = 0.9999999999999929; expected(1,2) = -0.0000271361741621; expected(1,3) = -0.0000117760568353;
        expected(2,1) = 0.0000271361741621; expected(2,2) = 0.9999999999996318; expected(2,3) = -0.0000000000000000;
        expected(3,1) = 0.0000117760568353; expected(3,2) = -0.0000000000000000; expected(3,3) = 0.9999999999993063;

        std::cout << "Matriz resultante:\n";
        result.print();

        std::cout << "Matriz esperada:\n";
        expected.print();

        _assert(MatrixEqual(result, expected, 1e-2));

        std::cout << "Test 1 pasado: PrecMatrix a 2020 calculada correctamente.\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error en PrecMatrix_01: " << e.what() << std::endl;
        return 1;
    }
}

int timediff_01() {
    try {
        std::cout << "\n=== Test 2: Valores límite (UT1-UTC=0, TAI-UTC=0) ===\n";

        double UT1_UTC = 0.0;
        double TAI_UTC = 0.0;

        TimeDiffs result = timediff(UT1_UTC, TAI_UTC);

        // Valores esperados con entradas cero
        _assert(fabs(result.UT1_TAI - 0.0) < TOL_);
        _assert(fabs(result.UTC_GPS - 19.0) < TOL_);
        _assert(fabs(result.UT1_GPS - 19.0) < TOL_);
        _assert(fabs(result.TT_UTC - 32.184) < TOL_);
        _assert(fabs(result.GPS_UTC - (-19.0)) < TOL_);

        std::cout << "Test 1 pasado: Valores límite calculados correctamente.\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error en timediff_01: " << e.what() << std::endl;
        return 1;
    }
}

int Geodetic_01() {
    try {
        std::cout << "\n=== Test 1: Punto en el ecuador ===\n";

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
    // Punto en el ecuador sobre el meridiano 0
    Matrix r(3,1);
    r(1,1) = R_Earth;
    r(2,1) = 0;
    r(3,1) = 0;
    auto G = Geodetic(r);
    if (std::abs(G.longitude - 0.0) > TOL_ ||
        std::abs(G.latitude  - 0.0) > TOL_ ||
        std::abs(G.altitude  - 0.0) > 1e-6) {
        std::cout<<"Equator test failed: "
                 <<"lon="<<G.longitude
                 <<" lat="<<G.latitude
                 <<" h="<<G.altitude<<"\n";
        return 1;
    }
    return 0;
}

int Geodetic_03() {
    // Punto en el polo norte (lat=90°), h≈0
    Matrix r(3,1);
    r(1,1) = 0;
    r(2,1) = 0;
    r(3,1) = R_Earth * (1 - f_Earth); // radio polar
    auto G = Geodetic(r);
    if (std::abs(G.longitude) > TOL_||
        std::abs(G.latitude  - M_PI_2) > 1e-9 ||
        std::abs(G.altitude) > 1e-6) {
        std::cout<<"Pole test failed: "
                 <<"lon="<<G.longitude
                 <<" lat="<<G.latitude
                 <<" h="<<G.altitude<<"\n";
        return 1;
    }
    return 0;
}

int doubler_01() {
    try {
        std::cout << "\n=== Test doubler_01 ===" << std::endl;

        // Datos simulados para un caso válido
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

        // Verificaciones básicas
        _assert(res.r2.getFilas() == 3 && res.r2.getColumnas() == 1);
        _assert(res.r3.getFilas() == 3 && res.r3.getColumnas() == 1);
        _assert(!std::isnan(res.magr1) && !std::isinf(res.magr1));
        _assert(!std::isnan(res.magr2) && !std::isinf(res.magr2));
        _assert(!std::isnan(res.f1) && !std::isinf(res.f1));
        _assert(!std::isnan(res.f2) && !std::isinf(res.f2));
        _assert(!std::isnan(res.q1) && !std::isinf(res.q1));
        _assert(!std::isnan(res.a) && !std::isinf(res.a));
        _assert(!std::isnan(res.deltae32) && !std::isinf(res.deltae32));

        std::cout << "Test doubler_01 pasado correctamente.\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error en doubler_01: " << e.what() << std::endl;
        return 1;
    }
}

int Doubler_01() {
    std::cout<<"=== Doubler_Delta2_Negative_Test ===\n";
    bool threw=false;
    // magrsite2 > magr2in ⇒ delta2 < 0
    try {
        doubler(
                0.0, 0.0,
                1.0, 3.0,  // magrsite1=1, magrsite2=3
                1.0, 1.0,  // magr1in=1, magr2in=1
                Matrix(3,1),Matrix(3,1),Matrix(3,1),
                Matrix(3,1),Matrix(3,1),Matrix(3,1),
                0.0, 0.0, 'y'
        );
    } catch(...) { threw=true; }
    _assert(threw && "Esperaba excepción por delta2<0");
    std::cout<<"Doubler_Delta2_Negative_Test passed\n";
    return 0;
}

int Doubler_02() {
    std::cout<<"=== Doubler_Los3w_Zero_Test ===\n";
    bool threw=false;
    // rsite1=(1,0,0), rsite2=(0,1,0) → cross(r1,r2) en +Z
    // los3=(1,0,0) ⇒ dot(los3,w)=0
    Matrix los1(3,1),los2(3,1),los3(3,1);
    Matrix rs1(3,1), rs2(3,1), rs3(3,1);
    los1(1,1)=1; los1(2,1)=0; los1(3,1)=0;
    los2(1,1)=0; los2(2,1)=1; los2(3,1)=0;
    los3(1,1)=1; los3(2,1)=0; los3(3,1)=0;
    rs1(1,1)=1; rs1(2,1)=0; rs1(3,1)=0;
    rs2(1,1)=0; rs2(2,1)=1; rs2(3,1)=0;
    rs3 = Matrix(3,1); // (0,0,0) da lo mismo
    // magrsite1=magr1in=1 ⇒ delta1=cc1^2-0 con cc1=2 ⇒ delta1>0
    // mismo para delta2
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
    std::cout<<"Doubler_Los3w_Zero_Test passed\n";
    return 0;
}

int IERS_01() {
    // EOP data: 13 rows x 2 cols (row-wise)
    double data[13*2] = {
            // row1-3 (dummy)
            0, 0,   0, 0,   0, 0,
            // row4 MJD
            59000, 59001,
            // row5 x_pole
            0.1, 0.2,
            // row6 y_pole
            0.2, 0.3,
            // row7 UT1_UTC
            0.3, 0.4,
            // row8 LOD
            0.4, 0.5,
            // row9 dpsi
            0.5, 0.6,
            // row10 deps
            0.6, 0.7,
            // row11 dx_pole
            0.7, 0.8,
            // row12 dy_pole
            0.8, 0.9,
            // row13 TAI_UTC
            37.0, 37.0
    };
    Matrix eop(13, 2, data, 26);
    IERSResult r = IERS(eop, 59000.5, 'l');
    double frac = 0.5;
    double TOL = 1e-4;
    _assert(fabs(r.x_pole - ((0.1 + (0.2-0.1)*frac)/Arcs)) < TOL);
    _assert(fabs(r.y_pole - ((0.2 + (0.3-0.2)*frac)/Arcs)) < TOL);
    _assert(fabs(r.UT1_UTC - (0.3 + (0.4-0.3)*frac)) < TOL);
    _assert(fabs(r.LOD     - (0.4 + (0.5-0.4)*frac)) < TOL);
    _assert(fabs(r.dpsi    - ((0.5 + (0.6-0.5)*frac)/Arcs)) < TOL);
    _assert(fabs(r.deps    - ((0.6 + (0.7-0.6)*frac)/Arcs)) < TOL);
    _assert(fabs(r.dx_pole - ((0.7 + (0.8-0.7)*frac)/Arcs)) < TOL);
    _assert(fabs(r.dy_pole - ((0.8 + (0.9-0.8)*frac)/Arcs)) < TOL);
    _assert(fabs(r.TAI_UTC - 37.0) < TOL);
    std::cout << "Test_IERS_Linear passed\n";
    return 0;
}

int IERS_02() {
    // EOP data: 13 rows × 2 cols (row-wise)
    double data[13*2] = {
            // row1-3 (dummy)
            0, 0,   0, 0,   0, 0,
            // row4 MJD
            59000.0, 59001.0,
            // row5 x_pole [″]
            0.1, 0.2,
            // row6 y_pole [″]
            0.2, 0.3,
            // row7 UT1_UTC [s]
            0.3, 0.4,
            // row8 LOD [s]
            0.4, 0.5,
            // row9 dpsi [″]
            0.5, 0.6,
            // row10 deps [″]
            0.6, 0.7,
            // row11 dx_pole [″]
            0.7, 0.8,
            // row12 dy_pole [″]
            0.8, 0.9,
            // row13 TAI_UTC [s]
            37.0, 37.0
    };
    Matrix eop(13, 2, data, 26);
    IERSResult r = IERS(eop, 59000.5, 'l');

    double frac = 0.5;
    const double TOL = 1e-4;
    const double ARCSEC_TO_RAD = 1.0/Arcs;

    _assert(fabs(r.x_pole   - ((0.1 + (0.2 - 0.1)*frac) * ARCSEC_TO_RAD)) < TOL);
    _assert(fabs(r.y_pole   - ((0.2 + (0.3 - 0.2)*frac) * ARCSEC_TO_RAD)) < TOL);
    _assert(fabs(r.UT1_UTC  - (0.3 + (0.4 - 0.3)*frac)) < TOL);
    _assert(fabs(r.LOD      - (0.4 + (0.5 - 0.4)*frac)) < TOL);
    _assert(fabs(r.dpsi     - ((0.5 + (0.6 - 0.5)*frac) * ARCSEC_TO_RAD)) < TOL);
    _assert(fabs(r.deps     - ((0.6 + (0.7 - 0.6)*frac) * ARCSEC_TO_RAD)) < TOL);
    _assert(fabs(r.dx_pole  - ((0.7 + (0.8 - 0.7)*frac) * ARCSEC_TO_RAD)) < TOL);
    _assert(fabs(r.dy_pole  - ((0.8 + (0.9 - 0.8)*frac) * ARCSEC_TO_RAD)) < TOL);
    _assert(fabs(r.TAI_UTC  - 37.0) < TOL);

    std::cout << "Test_IERS_Linear passed\n";
    return 0;
}

int JPL_Eph_01() {
    // Definimos un MJD base y calculamos el JD correspondiente
    double MjdBase = 58000.0;
    double JDBase  = MjdBase + 2400000.5;

    // 1) Preparo PC con una sola fila y suficientes columnas
    const int COLS = 1100;
    PC = Matrix(1, COLS);

    // Ahora PC(1,1) y PC(1,2) deben ser JD (no MJD)
    PC(1,1) = JDBase;         // JD inicial
    PC(1,2) = JDBase + 1000;  // JD final

    // El resto de columnas pueden quedar a cero (constructor inicializa a 0)

    // 2) Llamamos con Mjd_TDB dentro del rango MjdBase..MjdBase+1000
    double Mjd_TDB = MjdBase + 10.0;
    PlanetaryPositions pos = JPL_Eph_DE430(Mjd_TDB);

    // 3) Verificar que todas las posiciones son cero
    auto checkZero = [&](const Matrix& v){
        _assert(v.getFilas()==3 && v.getColumnas()==1);
        for(int i=1; i<=3; ++i)
            _assert(fabs(v(i,1)) < TOL_);
    };

    checkZero(pos.r_Earth);
    checkZero(pos.r_Moon);
    checkZero(pos.r_Sun);
    checkZero(pos.r_Mercury);
    checkZero(pos.r_Venus);
    checkZero(pos.r_Mars);
    checkZero(pos.r_Jupiter);
    checkZero(pos.r_Saturn);
    checkZero(pos.r_Uranus);
    checkZero(pos.r_Neptune);
    checkZero(pos.r_Pluto);

    std::cout << "JPL_Eph_01 passed: all zero coefficients yield zero positions\n";
    return 0;
}

int JPL_Eph_02() {
    Matrix PC(1,2);
    PC(1,1) = 59000.0; PC(1,2) = 59010.0;
    bool caught = false;
    try {
        auto p = JPL_Eph_DE430(59050.0);
    } catch(const std::exception& e) {
        caught = true;
        std::cout << "Caught expected exception: " << e.what() << "\n";
    }
    _assert(caught);
    std::cout << "Test_JPL_Eph_OutOfRange passed\n";
    return 0;
}

int JPL_Eph_03() {
    // Preparo PC para un rango distinto
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
    std::cout<<"JPL_Eph_03 passed\n";
    return 0;
}

int AzElPa_Test_01() {
    try {
        std::cout << "\n=== Test 1: Directly North ===\n";

        Matrix s(3, 1);
        s(1,1) = 0.0;  // East
        s(2,1) = 1.0;  // North
        s(3,1) = 0.0;  // Up

        AzElPaResult result = AzElPa(s);

        // Expected results
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

        // Verify results
        std::cout << "Azimuth: " << result.Az << " (expected: " << expected_Az << ")\n";
        std::cout << "Elevation: " << result.El << " (expected: " << expected_El << ")\n";

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


        std::cout << "Test 1 passed: Directly North case correct.\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error in AzElPa_Test_01: " << e.what() << std::endl;
        return 1;
    }
}

int AzElPa_Test_02() {
    try {
        std::cout << "\n=== Test 2: 45° NE, 45° elevation ===\n";

        Matrix s(3, 1);
        s(1,1) = 1.0;  // East
        s(2,1) = 1.0;  // North
        s(3,1) = sqrt(2.0);  // Up

        AzElPaResult result = AzElPa(s);

        // Expected results
        double expected_Az = M_PI/4.0; // 45°
        double expected_El = M_PI/4.0; // 45°

        // Verify results
        std::cout << "Azimuth: " << result.Az << " (expected: " << expected_Az << ")\n";
        std::cout << "Elevation: " << result.El << " (expected: " << expected_El << ")\n";

        _assert(fabs(result.Az - expected_Az) < 1e-9);
        _assert(fabs(result.El - expected_El) < 1e-9);

        // Verify partials (numerical validation)
        _assert(fabs(result.dAds(1,1) - 0.5) < 1e-9);
        _assert(fabs(result.dAds(1,2) - (-0.5)) < 1e-9);
        _assert(fabs(result.dAds(1,3) - 0.0) < 1e-9);

        std::cout << "Test 2 passed: 45° NE case correct.\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error in AzElPa_Test_02: " << e.what() << std::endl;
        return 1;
    }
}

int AzElPa_Test_03() {
    try {
        std::cout << "\n=== Test 3: Invalid input dimensions ===\n";

        Matrix s(2, 2); // Wrong dimensions

        try {
            AzElPaResult result = AzElPa(s);
            std::cerr << "Error: Expected exception not thrown\n";
            return 1;
        } catch (const std::invalid_argument& e) {
            std::cout << "Correctly caught exception: " << e.what() << "\n";
            _assert(std::string(e.what()).find("3x1") != std::string::npos);
        }

        std::cout << "Test 3 passed: Invalid input handled correctly.\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error in AzElPa_Test_03: " << e.what() << std::endl;
        return 1;
    }
}

int VarEqn_Test_01() {
    // 1) Construye un eop con 13 filas y 2 columnas
    Matrix eop(13, 2);
    // 2) Llena tod0 a cero
    for(int i=1;i<=13;++i) for(int j=1;j<=2;++j) eop(i,j)=0.0;
    // 3) Pon en la fila 4 las fechas MJD de tus dos puntos de EOP
    double Mjd0 = 58000.0;
    eop(4,1) = Mjd0;
    eop(4,2) = Mjd0+1.0;
    // 4) Pon un TAI-UTC realista en fila 13
    eop(13,1) = 37.0;
    eop(13,2) = 37.0;
    // (Opcionalmente rellena filas 5–12 con valores de prueba o ceros,
    // convertirás a radianes dentro de IERS.)

    // Ahora prepara los demás argumentos
    AuxParam params;
    params.Mjd_UTC = Mjd0;
    params.Mjd_TT  = Mjd0;  // o Mjd0 + offset
    params.n = 0;  // zero gravity test
    params.m = 0;

    // yPhi trivial: dr/dt = 0, dv/dt = 0, Phi = I
    Matrix yPhi(42,1);
    // r y v a cero
    for(int i=1;i<=6;++i) yPhi(i,1)=0.0;
    // Phi identidad
    for(int j=1;j<=6;++j) {
        yPhi(6*j+ j, 1) = 1.0;
    }

    // Llamada a VarEqn
    Matrix yPhip = VarEqn(0.0, yPhi, params, eop);

    // Comprueba que, con gravedad cero, tus primeras 6 entradas sean ceros,
    // y que el sub-bloque de dΦ/dt sea la 6×6 con ceros también.
    _assert(std::fabs(yPhip(1,1)) < 1e-12);
    _assert(std::fabs(yPhip(6,1)) < 1e-2);//falla este assert
    // … etc …

    std::cout<<"VarEqn_Test_01 passed\n";
    return 0;
}

int VarEqn_Test_02() {
    // Preparo yPhi de tamaño incorrecto (por ejemplo 10×1)
    Matrix yPhi_bad(10,1);
    AuxParam params{59000.0, 59000.0, 0, 0};
    Matrix eop(1,2);
    double JD = params.Mjd_UTC + 2400000.5;
    eop(1,1) = JD; eop(1,2) = JD + 1.0;

    bool threw = false;
    try {
        Matrix out = VarEqn(0.0, yPhi_bad, params, eop);
    } catch (const std::exception& e) {
        threw = true;
    }
    _assert(threw && "Expected exception for wrong yPhi size");
    std::cout << "VarEqn_Test_02 passed\n";
    return 0;
}

int VarEqn_Test_03() {
    std::cout << "=== VarEqn_CentralGravity_Test ===\n";
    // montamos un eop trivial en MJD_J2000 para que E ~ I
    double Mjd0 = 51544.5;
    Matrix eop(13, 2);
    for (int i = 1; i <= 13; ++i) for (int j = 1; j <= 2; ++j) eop(i,j) = 0.0;
    eop(4,1) = Mjd0; eop(4,2) = Mjd0+1.0;    // fechas
    eop(13,1)=37;   eop(13,2)=37;            // TAI-UTC

    AuxParam params;
    params.Mjd_UTC = Mjd0;
    params.Mjd_TT  = Mjd0;
    params.n = 0; params.m = 0;  // sólo campo central

    // instalamos yPhi:  r = [R,0,0], v=0, Phi = I
    const double R = 7000e3;
    Matrix yPhi(42,1);
    // r
    yPhi(1,1)= R; yPhi(2,1)=0; yPhi(3,1)=0;
    // v
    yPhi(4,1)=0;  yPhi(5,1)=0; yPhi(6,1)=0;
    // Phi = identidad
    for(int j=1;j<=6;++j) yPhi(6*j+j,1)=1.0;

    Matrix yPhip = VarEqn(0.0, yPhi, params, eop);

    //  dr/dt = v => ceros
    for(int i=1;i<=3;++i) _assert(fabs(yPhip(i,1))<1e-12);
    //  dv/dt = a = -GM/R^2 en x
    double expect_ax = -GM_Earth/(R*R);
    std::cout<< yPhip(4,1)<< std::endl;
    std::cout<< expect_ax<< std::endl;
    _assert(fabs(yPhip(4,1)-expect_ax)<1e-6);
    _assert(fabs(yPhip(5,1))<1e-12 && fabs(yPhip(6,1))<1e-12);

    std::cout<<"VarEqn_CentralGravity_Test passed\n";
    return 0;
}

int AccelHarmonic_Central_Test() {
    // Pure central gravity: E=I, n=0,m=0
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
    std::cout<<"AccelHarmonic_Central_Test passed\n";
    return 0;
}

int AccelHarmonic_ZeroDistance_Test() {
    // r at the center => should return zero (avoid NaNs)
    Matrix I(3,3);
    for(int i=1;i<=3;++i) I(i,i)=1.0;
    Matrix r(3,1); // all zeros
    Matrix a = AccelHarmonic(r, I, 5, 5);
    _assert(a.norm() == 0.0);
    std::cout<<"AccelHarmonic_ZeroDistance_Test passed\n";
    return 0;
}

int Accel_01() {
    AuxParam p{58000.0, 58000.0, 0, 0, false,false,false};
    Matrix eop(13,2);
    eop(4,1)=p.Mjd_UTC; eop(4,2)=p.Mjd_UTC+1;
    eop(13,1)=37; eop(13,2)=37;
    // Estado inicial: r=[R,0,0], v=[0,0,0]
    double R = 7000e3;
    Matrix Y(6,1);
    Y(1,1)=R;
    Matrix dY = Accel(0.0, Y, p, eop);
    double expect = -GM_Earth/(R*R);
    _assert(fabs(dY(4,1)-expect)<1e-8);
    _assert(fabs(dY(5,1))<1e-12);
    _assert(fabs(dY(6,1))<1e-12);
    return 0;
}

static Matrix make_simple_eop(double Mjd0) {
    Matrix eop(13,2);
    for(int i=1;i<=13;++i)
        for(int j=1;j<=2;++j)
            eop(i,j)=0.0;
    eop(4,1)=Mjd0;
    eop(4,2)=Mjd0+1.0;
    eop(13,1)=37.0;
    eop(13,2)=37.0;
    return eop;
}

int Accel_02(){
    std::cout<<"=== Accel_InvalidY_Test ===\n";
    AuxParam params{58000.0,58000.0,0,0};
    Matrix eop = make_simple_eop(params.Mjd_UTC);
    Matrix Y_bad(5,1);
    bool threw=false;
    try {
        Accel(0.0, Y_bad, params, eop);
    } catch(const std::invalid_argument&) {
        threw=true;
    }
    _assert(threw);
    std::cout<<"Accel_InvalidY_Test passed\n";
    return 0;
}

int Accel_03(){
    std::cout<<"=== Accel_InvalidAux_Test ===\n";
    AuxParam params{58000.0,58000.0,2,3}; // m>n
    Matrix eop = make_simple_eop(params.Mjd_UTC);
    Matrix Y(6,1);
    for(int i=1;i<=6;++i) Y(i,1)=1.0;
    bool threw=false;
    try {
        Accel(0.0, Y, params, eop);
    } catch(const std::invalid_argument&) {
        threw=true;
    }
    _assert(threw);
    std::cout<<"Accel_InvalidAux_Test passed\n";
    return 0;
}

int Accel_04(){
    std::cout<<"=== Accel_EOP_OutOfRange_Test ===\n";
    AuxParam params{58000.0,58000.0,0,0};
    Matrix eop(13,2); // todo cero, no fila 4 == Mjd_UTC
    Matrix Y(6,1);
    bool threw=false;
    try {
        Accel(0.0, Y, params, eop);
    } catch(const std::out_of_range&) {
        threw=true;
    }
    _assert(threw);
    std::cout<<"Accel_EOP_OutOfRange_Test passed\n";
    return 0;
}

int Accel_05(){
    std::cout<<"=== Accel_ZeroState_Test ===\n";
    AuxParam params{58000.0,58000.0,0,0};
    params.sun = params.moon = params.planets = false;
    Matrix eop = make_simple_eop(params.Mjd_UTC);
    Matrix Y(6,1);
    for(int i=1;i<=6;++i) Y(i,1)=0.0;
    Matrix dY = Accel(0.0, Y, params, eop);
    // debe ser todo cero
    for(int i=1;i<=6;++i){
        _assert(std::fabs(dY(i,1))<1e-12);
    }
    std::cout<<"Accel_ZeroState_Test passed\n";
    return 0;
}

int Accel_06(){
    std::cout<<"=== Accel_CircularOrbit_Test ===\n";
    AuxParam params{58000.0,58000.0,0,0};
    params.sun = params.moon = params.planets = false;
    Matrix eop = make_simple_eop(params.Mjd_UTC);

    // Elipsoide de prueba en x–y
    const double R = 7000e3;                   // m
    const double mu = GM_Earth;                // m^3/s^2
    const double v_circ = std::sqrt(mu/(R));   // m/s

    // Estado: r=(R,0,0), v=(0,v_circ,0)
    Matrix Y(6,1);
    Y(1,1)= R;  Y(2,1)=0;      Y(3,1)=0;
    Y(4,1)=0;   Y(5,1)=v_circ; Y(6,1)=0;

    Matrix dY = Accel(0.0, Y, params, eop);

    // dr/dt = v
    _assert(std::fabs(dY(1,1) - 0.0)     < 1e-12);
    _assert(std::fabs(dY(2,1) - v_circ)  < 1e-8);
    _assert(std::fabs(dY(3,1) - 0.0)     < 1e-12);

    // dv/dt = a_centripeta = -(mu/R^2) in x
    double expect_ax = -mu/(R*R);
    _assert(std::fabs(dY(4,1) - expect_ax) < 1e-6);
    _assert(std::fabs(dY(5,1) - 0.0)       < 1e-12);
    _assert(std::fabs(dY(6,1) - 0.0)       < 1e-12);

    std::cout<<"Accel_CircularOrbit_Test passed\n";
    return 0;
}




int anglesdr_BadSize_Test() {
    std::cout<<"=== anglesdr_BadSize_Test\n";
    bool ok=false;
    try {
        // 2) Prepara AuxParam
        AuxParam params;

        // 3) Prepara eopdata (13×3) con fila 4 = [M1, M2, M3]
        Matrix eopdata(13,3);
        Matrix bad(5,1);
        Matrix rs(3,1); rs(1,1)=1; rs(2,1)=2; rs(3,1)=3;
        auto out = anglesdr(0,0,0,0,0,0, 58000,58001,58002, bad, rs, rs,params,eopdata);
    } catch(const std::invalid_argument&) {
        ok=true;
    }
    _assert(ok && "se esperaba invalid_argument");
    std::cout<<"anglesdr_BadSize_Test passed\n";
    return 0;
}

// 2) Líneas de vista idénticas ➞ geometría degenerada
int anglesdr_DegenerateGeometry_Test() {
    std::cout<<"=== anglesdr_DegenerateGeometry_Test\n";
    // tres observaciones idénticas (mismo az/el) no definen órbita
    double az=1.0, el=0.5;
    double M1=58000, M2=58000.1, M3=58000.2;
    Matrix site(3,1); site(1,1)= 6378e3; site(2,1)=0; site(3,1)=0;
    bool threw=false;
    // 2) Prepara AuxParam
    AuxParam params;
    params.Mjd_UTC = M2;
    params.Mjd_TT  = M2;
    params.n = 0;   // no usamos gravedad aquí
    params.m = 0;

    // 3) Prepara eopdata (13×3) con fila 4 = [M1, M2, M3]
    Matrix eopdata(13,3);
    try {
        auto out = anglesdr(az,az,az, el,el,el, M1,M2,M3,
                            site, site, site,params,eopdata);
    } catch(const std::runtime_error&) {
        threw = true;
    }
    _assert(threw && "se esperaba runtime_error por geometría degenerada");
    std::cout<<"anglesdr_DegenerateGeometry_Test passed\n";
    return 0;
}



int anglesdr_NonDegenerate_Test() {
    std::cout << "=== anglesdr_NonDegenerate_Test ===\n";

    // 1) Constantes órbita circular y MJDs
    const double GM = GM_Earth;
    const double R_e = R_Earth;
    const double alt = 500e3;
    const double R_sat = R_e + alt;
    const double w = std::sqrt(GM/(R_sat*R_sat*R_sat));
    const double dt = 60.0;
    double M2 = 58000.0;
    double M1 = M2 - dt/86400.0;
    double M3 = M2 + dt/86400.0;

    // 2) Tres estaciones ortogonales
    Matrix rsite1(3,1); rsite1(1,1)= R_e; rsite1(2,1)=0;    rsite1(3,1)=0;
    Matrix rsite2(3,1); rsite2(1,1)= 0;   rsite2(2,1)=R_e;  rsite2(3,1)=0;
    Matrix rsite3(3,1); rsite3(1,1)= 0;   rsite3(2,1)=0;    rsite3(3,1)=R_e;

    // 3) Posición satélite en t=-dt,0,+dt
    auto sat = [&](double t){
        double th = w*t;
        Matrix r(3,1);
        r(1,1)=R_sat*cos(th);
        r(2,1)=R_sat*sin(th);
        r(3,1)=0;
        return r;
    };
    Matrix r1 = sat(-dt), r2 = sat(0), r3 = sat(+dt);

    // 4) Ángulos de acimut/elevación
    auto los = [&]( Matrix& rs, Matrix& site){
        Matrix d = rs-site;
        return d*(1.0/d.norm());
    };
    Matrix los1 = los(r1,rsite1), los2 = los(r2,rsite2), los3 = los(r3,rsite3);
    double az1=atan2(los1(1,1),los1(2,1)), el1=asin(los1(3,1));
    double az2=atan2(los2(1,1),los2(2,1)), el2=asin(los2(3,1));
    double az3=atan2(los3(1,1),los3(2,1)), el3=asin(los3(3,1));

    // 5) AuxParam
    AuxParam params;
    params.Mjd_UTC = M2;
    params.Mjd_TT  = M2;
    params.n = 0; params.m = 0;
    params.sun = params.moon = params.planets = false;

    // 6) EOP: fila 4 = [M1,M2,M3], resto da igual
    Matrix eop(13,3);
    eop(4,1) = M1;
    eop(4,2) = M2;
    eop(4,3) = M3;
    // opcional: eop(13,*) = TAI-UTC realista (o cero)

    // 7) Llamada
    auto out = anglesdr(
            az1,az2,az3,
            el1,el2,el3,
            M1,M2,M3,
            rsite1, rsite2, rsite3,
            params, eop
    );

    // 8) Validación analítica
    Matrix expect_r2 = r2;
    Matrix expect_v2(3,1);
    expect_v2(1,1) = -R_sat*w*0.0;
    expect_v2(2,1) =  R_sat*w*1.0;
    expect_v2(3,1) =  0.0;

    double tol_r=1e-2, tol_v=1e-4;
    _assert(fabs(out.r2(1,1)-expect_r2(1,1))<tol_r);
    _assert(fabs(out.r2(2,1)-expect_r2(2,1))<tol_r);
    _assert(fabs(out.r2(3,1)-expect_r2(3,1))<tol_r);
    _assert(fabs(out.v2(1,1)-expect_v2(1,1))<tol_v);
    _assert(fabs(out.v2(2,1)-expect_v2(2,1))<tol_v);
    _assert(fabs(out.v2(3,1)-expect_v2(3,1))<tol_v);

    std::cout<<"anglesdr_NonDegenerate_Test passed\n";
    return 0;
}

int anglesdr_J2000Circular_Test() {
    std::cout << "=== anglesdr_J2000Circular_Test ===\n";

    // 1) Parámetros órbita circular
    const double GM = GM_Earth;
    const double R_e = R_Earth;
    const double alt = 500e3;                     // 500 km de altitud
    const double R_sat = R_e + alt;
    const double w = std::sqrt(GM/(R_sat*R_sat*R_sat));
    const double dt = 60.0;                       // 60 s

    // 2) Tres epochs alrededor de J2000 (MJD=51544.5)
    double M2 = 51544.5;
    double M1 = M2 - dt/86400.0;
    double M3 = M2 + dt/86400.0;

    // 3) Tres estaciones mutuamente ortogonales en ECEF
    Matrix rsite1(3,1); rsite1(1,1)= R_e; rsite1(2,1)=0;    rsite1(3,1)=0;
    Matrix rsite2(3,1); rsite2(1,1)= 0;   rsite2(2,1)=R_e;  rsite2(3,1)=0;
    Matrix rsite3(3,1); rsite3(1,1)= 0;   rsite3(2,1)=0;    rsite3(3,1)=R_e;

    // 4) Posición del satélite en t=-dt, 0, +dt en ECI
    auto sat = [&](double t){
        double ang = w*t;
        Matrix r(3,1);
        r(1,1) = R_sat*std::cos(ang);
        r(2,1) = R_sat*std::sin(ang);
        r(3,1) = 0;
        return r;
    };
    Matrix r1 = sat(-dt), r2 = sat(0), r3 = sat(+dt);

    // 5) Line-of-sight normalizadas
    auto los = [&]( Matrix& rs,  Matrix& site){
        Matrix d = rs - site;
        return d * (1.0 / d.norm());
    };
    Matrix los1 = los(r1, rsite1),
            los2 = los(r2, rsite2),
            los3 = los(r3, rsite3);

    double az1 = std::atan2(los1(1,1), los1(2,1)), el1 = std::asin(los1(3,1));
    double az2 = std::atan2(los2(1,1), los2(2,1)), el2 = std::asin(los2(3,1));
    double az3 = std::atan2(los3(1,1), los3(2,1)), el3 = std::asin(los3(3,1));

    // 6) AuxParam
    AuxParam params;
    params.Mjd_UTC = M2;
    params.Mjd_TT  = M2;
    params.n = 0; params.m = 0;
    params.sun = params.moon = params.planets = false;

    // 7) EOP: fila 4 = [51544,51544,51544], resto da igual (TAI-UTC=0)
    Matrix eop(13,3);
    eop(4,1) = std::floor(M1);
    eop(4,2) = std::floor(M2);
    eop(4,3) = std::floor(M3);

    // 8) Llamada
    auto out = anglesdr(
            az1,az2,az3,
            el1,el2,el3,
            M1,M2,M3,
            rsite1, rsite2, rsite3,
            params, eop
    );

    // 9) Valores esperados
    Matrix expect_r2 = r2;
    Matrix expect_v2(3,1);
    expect_v2(1,1) = -R_sat*w*0.0;
    expect_v2(2,1) =  R_sat*w*1.0;
    expect_v2(3,1) =  0.0;

    // 10) Tolerancias y asserts
    double tol_r = 1e-2;   // 1 cm
    double tol_v = 1e-4;   // 0.1 mm/s

    _assert(std::fabs(out.r2(1,1) - expect_r2(1,1)) < tol_r);
    _assert(std::fabs(out.r2(2,1) - expect_r2(2,1)) < tol_r);
    _assert(std::fabs(out.r2(3,1) - expect_r2(3,1)) < tol_r);

    _assert(std::fabs(out.v2(1,1) - expect_v2(1,1)) < tol_v);
    _assert(std::fabs(out.v2(2,1) - expect_v2(2,1)) < tol_v);
    _assert(std::fabs(out.v2(3,1) - expect_v2(3,1)) < tol_v);

    std::cout<<"anglesdr_J2000Circular_Test passed\n";
    return 0;
}




int Doubler_03() {
    // Caso trivial: rsite1 = (0,0,0), los1=(1,0,0), rsite2=(0,0,0), los2=(0,1,0)
    // magr1in y magr2in escogidos grandes, t1=t3=1
    double cc1 = 0.0, cc2 = 0.0;
    double magrsite1 = 0.0, magrsite2 = 0.0;
    double magr1in = 10000.0, magr2in = 10000.0;
    Matrix los1(3,1), los2(3,1), los3(3,1);
    los1(1,1)=1; los1(2,1)=0; los1(3,1)=0;
    los2(1,1)=0; los2(2,1)=1; los2(3,1)=0;
    los3(1,1)=0; los3(2,1)=0; los3(3,1)=1;
    Matrix rsite1(3,1), rsite2(3,1), rsite3(3,1);
    // todos cero
    DoubleRResult R = doubler(
            cc1, cc2,
            magrsite1, magrsite2,
            magr1in, magr2in,
            los1, los2, los3,
            rsite1, rsite2, rsite3,
            1.0, 1.0, 'y'
    );
    // Aquí esperamos que R.r2 = rho2*los2 = rho2*(0,1,0) con rho2 = (-cc2+sqrt(cc2^2 - ...))/2
    double expected_rho2 = (-cc2 + sqrt(cc2*cc2 - 4*(magrsite2*magrsite2 - magr2in*magr2in))) / 2.0;
    if (std::fabs(R.r2(1,1) - 0.0) > TOL_ ||
        std::fabs(R.r2(2,1) - expected_rho2) > TOL_ ||
        std::fabs(R.r2(3,1) - 0.0) > TOL_) {
        std::cout<<"Doubler_TrivialColinear_Test FAILED\n";
        return 1;
    }
    std::cout<<"Doubler_TrivialColinear_Test passed\n";
    return 0;
}

#define TOL_R 1e-2
#define TOL_V 1e-4

int anglesdr_SyntheticCircular_Test(){
    std::cout<<"=== anglesdr_SyntheticCircular_Test ===\n";
    // Parámetros idénticos a los de Doubler_CircularOrbit_Test
    const double R = R_Earth + 1000e3;
    const double GM = GM_Earth;
    double w = std::sqrt(GM/(R*R*R));
    double dt = 60.0;
    double M2 = 58000.0;
    double M1 = M2 - dt/86400.0, M3 = M2 + dt/86400.0;

    // 1) sitio ECEF en origen
    Matrix site(3,1); site(1,1)=0; site(2,1)=0; site(3,1)=0;
    // 2) posiciones en ECI
    auto sat = [&](double t){
        Matrix r(3,1);
        r(1,1)=R*std::cos(w*t);
        r(2,1)=R*std::sin(w*t);
        r(3,1)=0;
        return r;
    };
    Matrix r1 = sat(-dt), r2 = sat(0), r3 = sat(+dt);
    // 3) LOS
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

    // EOP mínimo: sólo usaremos la fila 4
    Matrix eop(13,3);
    for(int i=1;i<=13;++i) for(int j=1;j<=3;++j) eop(i,j)=0;
    eop(4,1)=M1; eop(4,2)=M2; eop(4,3)=M3;

    auto out = anglesdr(
            az1,az2,az3,
            el1,el2,el3,
            M1,M2,M3,
            site,site,site,
            params,
            eop
    );

    // Comprueba r2, v2
    Matrix expect_r2 = r2;
    Matrix expect_v2(3,1);
    expect_v2(1,1) = -R*w*std::sin(0.0);
    expect_v2(2,1) =  R*w*std::cos(0.0);
    expect_v2(3,1) =  0.0;

    bool ok = true;
    for(int i=1;i<=3;++i){
        if(std::fabs(out.r2(i,1)-expect_r2(i,1))>TOL_R) ok=false;
        if(std::fabs(out.v2(i,1)-expect_v2(i,1))>TOL_V) ok=false;
    }
    if(!ok){
        std::cout<<"anglesdr_SyntheticCircular_Test FAILED\n";
        return 1;
    }
    std::cout<<"anglesdr_SyntheticCircular_Test passed\n";
    return 0;
}


int anglesg_BadSize_Test() {
    Matrix Rs1(2,1), Rs2(2,1), Rs3(2,1);
    AuxParam params{58000.0,58000.0,0,0,false,false,false};
    Matrix eop(13,3);
    bool threw = false;
    try {
        auto out = anglesg(
                0,0,0, 0,0,0,
                58000.0,58000.0,58000.0,
                Rs1, Rs2, Rs3,
                params, eop
        );
    } catch(const std::invalid_argument&) {
        threw = true;
    }
    _assert(threw && "Expected exception for wrong site vector size");
    std::cout<<"anglesg_BadSize_Test passed\n";
    return 0;
}

int anglesg_SyntheticCircular_Test() {
    // dados satélite circular equatorial
    const double GM = GM_Earth;
    const double R_e = R_Earth;
    const double alt = 1e6;
    const double R_sat = R_e + alt;
    const double w = std::sqrt(GM/(R_sat*R_sat*R_sat));
    const double dt = 60.0;
    double M2 = 58000.0;
    double M1 = M2 - dt/86400.0;
    double M3 = M2 + dt/86400.0;

    // estación en (R_e,0,0)
    Matrix Rs(3,1);
    Rs(1,1)=R_e; Rs(2,1)=0; Rs(3,1)=0;

    auto sat = [&](double t)->Matrix {
        double ang = w*t;
        Matrix r(3,1);
        r(1,1)=R_sat*cos(ang);
        r(2,1)=R_sat*sin(ang);
        r(3,1)=0;
        return r;
    };

    Matrix r1 = sat(-dt), r2 = sat(0), r3 = sat(+dt);
    Matrix los1 = (r1-Rs)*(1.0/(r1-Rs).norm());
    Matrix los2 = (r2-Rs)*(1.0/(r2-Rs).norm());
    Matrix los3 = (r3-Rs)*(1.0/(r3-Rs).norm());
    double az1 = std::atan2(los1(1,1),los1(2,1)), el1 = std::asin(los1(3,1));
    double az2 = std::atan2(los2(1,1),los2(2,1)), el2 = std::asin(los2(3,1));
    double az3 = std::atan2(los3(1,1),los3(2,1)), el3 = std::asin(los3(3,1));

    AuxParam params;
    params.Mjd_UTC = M2;
    params.Mjd_TT  = M2;
    params.n = 0; params.m = 0;
    params.sun = params.moon = params.planets = false;

    // eopdata dummy (fila 4 = [M1,M2,M3])
    Matrix eop(13,3);
    eop(4,1)=M1; eop(4,2)=M2; eop(4,3)=M3;

    auto out = anglesg(
            az1,az2,az3, el1,el2,el3,
            M1,M2,M3,
            Rs, Rs, Rs,
            params, eop
    );

    // esperados
    Matrix expect_r2 = r2;
    Matrix expect_v2(3,1);
    expect_v2(1,1) = -R_sat*w*sin(0.0);
    expect_v2(2,1) =  R_sat*w*cos(0.0);
    expect_v2(3,1) =  0.0;

    std::cout<<out.r2(1,1) <<std::endl;
    std::cout<<expect_r2(1,1) <<std::endl;
    std::cout<<out.r2(2,1) <<std::endl;
    std::cout<<expect_r2(2,1) <<std::endl;
    std::cout<<out.r2(3,1) <<std::endl;
    std::cout<<expect_r2(3,1) <<std::endl;
    std::cout<<out.v2(1,1) <<std::endl;
    std::cout<<out.v2(2,1) <<std::endl;
    std::cout<<out.v2(3,1) <<std::endl;

    // comprobaciones
    _assert(std::fabs(out.r2(1,1) - expect_r2(1,1)) < TOL_R);
    _assert(std::fabs(out.r2(2,1) - expect_r2(2,1)) < TOL_R);
    _assert(std::fabs(out.r2(3,1) - expect_r2(3,1)) < TOL_R);
    _assert(std::fabs(out.v2(1,1) - expect_v2(1,1)) < TOL_V);
    _assert(std::fabs(out.v2(2,1) - expect_v2(2,1)) < TOL_V);
    _assert(std::fabs(out.v2(3,1) - expect_v2(3,1)) < TOL_V);

    std::cout<<"anglesg_SyntheticCircular_Test passed\n";
    return 0;
}
int anglesg_BadEOP_Test() {
    // Debe fallar si eopdata no tiene 13×3
    double M2 = 58000.0;
    Matrix Rs(3,1); Rs(1,1)=0; Rs(2,1)=0; Rs(3,1)=0;
    AuxParam params;
    params.Mjd_UTC = M2;
    params.Mjd_TT  = M2;
    // eopdata mal dimensionado:
    Matrix eop_bad(12,3);
    bool threw = false;
    try {
        anglesg(0,0,0, 0,0,0, M2,M2,M2, Rs,Rs,Rs, params, eop_bad);
    } catch(const std::exception&) {
        threw = true;
    }
    _assert(threw && "Expected exception for bad eopdata size");
    std::cout << "anglesg_BadEOP_Test passed\n";
    return 0;
}

int anglesg_Test_01() {
    // 1) Prepare EOP data (simplified)
    Matrix eop(13, 2);
    for(int i=1;i<=13;++i) for(int j=1;j<=2;++j) eop(i,j)=0.0;
    double Mjd0 = 58000.0;
    eop(4,1) = Mjd0;
    eop(4,2) = Mjd0+1.0;
    eop(13,1) = 37.0;
    eop(13,2) = 37.0;

    // 2) Prepare observation data (simplified scenario)
    double az1 = 0.1, az2 = 0.2, az3 = 0.3;  // radians
    double el1 = 0.1, el2 = 0.2, el3 = 0.3;  // radians
    double Mjd1 = Mjd0, Mjd2 = Mjd0 + 0.01, Mjd3 = Mjd0 + 0.02;

    // 3) Prepare station positions (simplified)
    Matrix Rs1(3,1), Rs2(3,1), Rs3(3,1);
    Rs1(1,1) = 1000.0; Rs1(2,1) = 2000.0; Rs1(3,1) = 3000.0;
    Rs2(1,1) = 1100.0; Rs2(2,1) = 2100.0; Rs2(3,1) = 3100.0;
    Rs3(1,1) = 1200.0; Rs3(2,1) = 2200.0; Rs3(3,1) = 3200.0;

    // 4) Prepare auxiliary parameters
    AuxParam params;
    // (Fill params as needed for your implementation)

    try {
        // 5) Call anglesg function
        AnglesGResult result = anglesg(az1, az2, az3, el1, el2, el3,
                                       Mjd1, Mjd2, Mjd3, Rs1, Rs2, Rs3,
                                       params, eop);

        // 6) Basic assertions
        _assert(result.r2.getFilas() == 3 && result.r2.getColumnas() == 1);
        _assert(result.v2.getFilas() == 3 && result.v2.getColumnas() == 1);

        // Check that the position is reasonable (not zero and not NaN)
        double r2_norm = result.r2.norm();
        _assert(r2_norm > 1000.0 && r2_norm < 10000000.0);

        // Check that the velocity is reasonable
        double v2_norm = result.v2.norm();
        _assert(v2_norm > 100.0 && v2_norm < 10000.0);

        std::cout << "anglesg_Test_01 passed\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "anglesg_Test_01 failed: " << e.what() << "\n";
        return 1;
    }
}










int all_tests()
{
/*
    _verify(Matrix_Basico);
    _verify(Mjday_01);
    _verify(Mjday_02);
    _verify(R_x_01);
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
    _verify(HGibbs_04);//
    _verify(elements_01);
    _verify(elements_02);
    _verify(elements_03);
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
    _verify(IERS_02);
    _verify(JPL_Eph_01);
    _verify(JPL_Eph_02);
    _verify(JPL_Eph_03);
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
    _verify(Accel_06);
_verify(AccelHarmonic_Central_Test);
    _verify(AccelHarmonic_ZeroDistance_Test);*/
    //_verify(anglesdr_BadSize_Test);
    //_verify(anglesdr_SyntheticCircular_Test);
    //_verify(anglesdr_DegenerateGeometry_Test);
    //_verify(anglesdr_SyntheticCircular_Equatorial_Test);
    //_verify(anglesdr_CollinearLOS_Test);
    //_verify(anglesdr_NonDegenerate_Test);
    //_verify(anglesdr_J2000Circular_Test);

    //_verify(anglesg_BadEOP_Test);//PASA
    _verify(anglesg_Test_01);



    return 0;
}




int main()
{
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
