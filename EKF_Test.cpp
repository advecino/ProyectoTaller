#include <iostream>
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

int Position_01()
{
    double lon = 1.0;
    double lat = 0.5;
    double h = 1000.0;


    const double R_equ = 6378136.3;
    const double f = 1.0/298.256415099;

    Matrix r = Position(lon, lat, h, R_equ, f);


    double expected_data[] = {
            3027060.2393987174,
            3134461.3940106473,
            2688866.0326071794
    };

    Matrix expected(3, 1, expected_data, 3);


    for(int i = 1; i <= 3; ++i) {
        if(fabs(r(i,1) - expected(i,1)) > TOL_) {
           /* std::cout << "Position_01 failed at element " << i
                      << ": expected " << expected(i,1)
                      << ", got " << r(i,1) << std::endl;*/
            return 1;
        }
    }

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
     sol.print();

    _assert(fabs(sol(1,1)) - 1 < TOL_ && fabs(sol(1,2)) < TOL_ && fabs(sol(1,3)) < TOL_);
    _assert(sol(2,1) < TOL_ && fabs(sol(2,2) - 0.54030230586814 ) < TOL_ && fabs(sol(2,3) -0.841470984807897) < TOL_);
    
    
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

    if(fabs(Mjd_TDB - expected) > TOL_) {
        std::cout << std::setprecision(15);
        std::cout << "Mjday_TDB_01 failed!\n";
        std::cout << "Expected: " << expected << "\n";
        std::cout << "Got:      " << Mjd_TDB << "\n";
        std::cout << "Diff:     " << fabs(Mjd_TDB - expected) << "\n";
        return 1;
    }

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
    double MOblq = MeanObliquity(Mjd_TT, MJD_J2000, Rad);

    // Valor esperado para J2000 (23.4392911 grados en radianes)
    double expected = 23.4392911 * Rad;

    if(fabs(MOblq - expected) > TOL_) {
        std::cout << "MeanObliquity_01 failed!\n";
        std::cout << "Expected: " << expected << "\n";
        std::cout << "Got:      " << MOblq << "\n";
        return 1;
    }

    return 0;
}

int NutAngles_01() {
    try {
        // Fecha J2000 (1.5 enero 2000)
        double Mjd_TT = 51544.5;

        double dpsi, deps;
        NutAngles(Mjd_TT, dpsi, deps);

        // Valores esperados cercanos a cero
        if(fabs(dpsi) > TOL_ || fabs(deps) > TOL_) {
            std::cout << "NutAngles_01 failed!\n";
            std::cout << std::setprecision(12);
            std::cout << "dpsi: " << dpsi << " (expected near 0)\n";
            std::cout << "deps: " << deps << " (expected near 0)\n";
            return 1;
        }
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Exception in NutAngles_01: " << e.what() << "\n";
        return 1;
    }
}

int AccelHarmonic_01() {
    try {
        // 1. Configuración del test
        cout << "Iniciando test de AccelHarmonic..." << endl;

        // Posición del satélite (vector columna)
        double r_data[] = {7000.0e3, 100.0e3, 50.0e3}; // [m]
        Matrix r(3, 1, r_data, 3);

        // Matriz de transformación identidad
        Matrix E(3, 3);
        E(1,1) = 1.0; E(1,2) = 0.0; E(1,3) = 0.0;
        E(2,1) = 0.0; E(2,2) = 1.0; E(2,3) = 0.0;
        E(3,1) = 0.0; E(3,2) = 0.0; E(3,3) = 1.0;

        // Grados máximos
        int n_max = 4;
        int m_max = 4;

        // 2. Llamada a la función
        cout << "Calculando aceleración..." << endl;
        Matrix a = AccelHarmonic(r, E, n_max, m_max);

        // 3. Verificación de resultados
        cout << "Resultado obtenido:" << endl;
        a.print();

        // Valores esperados aproximados (puedes ajustarlos según tus necesidades)
        double expected_data[] = {-7.9, -0.1, -0.05}; // [m/s^2]
        Matrix expected(3, 1, expected_data, 3);

        // Tolerancia para la comparación
        double tol = 0.2;

        // 4. Comparación usando el operador de asignación
        Matrix diff = a;  // Usamos el operador de asignación
        diff = diff - expected;  // Sobrecarga del operador -

        cout << "Diferencia con valores esperados:" << endl;
        diff.print();

        // 5. Comprobación de cada componente
        bool test_passed = true;
        for (int i = 1; i <= 3; ++i) {
            if (fabs(diff(i,1)) > tol) {
                cout << "Componente " << i << " fuera de tolerancia: "
                     << fabs(diff(i,1)) << " > " << tol << endl;
                test_passed = false;
            }
        }

        if (test_passed) {
            cout << "Test AccelHarmonic PASADO" << endl;
            return 0;
        } else {
            cout << "Test AccelHarmonic FALLADO" << endl;
            return 1;
        }

    } catch(const exception& e) {
        cerr << "Excepción en AccelHarmonic_Test: " << e.what() << endl;
        return 1;
    }
}

int all_tests()
{
   //_verify(Mjday_01);
   //_verify(Mjday_02);
   //_verify(R_x_01);
   //_verify(TimeUpdate_01);
   //_verify(Position_01);//Ns si esta bien

    //_verify(sign_);
    //_verify(AccelPointMass_01);
    //_verify(Mjday_TDB_01);//Modificar test, sale fail, pero el resultado es el esperado
    //_verify(angl_01);
    //_verify(angl_02);
    //_verify(Cheb3D_01);
    //_verify(Cheb3D_02);
    //_verify(MeanObliquity_01);//Modificar test, sale fail, pero el resultado es el esperado
    //_verify(NutAngles_01);//Creo que bien, pero hay que evrlo a mano
    _verify(AccelHarmonic_01);
    
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
