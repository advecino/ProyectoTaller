#include "../include/anglesg.h"
#include "../include/gibbs.h"
#include "../include/elements.h"
#include "../include/hgibbs.h"
#include <cmath>
#include <complex>
#include <stdexcept>

using std::sqrt;
using std::cos;
using std::sin;
using std::pow;

AnglesGResult anglesg(
        double az1, double az2, double az3,
        double el1, double el2, double el3,
        double Mjd1, double Mjd2, double Mjd3,
        Matrix &Rs1,
        Matrix &Rs2,
        Matrix &Rs3,
        AuxParam &params,
        Matrix &eopdata
) {
    // 1) topocentric line-of-sight
    Matrix L1(3, 1), L2(3, 1), L3(3, 1);
    L1(1, 1) = cos(el1) * sin(az1);
    L1(2, 1) = cos(el1) * cos(az1);
    L1(3, 1) = sin(el1);
    L2(1, 1) = cos(el2) * sin(az2);
    L2(2, 1) = cos(el2) * cos(az2);
    L2(3, 1) = sin(el2);
    L3(1, 1) = cos(el3) * sin(az3);
    L3(2, 1) = cos(el3) * cos(az3);
    L3(3, 1) = sin(el3);

    // 2) geodetic to body-fixed

    auto g2 = Geodetic(Rs2);
    auto g1 = Geodetic(Rs1);
    auto g3 = Geodetic(Rs3);
    Matrix M1 = LTC(g1.longitude, g1.latitude);
    Matrix M2 = LTC(g2.longitude, g2.latitude);
    Matrix M3 = LTC(g3.longitude, g3.latitude);
    Matrix Lb1 = M1.transpuesta() * L1;
    Matrix Lb2 = M1.transpuesta() * L2;
    Matrix Lb3 = M1.transpuesta() * L3;

    // 3) to inertial epochs
    auto toInertial = [&](double Mjd_obs, Matrix &Lb, Matrix &Rs) {
        IERSResult ire = IERS(eopdata, Mjd_obs, 'l');
        TimeDiffs td = timediff(ire.UT1_UTC, ire.TAI_UTC);
        double Mjd_TT = Mjd_obs + td.TT_UTC / 86400.0;
        double Mjd_UT1 = Mjd_TT + (ire.UT1_UTC - td.TT_UTC) / 86400.0;
        Matrix P = PrecMatrix(MJD_J2000, Mjd_TT);
        Matrix N = NutMatrix(Mjd_TT);
        Matrix T = N * P;
        Matrix E = PoleMatrix(ire.x_pole, ire.y_pole) * GHAMatrix(Mjd_UT1) * T;
        Lb = E.transpuesta() * Lb;
        Rs = E.transpuesta() * Rs;
    };
    toInertial(Mjd1, Lb1, Rs1);
    toInertial(Mjd2, Lb2, Rs2);
    toInertial(Mjd3, Lb3, Rs3);

    // 4) time intervals
    double tau1 = (Mjd1 - Mjd2) * 86400.0;
    double tau3 = (Mjd3 - Mjd2) * 86400.0;

    // 5) linear and cubic coefficients
    double a1 = tau3 / (tau3 - tau1);
    double a3 = -tau1 / (tau3 - tau1);
    double b1 = tau3 / (6 * (tau3 - tau1)) * ((tau3 - tau1) * (tau3 - tau1) - tau3 * tau3);
    double b3 = -tau1 / (6 * (tau3 - tau1)) * ((tau3 - tau1) * (tau3 - tau1) - tau1 * tau1);

    // 6) D matrix
    Matrix A(3, 3), B(3, 3);
    A.setColumn(1, Lb1);
    A.setColumn(2, Lb2);
    A.setColumn(3, Lb3);
    B.setColumn(1, Rs1);
    B.setColumn(2, Rs2);
    B.setColumn(3, Rs3);
    Matrix D = A.inversa() * B;

    // 7) compute d1s,d2s
    double Ccye = 2 * Matrix::dot(Lb2, Rs2);
    double d1s = D(2, 1) * a1 - D(2, 2) + D(2, 3) * a3;
    double d2s = D(2, 1) * b1 + D(2, 3) * b3;

    // 8) build polynomial coeffs
    double poly[9];
    poly[0] = 1.0;
    poly[1] = 0.0;
    poly[2] = -(d1s * d1s + d1s * Ccye + Rs2.norm() * Rs2.norm());
    poly[3] = 0.0;
    poly[4] = 0.0;
    poly[5] = -GM_Earth * (d2s * Ccye + 2 * d1s * d2s);
    poly[6] = 0.0;
    poly[7] = 0.0;
    poly[8] = -GM_Earth * GM_Earth * d2s * d2s;

    // 9) solve roots
    double zeror[8], zeroi[8];
    int m = real_poly_roots(poly, 8, zeror, zeroi);
    double bigr2 = -1.0;
    for (int i = 0; i < m; ++i) {
        if (fabs(zeroi[i]) < 1e-8 && zeror[i] > bigr2) bigr2 = zeror[i];
    }
    if (bigr2 <= 0) throw std::runtime_error("anglesg: no positive real root");

    //10) compute rhos
    double u = GM_Earth / (bigr2 * bigr2 * bigr2);
    double C_1 = a1 + b1 * u, C_2 = -1.0, C_3 = a3 + b3 * u;
    Matrix colC(3, 1);
    colC(1, 1) = C_1;
    colC(2, 1) = C_2;
    colC(3, 1) = C_3;
    Matrix tmp = D * colC * (-1.0);
    double rho1 = tmp(1, 1) / (a1 + b1 * u);
    double rho2 = -tmp(2, 1);
    double rho3 = tmp(3, 1) / (a3 + b3 * u);

    //11) iterative refinement
    double rhoold2;
    Matrix r1(3, 1), r2(3, 1), r3(3, 1), v2(3, 1);
    int ll = 0;
    do {
        if (ll++ > 8) break;
        rhoold2 = rho2;
        r1 = Rs1 + Lb1 * rho1;
        r2 = Rs2 + Lb2 * rho2;
        r3 = Rs3 + Lb3 * rho3;
        auto G = gibbs(r1, r2, r3);
        if (G.error != "ok" && G.copa < M_PI / 180.0) G = hgibbs(r1, r2, r3, Mjd1, Mjd2, Mjd3);
        v2 = G.v2;
        KeplerianElements ELEM;
        try {
            ELEM = elements(r2, v2);
        } catch (const std::invalid_argument &ex) {
            // órbita degenerada (circular o ecuatorial): salimos del refinamiento
            break;
        }
        double p = ELEM.p;
        double mag2 = r2.norm();
        double u_now = GM_Earth / (mag2 * mag2 * mag2);
        double rdot = Matrix::dot(r2, v2) / mag2;
        double udot = -3 * GM_Earth * rdot / (mag2 * mag2 * mag2 * mag2);
        double f1 = 1 - 0.5 * u_now * tau1 * tau1 - udot * tau1 * tau1 * tau1 / 6.0 -
                    u_now * u_now * tau1 * tau1 * tau1 * tau1 / 24.0 - udot * tau1 * tau1 * tau1 * tau1 * tau1 / 30.0;
        double g1 = tau1 - u_now * tau1 * tau1 * tau1 / 6.0 - udot * tau1 * tau1 * tau1 * tau1 / 12.0 -
                    u_now * u_now * tau1 * tau1 * tau1 * tau1 * tau1 / 120.0 -
                    udot * tau1 * tau1 * tau1 * tau1 * tau1 * tau1 / 120.0;
        double f3 = 1 - 0.5 * u_now * tau3 * tau3 - udot * tau3 * tau3 * tau3 / 6.0 -
                    u_now * u_now * tau3 * tau3 * tau3 * tau3 / 24.0 - udot * tau3 * tau3 * tau3 * tau3 * tau3 / 30.0;
        double g3 = tau3 - u_now * tau3 * tau3 * tau3 / 6.0 - udot * tau3 * tau3 * tau3 * tau3 / 12.0 -
                    u_now * u_now * tau3 * tau3 * tau3 * tau3 * tau3 / 120.0 -
                    udot * tau3 * tau3 * tau3 * tau3 * tau3 * tau3 / 120.0;
        C_1 = g3 / (f1 * g3 - f3 * g1);
        C_2 = -1.0;
        C_3 = -g1 / (f1 * g3 - f3 * g1);
        colC(1, 1) = C_1;
        colC(2, 1) = C_2;
        colC(3, 1) = C_3;
        tmp = D * colC * (-1.0);
        rho1 = tmp(1, 1) / (a1 + b1 * u);
        rho2 = -tmp(2, 1);
        rho3 = tmp(3, 1) / (a3 + b3 * u);
    } while (fabs(rho2 - rhoold2) > 1e-12);

    //12) final output
    AnglesGResult out;
    out.r2 = Rs2 + Lb2 * rho2;
    out.v2 = v2;
    return out;
}

