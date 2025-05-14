//
// Created by adria on 12/05/2025.
//

#include "anglesg.h"


#include "anglesg.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>

AnglesGResult anglesg(double az1, double az2, double az3,
                      double el1, double el2, double el3,
                      double Mjd1, double Mjd2, double Mjd3,
                      const Matrix& Rs1, const Matrix& Rs2, const Matrix& Rs3,
                      Matrix& eopdata) {

    // Validate inputs
    if (Rs1.getFilas() != 3 || Rs1.getColumnas() != 1 ||
        Rs2.getFilas() != 3 || Rs2.getColumnas() != 1 ||
        Rs3.getFilas() != 3 || Rs3.getColumnas() != 1) {
        throw std::invalid_argument("Site position vectors must be 3x1");
    }

    AnglesGResult result;

    // Line-of-sight vectors
    Matrix L1(3,1), L2(3,1), L3(3,1);
    L1(1,1) = cos(el1)*sin(az1); L1(2,1) = cos(el1)*cos(az1); L1(3,1) = sin(el1);
    L2(1,1) = cos(el2)*sin(az2); L2(2,1) = cos(el2)*cos(az2); L2(3,1) = sin(el2);
    L3(1,1) = cos(el3)*sin(az3); L3(2,1) = cos(el3)*cos(az3); L3(3,1) = sin(el3);

    // Convert site positions to geodetic coordinates
    GeodeticCoords geo1;
    GeodeticCoords geo2;
    GeodeticCoords geo3;
    geo1 = Geodetic(Rs1);
    geo2 = Geodetic(Rs2);
    geo3 = Geodetic(Rs3);

    // Get local tangent coordinate matrices
    Matrix M1 = LTC(geo1.longitude, geo1.latitude);
    Matrix M2 = LTC(geo2.longitude, geo2.latitude);
    Matrix M3 = LTC(geo3.longitude, geo3.latitude);

    // Transform to body-fixed system
    Matrix Lb1 = M1.transpuesta() * L1;
    Matrix Lb2 = M2.transpuesta() * L2;
    Matrix Lb3 = M3.transpuesta() * L3;

    // Transform to J2000 system for each time

    double Mjd_UTC, Mjd_TT, Mjd_UT1;

    // Time 1 transformation
    Mjd_UTC = Mjd1;
    IERSResult iers1 = IERS(eopdata, Mjd_UTC, 'l');
    TimeDiffs td1 = timediff(iers1.UT1_UTC, iers1.TAI_UTC);
    Mjd_TT = Mjd_UTC + td1.TT_UTC/86400.0;
    Mjd_UT1 = Mjd_TT + (iers1.UT1_UTC - td1.TT_UTC)/86400.0;
    Matrix E = PoleMatrix(iers1.x_pole, iers1.y_pole) * GHAMatrix(Mjd_UT1) *
            (NutMatrix(Mjd_TT) * PrecMatrix(MJD_J2000, Mjd_TT));
    Matrix Lm1 = E.transpuesta() * Lb1;
    Matrix Rs1_J2000 = E.transpuesta() * Rs1;

    // Time 2 transformation
    Mjd_UTC = Mjd2;
    IERSResult iers2 = IERS(eopdata, Mjd_UTC, 'l');
    TimeDiffs td2 = timediff(iers2.UT1_UTC, iers2.TAI_UTC);
    Mjd_TT = Mjd_UTC + td2.TT_UTC/86400.0;
    Mjd_UT1 = Mjd_TT + (iers2.UT1_UTC - td2.TT_UTC)/86400.0;
    E = PoleMatrix(iers2.x_pole, iers2.y_pole) * GHAMatrix(Mjd_UT1) *
        (NutMatrix(Mjd_TT) * PrecMatrix(MJD_J2000, Mjd_TT));
    Matrix Lm2 = E.transpuesta() * Lb2;
    Matrix Rs2_J2000 = E.transpuesta() * Rs2;

    // Time 3 transformation
    Mjd_UTC = Mjd3;
    IERSResult iers3 = IERS(eopdata, Mjd_UTC, 'l');
    TimeDiffs td3 = timediff(iers3.UT1_UTC, iers3.TAI_UTC);
    Mjd_TT = Mjd_UTC + td3.TT_UTC/86400.0;
    Mjd_UT1 = Mjd_TT + (iers3.UT1_UTC - td3.TT_UTC)/86400.0;
    E = PoleMatrix(iers3.x_pole, iers3.y_pole) * GHAMatrix(Mjd_UT1) *
        (NutMatrix(Mjd_TT) * PrecMatrix(MJD_J2000, Mjd_TT));
    Matrix Lm3 = E.transpuesta() * Lb3;
    Matrix Rs3_J2000 = E.transpuesta() * Rs3;

    // Time intervals
    double tau1 = (Mjd1 - Mjd2) * 86400.0;
    double tau3 = (Mjd3 - Mjd2) * 86400.0;

    // Gauss method coefficients
    double a1 = tau3 / (tau3 - tau1);
    double a3 = -tau1 / (tau3 - tau1);
    double b1 = tau3 / (6 * (tau3 - tau1)) * (pow(tau3 - tau1, 2) - pow(tau3, 2));
    double b3 = -tau1 / (6 * (tau3 - tau1)) * (pow(tau3 - tau1, 2) - pow(tau1, 2));

    // Matrix operations for Gauss method
    Matrix Lm(3,3);
    Lm.setColumn(1, Lm1);
    Lm.setColumn(2, Lm2);
    Lm.setColumn(3, Lm3);

    Matrix Rs(3,3);
    Rs.setColumn(1, Rs1_J2000);
    Rs.setColumn(2, Rs2_J2000);
    Rs.setColumn(3, Rs3_J2000);

    Matrix D = Lm.inversa() * Rs;

    double d1s = D(2,1)*a1 - D(2,2) + D(2,3)*a3;
    double d2s = D(2,1)*b1 + D(2,3)*b3;
    double Ccye = 2 * Matrix::dot(Lm2, Rs2_J2000);

    // Polynomial coefficients for r^8 equation
    std::vector<double> poly(10);
    poly[0] = 1.0;  // R2^8 term
    poly[1] = 0.0;
    poly[2] = -(pow(d1s, 2) + d1s*Ccye + pow(Rs2_J2000.norm(), 2));
    poly[3] = 0.0;
    poly[4] = 0.0;
    poly[5] = -GM_Earth * (d2s*Ccye + 2*d1s*d2s);
    poly[6] = 0.0;
    poly[7] = 0.0;
    poly[8] = -pow(GM_Earth, 2) * pow(d2s, 2);

    // Find real positive root (simplified root finding)
    double bigr2 = 0.0;
    // In practice, you'd use a proper root-finding algorithm here
    // This is a simplified approach for demonstration
    for (double r = 7000000.0; r < 100000000.0; r += 100000.0) {
        double val = poly[0]*pow(r,8) + poly[2]*pow(r,6) + poly[5]*pow(r,3) + poly[8];
        if (fabs(val) < 1e10) {  // Crude check for root
            bigr2 = r;
            break;
        }
    }

    if (bigr2 <= 0.0) {
        throw std::runtime_error("Failed to find valid orbit radius solution");
    }

    double u = GM_Earth / pow(bigr2, 3);
    double C1 = a1 + b1*u;
    double C2 = -1.0;
    double C3 = a3 + b3*u;

    Matrix temp = -1*D * Matrix({C1, C2, C3}, 3);
    double rho1 = temp(1,1) / (a1 + b1*u);
    double rho2 = -temp(2,1);
    double rho3 = temp(3,1) / (a3 + b3*u);

    // Iterative improvement
    double rhoold1 = rho1;
    double rhoold2 = rho2;
    double rhoold3 = rho3;
    int ll = 0;
    const int max_iter = 20;
    double tolerance = 1e-12;

    while ((fabs(rhoold2 - rho2) > tolerance) && ll < max_iter) {
        ll++;
        rhoold2 = rho2;

        Matrix r1 = Rs1_J2000 + rhoold1*Lm1;
        Matrix r2 = Rs2_J2000 + rhoold2*Lm2;
        Matrix r3 = Rs3_J2000 + rhoold3*Lm3;

        double magr1 = r1.norm();
        double magr2 = r2.norm();
        double magr3 = r3.norm();

        // Gibbs/HGibbs orbit determination
        GibbsResult gibbs_res = gibbs(r1, r2, r3);
        Matrix velocity(3,3);

        if (gibbs_res.error != "ok") {
            HGibbsResult hgibbs_res = hgibbs(r1, r2, r3, Mjd1, Mjd2, Mjd3);
            velocity = hgibbs_res.v2;
        } else {
            velocity = gibbs_res.v2;
        }

// Create combined state vector [r2; v2] for elements()
        Matrix y(6, 1);
// Position part (r2)
        for (int i = 1; i <= 3; i++) {
            y(i, 1) = r2(i, 1);
        }
// Velocity part
        for (int i = 1; i <= 3; i++) {
            y(i+3, 1) = velocity(i, 1);
        }

// Get orbital elements
        KeplerianElements elem = elements(y);

        double f3,g3,f1,g1;

        if (ll <= 8) {
            u = GM_Earth / pow(magr2, 3);
            double rdot = Matrix::dot(r2, gibbs_res.v2) / magr2;
            double udot = (-3 * GM_Earth * rdot) / pow(magr2, 4);

            double tausqr = pow(tau1, 2);
            f1 = 1 - 0.5*u*tausqr - (1.0/6)*udot*tausqr*tau1
                        - (1.0/24)*pow(u,2)*pow(tausqr,2)
                        - (1.0/30)*u*udot*pow(tausqr,2)*tau1;
            g1 = tau1 - (1.0/6)*u*tau1*tausqr - (1.0/12)*udot*pow(tausqr,2)
                        - (1.0/120)*pow(u,2)*pow(tausqr,2)*tau1
                        - (1.0/120)*u*udot*pow(tausqr,3);

            tausqr = pow(tau3, 2);
            f3 = 1 - 0.5*u*tausqr - (1.0/6)*udot*tausqr*tau3
                        - (1.0/24)*pow(u,2)*pow(tausqr,2)
                        - (1.0/30)*u*udot*pow(tausqr,2)*tau3;
            g3 = tau3 - (1.0/6)*u*tau3*tausqr - (1.0/12)*udot*pow(tausqr,2)
                        - (1.0/120)*pow(u,2)*pow(tausqr,2)*tau3
                        - (1.0/120)*u*udot*pow(tausqr,3);
        } else {
            double theta = angl(r1, r2);
            double theta1 = angl(r2, r3);

            double f1 = 1 - (magr1*(1 - cos(theta)) / elem.p);
            double g1 = (magr1*magr2*sin(-theta)) / sqrt(elem.p);
            double f3 = 1 - (magr3*(1 - cos(theta1)) / elem.p);
            double g3 = (magr3*magr2*sin(theta1)) / sqrt(elem.p);
        }

        C1 = g3 / (f1*g3 - f3*g1);
        C2 = -1.0;
        C3 = -g1 / (f1*g3 - f3*g1);

        double H1 = GM_Earth * tau3 / 12.0;
        double H3 = -GM_Earth * tau1 / 12.0;
        double H2 = H1 - H3;

        double G1 = -tau3 / (tau1 * (tau3 - tau1));
        double G3 = -tau1 / (tau3 * (tau3 - tau1));
        double G2 = G1 - G3;

        double D1 = G1 + H1/pow(magr1, 3);
        double D2 = G2 + H2/pow(magr2, 3);
        double D3 = G3 + H3/pow(magr3, 3);

        temp = -1*Matrix({D1, D2, D3}, 3) * Matrix({C1, C2, C3}, 3);
        rhoold1 = temp(1,1) / (a1 + b1*u);
        rhoold2 = -temp(1,1);
        rhoold3 = temp(1,1) / (a3 + b3*u);

        r1 = Rs1_J2000 + rhoold1*Lm1;
        r2 = Rs2_J2000 + rhoold2*Lm2;
        r3 = Rs3_J2000 + rhoold3*Lm3;

        if(!((fabs(rhoold2 - rho2) > tolerance) && ll < max_iter)){
            GibbsResult gibbs2 = gibbs(r1, r2, r3);
            result.r2 = Rs2_J2000 + rho2*Lm2;
            //result.v2 = gibbs(r1, r2, r3).v2;
        }
    }



    return result;
}