//
// Created by adria on 12/05/2025.
//

#include "anglesdr.h"
#include <cmath>
#include <stdexcept>

/*
%---------------------------------------------------------------------------
%
%  anglesdr.m
%
%  this function solves the problem of orbit determination using three
%  optical sightings.
%
%  inputs:
%    az1      - azimuth at t1               rad
%    az2      - azimuth at t2               rad
%    az3      - azimuth at t3               rad
%    el1      - elevation at t1             rad
%    el2      - elevation at t2             rad
%    el3      - elevation at t3             rad
%    Mjd1     - Modified julian date of t1
%    Mjd2     - Modified julian date of t2
%    Mjd3     - Modified julian date of t3
%    rsite1   - ijk site1 position vector   m
%    rsite2   - ijk site2 position vector   m
%    rsite3   - ijk site3 position vector   m
%
%  outputs:
%    r        - ijk position vector at t2   m
%    v        - ijk velocity vector at t2   m/s
%
% Last modified:   2015/08/12   M. Mahooti
%
%---------------------------------------------------------------------------
*/


AnglesDRResult anglesdr(double az1, double az2, double az3,
                        double el1, double el2, double el3,
                        double Mjd1, double Mjd2, double Mjd3,
                        const Matrix& rsite1, const Matrix& rsite2, const Matrix& rsite3,
                        Matrix& eopdata) {

    // Validate inputs
    if (rsite1.getFilas() != 3 || rsite1.getColumnas() != 1 ||
        rsite2.getFilas() != 3 || rsite2.getColumnas() != 1 ||
        rsite3.getFilas() != 3 || rsite3.getColumnas() != 1) {
        throw std::invalid_argument("Site position vectors must be 3x1");
    }

    AnglesDRResult result;
    const double tol = 1e-8 * R_Earth;  // From Sat_const
    const double pctchg = 0.005;
    char direct = 'y';

    // Initial guesses for magnitudes
    double magr1in = 1.1 * R_Earth;
    double magr2in = 1.11 * R_Earth;

    // Time differences
    double t1 = (Mjd1 - Mjd2) * 86400.0;
    double t3 = (Mjd3 - Mjd2) * 86400.0;

    // Line-of-sight vectors
    Matrix los1(3, 1), los2(3, 1), los3(3, 1);
    los1(1,1) = cos(el1)*sin(az1); los1(2,1) = cos(el1)*cos(az1); los1(3,1) = sin(el1);
    los2(1,1) = cos(el2)*sin(az2); los2(2,1) = cos(el2)*cos(az2); los2(3,1) = sin(el2);
    los3(1,1) = cos(el3)*sin(az3); los3(2,1) = cos(el3)*cos(az3); los3(3,1) = sin(el3);


    GeodeticCoords geodeticCoords1;
    geodeticCoords1 = Geodetic(rsite1);
    GeodeticCoords geodeticCoords2;
    geodeticCoords2 = Geodetic(rsite2);
    GeodeticCoords geodeticCoords3;
    geodeticCoords3 = Geodetic(rsite3);

    // Get local tangent coordinate matrices
    Matrix M1 = LTC(geodeticCoords1.longitude, geodeticCoords1.latitude);
    Matrix M2 = LTC(geodeticCoords2.longitude, geodeticCoords2.latitude);
    Matrix M3 = LTC(geodeticCoords3.longitude, geodeticCoords3.latitude);

    // Transform to body-fixed system
    los1 = M1.transpuesta() * los1;
    los2 = M2.transpuesta() * los2;
    los3 = M3.transpuesta() * los3;

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
    los1 = E.transpuesta() * los1;
    Matrix rsite1_J2000 = E.transpuesta() * rsite1;

    // Time 2 transformation
    Mjd_UTC = Mjd2;
    IERSResult iers2 = IERS(eopdata, Mjd_UTC, 'l');
    TimeDiffs td2 = timediff(iers2.UT1_UTC, iers2.TAI_UTC);
    Mjd_TT = Mjd_UTC + td2.TT_UTC/86400.0;
    Mjd_UT1 = Mjd_TT + (iers2.UT1_UTC - td2.TT_UTC)/86400.0;

    E = PoleMatrix(iers2.x_pole, iers2.y_pole) * GHAMatrix(Mjd_UT1) *
        (NutMatrix(Mjd_TT) * PrecMatrix(MJD_J2000, Mjd_TT));

    los2 = E.transpuesta() * los2;
    Matrix rsite2_J2000 = E.transpuesta() * rsite2;

    // Time 3 transformation
    Mjd_UTC = Mjd3;
    IERSResult iers3 = IERS(eopdata, Mjd_UTC, 'l');
    TimeDiffs td3 = timediff(iers3.UT1_UTC, iers3.TAI_UTC);
    Mjd_TT = Mjd_UTC + td3.TT_UTC/86400.0;
    Mjd_UT1 = Mjd_TT + (iers3.UT1_UTC - td3.TT_UTC)/86400.0;
    E = PoleMatrix(iers3.x_pole, iers3.y_pole) * GHAMatrix(Mjd_UT1) *
        (NutMatrix(Mjd_TT) * PrecMatrix(MJD_J2000, Mjd_TT));
    los3 = E.transpuesta() * los3;
    Matrix rsite3_J2000 = E.transpuesta() * rsite3;

    // Prepare for iteration
    double magr1old = 99999999.9;
    double magr2old = 99999999.9;
    double magrsite1 = rsite1_J2000.norm();
    double magrsite2 = rsite2_J2000.norm();
    double magrsite3 = rsite3_J2000.norm();
    double cc1 = 2.0 * Matrix::dot(los1, rsite1_J2000);
    double cc2 = 2.0 * Matrix::dot(los2, rsite2_J2000);
    int ktr = 0;

    // Main iteration loop
    while (fabs(magr1in - magr1old) > tol || fabs(magr2in - magr2old) > tol) {
        ktr++;

        DoubleRResult dr = doubler(cc1, cc2, magrsite1, magrsite2, magr1in, magr2in,
                                  los1, los2, los3, rsite1_J2000, rsite2_J2000, rsite3_J2000,
                                  t1, t3, direct);

        double f = 1.0 - dr.a/dr.magr2*(1.0 - cos(dr.deltae32));
        double g = t3 - sqrt(pow(dr.a,3)/GM_Earth)*(dr.deltae32 - sin(dr.deltae32));
        result.v2 = (dr.r3 - f*dr.r2)/g;

        // Partial derivatives calculation
        double magr1o = magr1in;
        magr1in = (1.0 + pctchg) * magr1in;
        double deltar1 = pctchg * magr1in;
        DoubleRResult dr_dr1 = doubler(cc1, cc2, magrsite1, magrsite2, magr1in, magr2in,
                                      los1, los2, los3, rsite1_J2000, rsite2_J2000, rsite3_J2000,
                                      t1, t3, direct);
        double pf1pr1 = (dr_dr1.f1 - dr.f1)/deltar1;
        double pf2pr1 = (dr_dr1.f2 - dr.f2)/deltar1;

        magr1in = magr1o;
        double magr2o = magr2in;
        magr2in = (1.0 + pctchg) * magr2in;
        double deltar2 = pctchg * magr2in;
        DoubleRResult dr_dr2 = doubler(cc1, cc2, magrsite1, magrsite2, magr1in, magr2in,
                                      los1, los2, los3, rsite1_J2000, rsite2_J2000, rsite3_J2000,
                                      t1, t3, direct);
        double pf1pr2 = (dr_dr2.f1 - dr.f1)/deltar2;
        double pf2pr2 = (dr_dr2.f2 - dr.f2)/deltar2;

        magr2in = magr2o;
        double delta = pf1pr1*pf2pr2 - pf2pr1*pf1pr2;
        double delta1 = pf2pr2*dr.f1 - pf1pr2*dr.f2;
        double delta2 = pf1pr1*dr.f2 - pf2pr1*dr.f1;

        deltar1 = -delta1/delta;
        deltar2 = -delta2/delta;

        magr1old = magr1in;
        magr2old = magr2in;
        magr1in += deltar1;
        magr2in += deltar2;
    }

    // Final calculation
    DoubleRResult final_dr = doubler(cc1, cc2, magrsite1, magrsite2, magr1in, magr2in,
                                    los1, los2, los3, rsite1_J2000, rsite2_J2000, rsite3_J2000,
                                    t1, t3, direct);

    double f = 1.0 - final_dr.a/final_dr.magr2*(1.0 - cos(final_dr.deltae32));
    double g = t3 - sqrt(pow(final_dr.a,3)/GM_Earth)*(final_dr.deltae32 - sin(final_dr.deltae32));
    result.r2 = final_dr.r2;
    result.v2 = (final_dr.r3 - f*final_dr.r2)/g;

    return result;
}