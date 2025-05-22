//
// Created by adria on 12/05/2025.
//

#include "../include/anglesdr.h"
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


AnglesDRResult anglesdr(
        double az1, double az2, double az3,
        double el1, double el2, double el3,
        double Mjd1, double Mjd2, double Mjd3,
        Matrix &rsite1,
        Matrix &rsite2,
        Matrix &rsite3,
        AuxParam &params,
        Matrix &eopdata
) {

    double magr1in = 1.1 * R_Earth;
    double magr2in = 1.11 * R_Earth;
    const char direct = 'y';

    const double tol = 1e-8 * R_Earth;
    const double pctchg = 0.005;


    double t1 = (Mjd1 - Mjd2) * 86400.0;
    double t3 = (Mjd3 - Mjd2) * 86400.0;


    Matrix los1(3, 1), los2(3, 1), los3(3, 1);
    los1(1, 1) = std::cos(el1) * std::sin(az1);
    los1(2, 1) = std::cos(el1) * std::cos(az1);
    los1(3, 1) = std::sin(el1);
    los2(1, 1) = std::cos(el2) * std::sin(az2);
    los2(2, 1) = std::cos(el2) * std::cos(az2);
    los2(3, 1) = std::sin(el2);
    los3(1, 1) = std::cos(el3) * std::sin(az3);
    los3(2, 1) = std::cos(el3) * std::cos(az3);
    los3(3, 1) = std::sin(el3);

    auto geo1 = Geodetic(rsite1);
    auto geo2 = Geodetic(rsite2);
    auto geo3 = Geodetic(rsite3);

    Matrix M1 = LTC(geo1.longitude, geo1.latitude);
    Matrix M2 = LTC(geo2.longitude, geo2.latitude);
    Matrix M3 = LTC(geo3.longitude, geo3.latitude);

    // body-fixed system
    los1 = M1.transpuesta() * los1;
    los2 = M2.transpuesta() * los2;
    los3 = M3.transpuesta() * los3;



    // mean of date system (J2000)
    double Mjd_UTC = Mjd1;
    IERSResult ires = IERS(eopdata, Mjd_UTC, 'l');
    TimeDiffs td = timediff(ires.UT1_UTC, ires.TAI_UTC);
    // TT and UT1 in Modified Julian Days:
    double Mjd_TT = Mjd_UTC + td.TT_UTC / 86400.0;
    double Mjd_UT1 = Mjd_TT + (ires.UT1_UTC - td.TT_UTC) / 86400.0;

    // Precession‐nutation + polar
    Matrix P = PrecMatrix(MJD_J2000, Mjd_TT);
    Matrix N = NutMatrix(Mjd_TT);
    Matrix T = N * P;
    Matrix E = PoleMatrix(ires.x_pole, ires.y_pole)
               * GHAMatrix(Mjd_UT1)
               * T;

    los1 = E.transpuesta() * los1;
    rsite1 = E.transpuesta() * rsite1;


    Mjd_UTC = Mjd2;
    IERSResult ires2 = IERS(eopdata, Mjd_UTC, 'l');
    TimeDiffs td2 = timediff(ires2.UT1_UTC, ires2.TAI_UTC);
    // TT and UT1 in Modified Julian Days:
    Mjd_TT = Mjd_UTC + td2.TT_UTC / 86400.0;
    Mjd_UT1 = Mjd_TT + (ires2.UT1_UTC - td2.TT_UTC) / 86400.0;

    // Precession‐nutation + polar
    P = PrecMatrix(MJD_J2000, Mjd_TT);
    N = NutMatrix(Mjd_TT);
    T = N * P;
    E = PoleMatrix(ires2.x_pole, ires2.y_pole)
        * GHAMatrix(Mjd_UT1)
        * T;

    los2 = E.transpuesta() * los2;
    rsite2 = E.transpuesta() * rsite2;


    Mjd_UTC = Mjd3;
    IERSResult ires3 = IERS(eopdata, Mjd_UTC, 'l');
    TimeDiffs td3 = timediff(ires3.UT1_UTC, ires3.TAI_UTC);
    // TT and UT1 in Modified Julian Days:
    Mjd_TT = Mjd_UTC + td3.TT_UTC / 86400.0;
    Mjd_UT1 = Mjd_TT + (ires3.UT1_UTC - td3.TT_UTC) / 86400.0;

    // Precession‐nutation + polar
    P = PrecMatrix(MJD_J2000, Mjd_TT);
    N = NutMatrix(Mjd_TT);
    T = N * P;
    E = PoleMatrix(ires3.x_pole, ires3.y_pole)
        * GHAMatrix(Mjd_UT1)
        * T;

    los3 = E.transpuesta() * los3;
    rsite3 = E.transpuesta() * rsite3;



    // 7) prepare looping
    double magr1old = 99999999.9, magr2old = 99999999.9;
    double magrsite1 = rsite1.norm();
    double magrsite2 = rsite2.norm();
    double magrsite3 = rsite3.norm();
    double cc1 = 2.0 * Matrix::dot(los1, rsite1);
    double cc2 = 2.0 * Matrix::dot(los2, rsite2);
    int ktr = 0;

    Matrix r2(3, 1), r3(3, 1), v2(3, 1);
    double f1, f2, q1, magr1, magr2, a, deltae32;

    // 8) Newton‐Raphson loop
    while (fabs(magr1in - magr1old) > tol || fabs(magr2in - magr2old) > tol) {
        ktr++;
        auto R1 = doubler(cc1, cc2,
                          magrsite1, magrsite2,
                          magr1in, magr2in,
                          los1, los2, los3,
                          rsite1, rsite2, rsite3,
                          t1, t3, direct);
        r2       = R1.r2;
        r3       = R1.r3;
        f1       !=0.0 ? R1.f1:0;
        f2       !=0.0 ? R1.f2:0;
        q1       !=0.0 ? R1.q1:0;
        magr1    = R1.magr1;
        magr2    = R1.magr2;
        a        !=0.0 ? R1.a:0;
        deltae32 = R1.deltae32;

        double f = 1.0 - a / R1.magr2 * (1.0 - cos(deltae32));
        f !=0.0 ? f:0.0;
        double g = t3 - sqrt(pow(a, 3) / GM_Earth) * (deltae32 - sin(deltae32));
        v2 = (r3 - f * r2) / g;

        double magr1o = magr1in;
        magr1in = (1.0 + pctchg) * magr1in;
        double deltar1 = pctchg * magr1in;
        auto R2 = doubler(cc1, cc2, magrsite1, magrsite2, magr1in, magr2in,
                          los1, los2, los3, rsite1, rsite2, rsite3, t1, t3, direct);
        double pf1pr1 = (R2.f1 - f1) / deltar1;
        double pf2pr1 = (R2.f2 - f2) / deltar1;

        magr1in = magr1o;
        deltar1 = pctchg * magr1in;
        double magr2o = magr2in;
        magr2in = (1.0 + pctchg) * magr2in;
        double deltar2 = pctchg * magr2in;
        auto R3 = doubler(cc1, cc2, magrsite1, magrsite2, magr1in, magr2in,
                          los1, los2, los3, rsite1, rsite2, rsite3, t1, t3, direct);
        double pf1pr2 = (R3.f1 - f1) / deltar2;
        double pf2pr2 = (R3.f2 - f2) / deltar2;

        magr2in = magr2o;
        deltar2 = pctchg * magr2in;

        double delta = pf1pr1 * pf2pr2 - pf2pr1 * pf1pr2;
        double delta1 = pf2pr2 * f1 - pf1pr2 * f2;
        double delta2 = pf1pr1 * f2 - pf2pr1 * f1;

        deltar1 = -delta1 / delta;
        deltar2 = -delta2 / delta;

        magr1old = magr1in;
        magr2old = magr2in;

        magr1in = magr1in + deltar1;
        magr2in = magr2in + deltar2;

    }

    auto R4 = doubler(cc1, cc2, magrsite1, magrsite2, magr1in, magr2in,
                      los1, los2, los3, rsite1, rsite2, rsite3, t1, t3, direct);
    r2        = R4.r2;
    r3        = R4.r3;
    magr2     = R4.magr2;     // so f uses the correct radius
    a         = R4.a;
    deltae32  = R4.deltae32;

    double f = 1.0 - a/magr2*(1.0 - std::cos(deltae32));
    double g = t3 - std::sqrt(std::pow(a,3)/GM_Earth)*(deltae32 - std::sin(deltae32));
    v2 = (r3 - f*r2) / g;

    AnglesDRResult out;
    out.r2 = r2;
    out.v2 = v2;
    return out;
}