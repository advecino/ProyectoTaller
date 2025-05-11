//
// Created by adria on 11/05/2025.
//

#include "../include/doubler.h"
#include <cmath>
#include <stdexcept>
#include <tuple>


DoubleRResult doubler(double cc1, double cc2,
                      double magrsite1, double magrsite2,
                      double magr1in, double magr2in,
                      const Matrix& los1, const Matrix& los2, const Matrix& los3,
                      const Matrix& rsite1, const Matrix& rsite2, const Matrix& rsite3,
                      double t1, double t3, char direct) {

    DoubleRResult result;

    double delta1 = cc1 * cc1 - 4 * (magrsite1*magrsite1 - magr1in*magr1in);
    double delta2 = cc2 * cc2 - 4 * (magrsite2*magrsite2 - magr2in*magr2in);
    if (delta1 < 0 || delta2 < 0)
        throw std::runtime_error("Error: raíz cuadrada de número negativo en rho1 o rho2.");

    double rho1 = (-cc1 + sqrt(delta1)) / 2.0;
    double rho2 = (-cc2 + sqrt(delta2)) / 2.0;

    Matrix r1 = los1 * rho1 + rsite1;
    result.r2 = los2 * rho2 + rsite2;

    result.magr1 = r1.norm();
    result.magr2 = result.r2.norm();

    Matrix w = Matrix(3,1);
    if (direct == 'y') {
        w = Matrix::cross(r1, result.r2) / (result.magr1 * result.magr2);
    } else {
        w = -1 * Matrix::cross(r1, result.r2) / (result.magr1 * result.magr2);
    }

    double los3w = Matrix::dot(los3, w);
    if (std::abs(los3w) < 1e-12)
        throw std::runtime_error("División por cero en cálculo de rho3.");

    double rho3 = -Matrix::dot(rsite3, w) / los3w;
    result.r3 = los3 * rho3 + rsite3;
    double magr3 = result.r3.norm();

    double cosdv21 = Matrix::dot(result.r2, r1) / (result.magr2 * result.magr1);
    double sindv21 = Matrix::cross(result.r2, r1).norm() / (result.magr2 * result.magr1);
    double dv21 = atan2(sindv21, cosdv21);

    double cosdv31 = Matrix::dot(result.r3, r1) / (magr3 * result.magr1);
    double sindv31 = sqrt(std::max(0.0, 1.0 - cosdv31 * cosdv31));
    double dv31 = atan2(sindv31, cosdv31);

    double cosdv32 = Matrix::dot(result.r3, result.r2) / (magr3 * result.magr2);
    double sindv32 = Matrix::cross(result.r3, result.r2).norm() / (magr3 * result.magr2);
    double dv32 = atan2(sindv32, cosdv32);

    double p, c1, c3;
    if (dv31 > M_PI) {
        c1 = (result.magr2 * sindv32) / (result.magr1 * sindv31);
        c3 = (result.magr2 * sindv21) / (magr3 * sindv31);
        p = (c1*result.magr1 + c3*magr3 - result.magr2) / (c1 + c3 - 1);
    } else {
        c1 = (result.magr1 * sindv31) / (result.magr2 * sindv32);
        c3 = (result.magr1 * sindv21) / (magr3 * sindv32);
        p = (c3*magr3 - c1*result.magr2 + result.magr1) / (-c1 + c3 + 1);
    }

    double ecosv1 = p/result.magr1 - 1;
    double ecosv2 = p/result.magr2 - 1;
    double ecosv3 = p/magr3 - 1;

    double esinv2;
    if (std::abs(dv21 - M_PI) > 1e-6) {
        esinv2 = (-cosdv21*ecosv2 + ecosv1) / sindv21;
    } else {
        esinv2 = (cosdv32*ecosv2 - ecosv3) / sindv31;
    }

    double e = sqrt(ecosv2*ecosv2 + esinv2*esinv2);
    result.a = p / (1 - e*e);

    double n, s, c;
    if (e < 0.99) {
        n = sqrt(GM_Earth / pow(result.a, 3));
        s = result.magr2/p * sqrt(1-e*e) * esinv2;
        c = result.magr2/p * (e*e + ecosv2);

        double sinde32 = magr3/sqrt(result.a*p)*sindv32 - magr3/p*(1 - cosdv32)*s;
        double cosde32 = 1 - result.magr2*magr3/(result.a*p)*(1 - cosdv32);
        result.deltae32 = atan2(sinde32, cosde32);

        double sinde21 = result.magr1/sqrt(result.a*p)*sindv21 + result.magr1/p*(1 - cosdv21)*s;
        double cosde21 = 1 - result.magr2*result.magr1/(result.a*p)*(1 - cosdv21);
        double deltae21 = atan2(sinde21, cosde21);

        double deltam32 = result.deltae32 + 2*s*pow(sin(result.deltae32/2), 2) - c*sin(result.deltae32);
        double deltam12 = -deltae21 + 2*s*pow(sin(deltae21/2), 2) + c*sin(deltae21);

        result.f1 = t1 - deltam12/n;
        result.f2 = t3 - deltam32/n;
    } else {
        n = sqrt(GM_Earth / pow(-result.a, 3));
        s = result.magr2/p * sqrt(e*e - 1) * esinv2;
        c = result.magr2/p * (e*e + ecosv2);

        double sindh32 = magr3/sqrt(-result.a*p)*sindv32 - magr3/p*(1 - cosdv32)*s;
        double sindh21 = result.magr1/sqrt(-result.a*p)*sindv21 + result.magr1/p*(1 - cosdv21)*s;

        double deltah32 = asinh(sindh32);
        double deltah21 = asinh(sindh21);

        double deltam32 = -deltah32 + 2*s*pow(sinh(deltah32/2), 2) + c*sinh(deltah32);
        double deltam12 = deltah21 + 2*s*pow(sinh(deltah21/2), 2) - c*sinh(deltah21);

        result.deltae32 = deltah32;
        result.f1 = t1 - deltam12/n;
        result.f2 = t3 - deltam32/n;
    }

    result.q1 = sqrt(result.f1*result.f1 + result.f2*result.f2);

    return result;
}