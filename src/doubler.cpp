//
// Created by adria on 11/05/2025.
//

#include "../include/doubler.h"
#include <cmath>
#include <stdexcept>



DoubleRResult doubler(double cc1, double cc2,
                      double magrsite1, double magrsite2,
                      double magr1in, double magr2in,
                      const Matrix& los1, const Matrix& los2, const Matrix& los3,
                      const Matrix& rsite1, const Matrix& rsite2, const Matrix& rsite3,
                      double t1, double t3, char direct) {
    DoubleRResult result;

    // 1) discriminantes
    double delta1 = cc1*cc1 - 4.0*(magrsite1*magrsite1 - magr1in*magr1in);
    double delta2 = cc2*cc2 - 4.0*(magrsite2*magrsite2 - magr2in*magr2in);
    if (delta1 < 0.0) throw std::runtime_error("Error: raíz cuadrada de número negativo en rho1");
    if (delta2 < 0.0) throw std::runtime_error("Error: raíz cuadrada de número negativo en rho2");

    // 2) raíces: en MATLAB usan la ‘+’ para coger la rama física
    double rho1 = (-cc1 + std::sqrt(delta1)) / 2.0;
    double rho2 = (-cc2 + std::sqrt(delta2)) / 2.0;

    // 3) vectores r1, r2
    Matrix r1 = los1 * rho1 + rsite1;
    result.r2 = los2 * rho2 + rsite2;

    // 4) magnitudes
    double magr1 = r1.norm();
    double magr2 = result.r2.norm();

    // 5) plano de órbita w
    Matrix w = Matrix::cross(r1, result.r2) * (1.0/(magr1*magr2));
    if (direct != 'y') w = w * (-1.0);

    // 6) rho3
    double los3w = Matrix::dot(los3, w);
    if (std::abs(los3w) < 1e-12) throw std::runtime_error("Error: división por cero en rho3");
    double rho3 = - Matrix::dot(rsite3, w) / los3w;
    result.r3 = los3 * rho3 + rsite3;
    double magr3 = result.r3.norm();

    // 7) ángulos geométricos
    double cos21 = Matrix::dot(result.r2, r1)/(magr2*magr1);
    double sin21 = Matrix::cross(result.r2, r1).norm()/(magr2*magr1);
    double dv21  = std::atan2(sin21, cos21);

    double cos31 = Matrix::dot(result.r3, r1)/(magr3*magr1);
    double sin31 = std::sqrt(std::max(0.0,1.0-cos31*cos31));
    double dv31  = std::atan2(sin31, cos31);

    double cos32 = Matrix::dot(result.r3, result.r2)/(magr3*magr2);
    double sin32 = Matrix::cross(result.r3, result.r2).norm()/(magr3*magr2);
    double dv32  = std::atan2(sin32, cos32);

    // 8) parámetros auxiliares p, c1, c3
    double p, c1, c3;
    if (dv31 > M_PI) {
        c1 = (magr2*sin32)/(magr1*sin31);
        c3 = (magr2*sin21)/(magr3*sin31);
        p  = (c1*magr1 + c3*magr3 - magr2)/(c1 + c3 - 1.0);
    } else {
        c1 = (magr1*sin31)/(magr2*sin32);
        c3 = (magr1*sin21)/(magr3*sin32);
        p  = (c3*magr3 - c1*magr2 + magr1)/(-c1 + c3 + 1.0);
    }

    double ecosv1 = p/magr1 - 1.0;
    double ecosv2 = p/magr2 - 1.0;
    double ecosv3 = p/magr3 - 1.0;

    double esinv2;
    if (std::abs(dv21 - M_PI) > 1e-6) {
        esinv2 = (-cos21 * ecosv2 + ecosv1)/sin21;
    } else {
        esinv2 = ( cos32 * ecosv2 - ecosv3)/sin31;
    }

    double e = std::sqrt(ecosv2*ecosv2 + esinv2*esinv2);
    result.a = p/(1.0 - e*e);

    // 9) anomalías excéntricas y tiempos Δe, f1, f2
    double n, s, c;
    if (e < 0.99) {
        n = std::sqrt(GM_Earth/std::pow(result.a,3));
        s = (magr2/p)*std::sqrt(1-e*e)*esinv2;
        c = (magr2/p)*(e*e + ecosv2);

        double sin32e = magr3/std::sqrt(result.a*p)*sin32 - magr3/p*(1.0-cos32)*s;
        double cos32e = 1.0 - magr2*magr3/(result.a*p)*(1.0-cos32);
        result.deltae32 = std::atan2(sin32e, cos32e);

        double sin21e = magr1/std::sqrt(result.a*p)*sin21 + magr1/p*(1.0-cos21)*s;
        double cos21e = 1.0 - magr2*magr1/(result.a*p)*(1.0-cos21);
        double deltae21 = std::atan2(sin21e, cos21e);

        double deltam32 = result.deltae32 + 2*s*std::pow(std::sin(result.deltae32/2),2) - c*std::sin(result.deltae32);
        double deltam12 = -deltae21      + 2*s*std::pow(std::sin(deltae21/2)     ,2) + c*std::sin(deltae21);

        result.f1 = t1 - deltam12/n;
        result.f2 = t3 - deltam32/n;
    } else {
        // caso hiperbólico… (idéntico al tuyo)
        n = std::sqrt(GM_Earth/std::pow(-result.a,3));
        s = (magr2/p)*std::sqrt(e*e - 1)*esinv2;
        c = (magr2/p)*(e*e + ecosv2);

        double sinh32h = magr3/std::sqrt(-result.a*p)*sin32 - magr3/p*(1.0-cos32)*s;
        double sinhh21 = magr1/std::sqrt(-result.a*p)*sin21 + magr1/p*(1.0-cos21)*s;

        double deltah32 = std::asinh(sinh32h);
        double deltah21 = std::asinh(sinhh21);

        double deltam32 = -deltah32      + 2*s*std::pow(std::sinh(deltah32/2),2) + c*std::sinh(deltah32);
        double deltam12 =  deltah21      + 2*s*std::pow(std::sinh(deltah21/2),2) - c*std::sinh(deltah21);

        result.deltae32 = deltah32;
        result.f1       = t1 - deltam12/n;
        result.f2       = t3 - deltam32/n;
    }

    result.q1 = std::sqrt(result.f1*result.f1 + result.f2*result.f2);
    return result;
}