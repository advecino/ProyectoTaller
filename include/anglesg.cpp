#include "anglesg.h"
#include "Geodetic.h"
#include "LTC.h"
#include "IERS.h"
#include "TimeDiff.h"
#include "PrecMatrix.h"
#include "NutMatrix.h"
#include "PoleMatrix.h"
#include "GHAMatrix.h"
#include "gibbs.h"
#include "hgibbs.h"
#include "elements.h"
#include "angl.h"
#include "global.h"
#include <vector>
#include <complex>
#include <stdexcept>
#include <cmath>

AnglesGResult::AnglesGResult()
        : r2(3,1), v2(3,1) {}

AnglesGResult::AnglesGResult(const Matrix& r, const Matrix& v)
        : r2(r), v2(v)
{
    if (r2.getFilas()!=3 || r2.getColumnas()!=1 ||
        v2.getFilas()!=3 || v2.getColumnas()!=1) {
        throw std::invalid_argument("AnglesGResult vectors must be 3×1");
    }
}

AnglesGResult anglesg(
        double az1, double az2, double az3,
        double el1, double el2, double el3,
        double Mjd1, double Mjd2, double Mjd3,
        const Matrix& Rs1_in,
        const Matrix& Rs2_in,
        const Matrix& Rs3_in
) {
    // Line‐of‐sight unit vectors
    Matrix L1(3,1), L2(3,1), L3(3,1);
    L1(1,1)=cos(el1)*sin(az1); L1(2,1)=cos(el1)*cos(az1); L1(3,1)=sin(el1);
    L2(1,1)=cos(el2)*sin(az2); L2(2,1)=cos(el2)*cos(az2); L2(3,1)=sin(el2);
    L3(1,1)=cos(el3)*sin(az3); L3(2,1)=cos(el3)*cos(az3); L3(3,1)=sin(el3);

    // Site geodetic coords & local tangent matrices
    auto g1 = Geodetic(Rs1_in);
    auto g2 = Geodetic(Rs2_in);
    auto g3 = Geodetic(Rs3_in);
    Matrix M1 = LTC(g1.longitude, g1.latitude);
    Matrix M2 = LTC(g2.longitude, g2.latitude);
    Matrix M3 = LTC(g3.longitude, g3.latitude);

    // Body‐fixed LOS
    Matrix Lb1 = M1.transpuesta() * L1;
    Matrix Lb2 = M2.transpuesta() * L2;
    Matrix Lb3 = M3.transpuesta() * L3;

    // Renombramos la función local para no chocar con std::transform
    auto doTransform = [&](double Mjd_UTC, const Matrix& Lb, const Matrix& Rs_in)
            -> std::pair<Matrix,Matrix>
    {
        IERSResult iers = IERS(eopdata, Mjd_UTC, 'l');
        TimeDiffs td = timediff(iers.UT1_UTC, iers.TAI_UTC);
        double Mjd_TT  = Mjd_UTC + td.TT_UTC/86400.0;
        double Mjd_UT1 = Mjd_TT + (iers.UT1_UTC - td.TT_UTC)/86400.0;
        Matrix E = PoleMatrix(iers.x_pole, iers.y_pole)
                   * GHAMatrix(Mjd_UT1)
                   * (NutMatrix(Mjd_TT) * PrecMatrix(MJD_J2000, Mjd_TT));
        return std::make_pair(E.transpuesta()*Lb,
                              E.transpuesta()*Rs_in);
    };

    // Ya no usamos structured bindings:
    std::pair<Matrix,Matrix> p1 = doTransform(Mjd1, Lb1, Rs1_in);
    Matrix Lm1 = p1.first;
    Matrix Rs1 = p1.second;

    std::pair<Matrix,Matrix> p2 = doTransform(Mjd2, Lb2, Rs2_in);
    Matrix Lm2 = p2.first;
    Matrix Rs2 = p2.second;

    std::pair<Matrix,Matrix> p3 = doTransform(Mjd3, Lb3, Rs3_in);
    Matrix Lm3 = p3.first;
    Matrix Rs3 = p3.second;

    // Ahora el resto del código puede usar Lm1, Lm2, Lm3, Rs1, Rs2, Rs3
    Matrix Lm(3,3), RsMat(3,3);
    for(int i=1;i<=3;++i){
        Lm(i,1)=Lm1(i,1); Lm(i,2)=Lm2(i,1); Lm(i,3)=Lm3(i,1);
        RsMat(i,1)=Rs1(i,1); RsMat(i,2)=Rs2(i,1); RsMat(i,3)=Rs3(i,1);
    }

    // Time intervals
    double tau1 = (Mjd1 - Mjd2)*86400.0;
    double tau3 = (Mjd3 - Mjd2)*86400.0;

    // Gauss coefficients
    double a1 =  tau3/(tau3 - tau1);
    double a3 = -tau1/(tau3 - tau1);
    double b1 =  tau3/(6*(tau3 - tau1))*((tau3 - tau1)*(tau3 - tau1) - tau3*tau3);
    double b3 = -tau1/(6*(tau3 - tau1))*((tau3 - tau1)*(tau3 - tau1) - tau1*tau1);

    // Build Lm and Rs matrices
    Matrix Rs(3,3);
    for(int i=1;i<=3;++i){
        Lm(i,1)=Lm1(i,1); Lm(i,2)=Lm2(i,1); Lm(i,3)=Lm3(i,1);
        Rs(i,1)=Rs1(i,1); Rs(i,2)=Rs2(i,1); Rs(i,3)=Rs3(i,1);
    }

    // Solve for distances via polynomial root‐finding (simplified)
    Matrix D = Lm.inversa() * Rs;
    double d1s = D(2,1)*a1 - D(2,2) + D(2,3)*a3;
    double d2s = D(2,1)*b1 + D(2,3)*b3;
    double Ccye = 2 * Matrix::dot(Lm2, Rs2);

    std::vector<double> poly(9);
    poly[0]= 1.0;
    poly[1]= 0.0;
    poly[2]= -(d1s*d1s + d1s*Ccye + Rs2.norm()*Rs2.norm());
    poly[3]= 0.0;
    poly[4]= 0.0;
    poly[5]= -GM_Earth*(d2s*Ccye + 2*d1s*d2s);
    poly[6]= 0.0;
    poly[7]= 0.0;
    poly[8]= -GM_Earth*GM_Earth*d2s*d2s;

    double bigr2 = 0.0;
    // crude scan for positive real root
    for(double r=6.37e6; r<1e8; r+=1e5){
        double v = poly[0]*std::pow(r,8)
                   + poly[2]*std::pow(r,6)
                   + poly[5]*std::pow(r,3)
                   + poly[8];
        if (std::fabs(v)<1e10) { bigr2=r; break; }
    }
    if(bigr2<=0) throw std::runtime_error("No valid root in anglesg");

    double u = GM_Earth/std::pow(bigr2,3);
    double C1= a1 + b1*u;
    double C2= -1.0;
    double C3= a3 + b3*u;
    Matrix tmp = -1*(D * Matrix({C1,C2,C3},3));
    double rho1 = tmp(1,1)/(a1+b1*u);
    double rho2 = -tmp(2,1);
    double rho3 =  tmp(3,1)/(a3+b3*u);

    // refine by one iteration
    Matrix r1 = Rs1 + rho1*Lm1;
    Matrix r2 = Rs2 + rho2*Lm2;
    Matrix r3 = Rs3 + rho3*Lm3;
    auto gib = gibbs(r1,r2,r3);
    Matrix v2 = (gib.error=="ok" ? gib.v2 : hgibbs(r1,r2,r3,Mjd1,Mjd2,Mjd3).v2);

    return AnglesGResult(r2, v2);
}
