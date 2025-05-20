#include "../include/hgibbs.h"
#include <cmath>

/*
%--------------------------------------------------------------------------
%
%  hgibbs.m
%
%  this function implements the herrick-gibbs approximation for orbit
%  determination, and finds the middle velocity vector for the 3 given
%  position vectors.
%
%  inputs:
%    r1          - ijk position vector #1         m
%    r2          - ijk position vector #2         m
%    r3          - ijk position vector #3         m
%    Mjd1        - julian date of 1st sighting    days from 4713 bc
%    Mjd2        - julian date of 2nd sighting    days from 4713 bc
%    Mjd3        - julian date of 3rd sighting    days from 4713 bc
%
%  outputs:
%    v2          - ijk velocity vector for r2     m/s
%    theta       - angl between vectors           rad
%    error       - flag indicating success        'ok',...
%
%--------------------------------------------------------------------------*/


GibbsResult hgibbs(const Matrix& r1, const Matrix& r2, const Matrix& r3,
                    double Mjd1, double Mjd2, double Mjd3) {
    GibbsResult result;
    const double tolangle =  1.0 * M_PI/180.0;


    result.theta = 0.0;
    result.theta1 = 0.0;
    result.error = "          ok";
    result.v2 = Matrix(3, 1);


    const double magr1 = r1.norm();
    const double magr2 = r2.norm();
    const double magr3 = r3.norm();


    const double dt21 = (Mjd2 - Mjd1) * 86400.0;
    const double dt31 = (Mjd3 - Mjd1) * 86400.0;
    const double dt32 = (Mjd3 - Mjd2) * 86400.0;

    const Matrix p = Matrix::cross(r2, r3);
    const Matrix pn = unit(p);
    const Matrix r1n = unit( r1 );
    result.copa = asin(Matrix::dot(pn, r1n));

    if (std::abs(Matrix::dot(r1n, pn)) > 0.017452406) {
        result.error = "not coplanar";
        return result;
    }


    result.theta = angl( r1,r2 );
    result.theta1 = angl( r2,r3 );


    if ((result.theta > tolangle) || (result.theta1 > tolangle)) {
        result.error = "   angl > 1ø";
        return result;
    }


    const double term1 = -dt32 * (1.0/(dt21*dt31) + GM_Earth/(12.0*magr1*magr1*magr1));
    const double term2 = (dt32-dt21) * (1.0/(dt21*dt32) + GM_Earth/(12.0*magr2*magr2*magr2));
    const double term3 = dt21 * (1.0/(dt32*dt31) + GM_Earth/(12.0*magr3*magr3*magr3));


    result.v2 = term1*r1 + term2*r2 + term3*r3;

    return result;
}