//
// Created by adria on 10/05/2025.
//asta aqui

#include "../include/gibbs.h"
#include "../include/Sat_const.h"
#include "../include/angl.h"
#include <cmath>

/*
%--------------------------------------------------------------------------
%
%  gibbs.m
%
%  this function performs the gibbs method of orbit determination. this
%  method determines the velocity at the middle point of the 3 given
%  position vectors.
%
%  inputs:
%    r1          - ijk position vector #1         m
%    r2          - ijk position vector #2         m
%    r3          - ijk position vector #3         m
%
%  outputs:
%    v2          - ijk velocity vector for r2     m/s
%    theta       - angl between vectors           rad
%    error       - flag indicating success        'ok',...
%
%--------------------------------------------------------------------------*/



GibbsResult gibbs(const Matrix& r1, const Matrix& r2, const Matrix& r3) {
    GibbsResult result;
    const double small = 1e-8;

    // Inicialización de resultados
    result.theta = 0.0;
    result.theta1 = 0.0;
    result.error = "          ok";
    result.v2 = Matrix(3, 1);

    // Calcular magnitudes
    const double magr1 = r1.norm();
    const double magr2 = r2.norm();
    const double magr3 = r3.norm();

    // Productos cruz
    Matrix p = Matrix::cross(r2, r3);
    Matrix q = Matrix::cross(r3, r1);
    Matrix w = Matrix::cross(r1, r2);

    // Vectores unitarios
    const Matrix pn = unit(p);
    const Matrix r1n = unit(r1);

    // Ángulo de coplanaridad
    result.copa = asin(Matrix::dot(pn, r1n));

    // Verificar coplanaridad
    if (std::abs(Matrix::dot(r1n, pn)) > 0.017452406) {
        result.error = "not coplanar";
        return result;
    }

    // Vectores intermedios
    const Matrix d = p + q + w;
    const double magd = d.norm();

    const Matrix n = magr1*p + magr2*q + magr3*w;
    const double magn = n.norm();

    const Matrix nn = unit(n);
    const Matrix dn = unit(d);
/*
    % -------------------------------------------------------------
    % determine if  the orbit is possible. both d and n must be in
    % the same direction, and non-zero.
    % -------------------------------------------------------------*/
    if ((std::abs(magd) < small) || (std::abs(magn) < small) ||
        (Matrix::dot(nn, dn) < small)) {
        result.error = "  impossible";
        return result;
    }else{
        result.theta = angl(r1,r2);
        result.theta1 = angl(r2,r3);

        //% ----------- perform gibbs method to find v2 -----------
        const double r1mr2 = magr1 - magr2;
        const double r3mr1 = magr3 - magr1;
        const double r2mr3 = magr2 - magr3;

        const Matrix s = (r1mr2*r3 + r3mr1*r2 + r2mr3*r1);
        const Matrix b = Matrix::cross(d, r2);
        const double l = sqrt(GM_Earth / (magd * magn));
        const double tover2 = l / magr2;
        result.v2 = tover2 * b + l * s;
    }

    return result;
}


