//
// Created by adria on 12/05/2025.
//

#ifndef PROYECTOTALLER_ANGLESG_H
#define PROYECTOTALLER_ANGLESG_H


#include "Matrix.h"
#include "Sat_const.h"
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
#include <vector>
#include <complex>

struct AnglesGResult {
    Matrix r2;  // Position vector at t2 [m]
    Matrix v2;  // Velocity vector at t2 [m/s]

    AnglesGResult() : r2(3,1), v2(3,1) {}
    AnglesGResult(const Matrix& r, const Matrix& v) : r2(r), v2(v) {
        if (r2.getFilas() != 3 || r2.getColumnas() != 1 ||
            v2.getFilas() != 3 || v2.getColumnas() != 1) {
            throw std::invalid_argument("Vectors must be 3x1");
        }
    }
};

/**
 * @brief Solves orbit determination using three optical sightings (Gauss method)
 *
 * @param az1,az2,az3 Azimuth angles [rad]
 * @param el1,el2,el3 Elevation angles [rad]
 * @param Mjd1,Mjd2,Mjd3 Modified Julian Dates
 * @param Rs1,Rs2,Rs3 Site position vectors [m]
 * @param eopdata Earth Orientation Parameters data
 * @return AnglesGResult Contains position and velocity vectors at t2
 */
AnglesGResult anglesg(double az1, double az2, double az3,
                      double el1, double el2, double el3,
                      double Mjd1, double Mjd2, double Mjd3,
                      const Matrix& Rs1, const Matrix& Rs2, const Matrix& Rs3,
                      Matrix& eopdata);


#endif //PROYECTOTALLER_ANGLESG_H
