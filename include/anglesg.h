//
// Created by adria on 12/05/2025.
//

#ifndef PROYECTOTALLER_ANGLESG_H
#define PROYECTOTALLER_ANGLESG_H


#include "Matrix.h"
#include "Sat_const.h"

/**
 * @brief Result of Gauss angles-only orbit determination
 */
struct AnglesGResult {
    Matrix r2;  ///< Position at t2 [m]
    Matrix v2;  ///< Velocity at t2 [m/s]
    AnglesGResult();
    AnglesGResult(const Matrix& r, const Matrix& v);
};

/**
 * @brief Solves orbit determination from three optical sightings (Gauss method)
 * @param az1,az2,az3  Azimuths at t1,t2,t3 [rad]
 * @param el1,el2,el3  Elevations at t1,t2,t3 [rad]
 * @param Mjd1,Mjd2,Mjd3  Modified Julian Dates of observations
 * @param Rs1,Rs2,Rs3    Site position vectors in ICRF [m]
 * @param eopdata        Earth‐orientation parameters matrix
 * @return Position and velocity at t2
 */
AnglesGResult anglesg(
        double az1, double az2, double az3,
        double el1, double el2, double el3,
        double Mjd1, double Mjd2, double Mjd3,
        const Matrix& Rs1,
        const Matrix& Rs2,
        const Matrix& Rs3
);


#endif //PROYECTOTALLER_ANGLESG_H
