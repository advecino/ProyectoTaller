//
// Created by adria on 12/05/2025.
//

#ifndef PROYECTOTALLER_ANGLESDR_H
#define PROYECTOTALLER_ANGLESDR_H


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
#include "doubler.h"
#include "AuxParam.h"
#include <cstring>

struct AnglesDRResult {
    Matrix r2;  // Position vector at t2 [m]
    Matrix v2;  // Velocity vector at t2 [m/s]
    AnglesDRResult() : r2(3, 1), v2(3, 1){};
};

/**
 * @brief Orbit determination from three optical sightings, “double‐r” method.
 * @param az1,az2,az3   azimuths at t1/t2/t3 [rad]
 * @param el1,el2,el3   elevations at t1/t2/t3 [rad]
 * @param Mjd1,Mjd2,Mjd3  observation MJDs
 * @param rsite1,rsite2,rsite3  site ECEF positions [m]
 * @param params       auxiliary params (Mjd_UTC, Mjd_TT, n, m, flags)
 * @param eopdata      Earth orientation data matrix
 * @return r2,v2       position & velocity at t2
 */
AnglesDRResult anglesdr(
        double az1, double az2, double az3,
        double el1, double el2, double el3,
        double Mjd1, double Mjd2, double Mjd3,
        Matrix& rsite1,
        Matrix& rsite2,
        Matrix& rsite3,
        AuxParam& params,
        Matrix& eopdata
);

#endif //PROYECTOTALLER_ANGLESDR_H

