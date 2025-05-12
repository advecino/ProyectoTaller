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
#include <cstring>

struct AnglesDRResult {
    Matrix r2;  // Position vector at t2 [m]
    Matrix v2;  // Velocity vector at t2 [m/s]
    AnglesDRResult() : r2(3, 1), v2(3, 1){};
};

/**
 * @brief Solves orbit determination using three optical sightings
 *
 * @param az1,az2,az3 Azimuth angles [rad]
 * @param el1,el2,el3 Elevation angles [rad]
 * @param Mjd1,Mjd2,Mjd3 Modified Julian Dates
 * @param rsite1,rsite2,rsite3 Site position vectors [m]
 * @param eopdata Earth Orientation Parameters data
 * @return AnglesDRResult Contains position and velocity vectors at t2
 */
AnglesDRResult anglesdr(double az1, double az2, double az3,
                        double el1, double el2, double el3,
                        double Mjd1, double Mjd2, double Mjd3,
                        const Matrix& rsite1, const Matrix& rsite2, const Matrix& rsite3,
                        const Matrix& eopdata);

#endif //PROYECTOTALLER_ANGLESDR_H

