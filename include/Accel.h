//
// Created by adria on 12/05/2025.
//

#ifndef PROYECTOTALLER_ACCEL_H
#define PROYECTOTALLER_ACCEL_H

#include "Matrix.h"

#include "IERS.h"
#include "TimeDiff.h"
#include "Mjday_TDB.h"
#include "JPL_Eph_DE430.h"
#include "AccelHarmonic.h"
#include "AccelPointMass.h"
#include "PrecMatrix.h"
#include "NutMatrix.h"
#include "PoleMatrix.h"
#include "GHAMatrix.h"
#include "Sat_const.h"


Matrix Accel(double x, const Matrix& Y, const Matrix& eopdata);

#endif //PROYECTOTALLER_ACCEL_H
