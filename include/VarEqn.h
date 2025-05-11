//
// Created by adria on 12/05/2025.
//

#ifndef PROYECTOTALLER_VAREQN_H
#define PROYECTOTALLER_VAREQN_H


#include "Matrix.h"
#include "IERS.h"
#include "TimeDiff.h"
#include "PrecMatrix.h"
#include "NutMatrix.h"
#include "PoleMatrix.h"
#include "GHAMatrix.h"
#include "AccelHarmonic.h"
#include "G_AccelHarmonic.h"
#include "Sat_const.h"
#include <vector>

/**
 * @brief Computes the variational equations (derivative of state vector and state transition matrix)
 *
 * @param x Time since epoch in [s]
 * @param yPhi (6+36)-dim vector comprising state vector (y) and state transition matrix (Phi)
 * @param eopdata Earth Orientation Parameters data
 * @return Matrix Derivative of yPhi (42x1)
 */
Matrix VarEqn(double x, const Matrix& yPhi, const Matrix& eopdata);

#endif //PROYECTOTALLER_VAREQN_H
