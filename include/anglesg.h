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
#include "rpoly.h"
#include "AuxParam.h"

struct AnglesGResult {
    Matrix r2;  ///< Posición en t2 [m]
    Matrix v2;  ///< Velocidad en t2 [m/s]
    AnglesGResult(): r2(3,1), v2(3,1) {}
};

/**
 * @brief Determinación de órbita usando tres observaciones ópticas.
 */
AnglesGResult anglesg(
        double az1, double az2, double az3,
        double el1, double el2, double el3,
        double Mjd1, double Mjd2, double Mjd3,
        Matrix& Rs1,
        Matrix& Rs2,
        Matrix& Rs3,
        AuxParam& params,
        Matrix& eopdata
);


#endif //PROYECTOTALLER_ANGLESG_H
