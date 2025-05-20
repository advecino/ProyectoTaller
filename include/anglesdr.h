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
    Matrix r2;  ///<  Vector de posición en t2 [m]
    Matrix v2;  ///<  Vector de velocidad en t2 [m/s]
    AnglesDRResult() : r2(3, 1), v2(3, 1){};
};

/**
 * @brief Determinación de órbita usando tres observaciones ópticas.
 * @param az1, az2, az3 Acimutes en t1, t2 y t3 [rad]
 * @param el1, el2, el3 Elevaciones en t1, t2 y t3 [rad]
 * @param Mjd1, Mjd2, Mjd3 Fechas julianas modificadas de observación
 * @param rsite1, rsite2, rsite3 Vectores de posición del sitio en ECEF [m]
 * @param params Parámetros auxiliares (fechas, orden de armónicos, etc.)
 * @param eopdata Matriz de parámetros de orientación terrestre
 * @return Estructura con r2 y v2 (posición y velocidad en t2)
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

