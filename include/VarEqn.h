//
// Created by adria on 12/05/2025.
//

#ifndef PROYECTOTALLER_VAREQN_H
#define PROYECTOTALLER_VAREQN_H


#include "Matrix.h"
#include "AuxParam.h"  // <<< sólo aquí

/**
 * @brief Derivada conjunta de estado y matriz de transición.
 * @param t     Tiempo desde epoch [s].
 * @param yPhi  (42×1) [r; v; vec(Φ)].
 * @param params Auxiliares (fechas, grado/orden).
 * @param eop    EOP data para IERS.
 * @return yPhip (42×1).
 */
Matrix VarEqn(
        double t,
        const Matrix& yPhi,
        const AuxParam& params,
        Matrix& eopdata
);


#endif //PROYECTOTALLER_VAREQN_H

