//
// Created by adria on 12/05/2025.
//

#ifndef PROYECTOTALLER_VAREQN_H
#define PROYECTOTALLER_VAREQN_H


#include "Matrix.h"
#include "AuxParam.h"  // <<< sólo aquí

/**
 * @brief Derivada conjunta de estado y matriz de transición.
 * @param t Tiempo transcurrido desde la época inicial (en segundos).
 * @param yPhi Vector de tamaño 42 con el estado y la matriz de transición.
 * @param params Parámetros auxiliares (fecha, grado y orden armónico, etc.).
 * @param eopdata Matriz con datos EOP del IERS para correcciones.
 * @return Vector de derivadas (42×1), es decir, d[yPhi]/dt.
 */
Matrix VarEqn(double t,const Matrix& yPhi,const AuxParam& params,Matrix& eopdata);


#endif //PROYECTOTALLER_VAREQN_H

