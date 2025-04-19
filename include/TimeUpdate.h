//
// Created by adria on 19/04/2025.
//

#ifndef PROYECTOTALLER_TIMEUPDATE_H
#define PROYECTOTALLER_TIMEUPDATE_H

#include "Matrix.h"

/**
 * @brief Realiza la actualización temporal de la matriz de covarianza
 * @param P Matriz de covarianza (se actualiza in-place)
 * @param Phi Matriz de transición de estado
 * @param Qdt Puntero a matriz de ruido del proceso (opcional, nullptr para omitir)
 */
void TimeUpdate(Matrix& P, const Matrix& Phi, const Matrix* Qdt = nullptr);




#endif //PROYECTOTALLER_TIMEUPDATE_H
