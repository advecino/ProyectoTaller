#ifndef PROYECTOTALLER_TIMEUPDATE_H
#define PROYECTOTALLER_TIMEUPDATE_H

#include "Matrix.h"

/**
 * @brief Actualiza la matriz de covarianza mediante la ecuación de predicción.
 * @param P Matriz de covarianza del estado (se modifica directamente).
 * @param Phi Matriz de transición de estado.
 * @param Qdt Puntero a la matriz de covarianza del ruido del proceso (opcional, puede ser nullptr).
 */
void TimeUpdate(Matrix& P, const Matrix& Phi, const Matrix* Qdt = nullptr);

#endif //PROYECTOTALLER_TIMEUPDATE_H
