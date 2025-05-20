#ifndef PROYECTOTALLER_MEASUPDATE_H
#define PROYECTOTALLER_MEASUPDATE_H


#include "Matrix.h"



/**
 * @brief Actualización de medida del filtro de Kalman.
 *
 * @param x Vector de estado (n×1), actualizado en la función.
 * @param P Matriz de covarianza del estado (n×n), actualizada.
 * @param K Matriz de ganancia de Kalman (n×m), devuelta.
 * @param z Vector de medidas (m×1).
 * @param g Vector de medidas predichas (m×1).
 * @param s Vector de desviaciones estándar de medida (m×1).
 * @param G Matriz de sensibilidad de la medida (m×n).
 * @param n Dimensión del estado.
 */
void MeasUpdate(Matrix& x,Matrix& P,Matrix& K, Matrix& z,const Matrix& g,const Matrix& s,const Matrix& G,int n);


#endif //PROYECTOTALLER_MEASUPDATE_H
