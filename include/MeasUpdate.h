//
// Created by adria on 20/04/2025.
//

#ifndef PROYECTOTALLER_MEASUPDATE_H
#define PROYECTOTALLER_MEASUPDATE_H


#include "Matrix.h"



/**
 * @brief Performs the Kalman measurement update.
 *
 * @param x     State vector (n×1), updated in place.
 * @param P     State covariance matrix (n×n), updated in place.
 * @param K     Kalman gain matrix (n×m), returned.
 * @param z     Measurement vector (m×1).
 * @param g     Predicted measurement vector (m×1).
 * @param s     Measurement standard-deviation vector (m×1).
 * @param G     Measurement sensitivity matrix (m×n).
 * @param n     Dimension of the state.
 */
void MeasUpdate(
        Matrix& x,
        Matrix& P,
        Matrix& K,
        Matrix& z,
        const Matrix& g,
        const Matrix& s,
        const Matrix& G,
        int n
);


#endif //PROYECTOTALLER_MEASUPDATE_H
