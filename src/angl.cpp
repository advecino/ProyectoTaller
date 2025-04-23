/*%--------------------------------------------------------------------------
%
%  inputs:
%    vec1         - vector 1
%    vec2         - vector 2
%
%  output:
%    theta        - angle between the two vectors  -pi to pi
%
%--------------------------------------------------------------------------*/
/**
 * @file angl.cpp
 * @brief Implementación de la función angl
 * @author [Tu Nombre]
 * @date [Fecha]
 */

#include "../include/angl.h"
#include <cmath>

double angl(const Matrix& vec1, const Matrix& vec2)
{
    const double small = 0.00000001;
    const double undefined = 999999.1;

    double magv1 = 0.0;
    double magv2 = 0.0;
    double dot_product = 0.0;

    for(int i = 1; i <= 3; ++i) {
        magv1 += vec1(i,1) * vec1(i,1);
        magv2 += vec2(i,1) * vec2(i,1);
        dot_product += vec1(i,1) * vec2(i,1);
    }

    magv1 = sqrt(magv1);
    magv2 = sqrt(magv2);

    if(magv1 * magv2 > small * small) {
        double temp = dot_product / (magv1 * magv2);

        if(temp > 1.0) temp = 1.0;
        if(temp < -1.0) temp = -1.0;

        return acos(temp);
    }
    else {
        return undefined;
    }
}