#include "../include/Position.h"
#include <cmath>

Matrix Position(double lon, double lat, double h, double R_equ, double f)
{
    const double e2 = f * (2.0 - f);

    const double CosLat = cos(lat);
    const double SinLat = sin(lat);

    const double N = R_equ / sqrt(1.0 - e2 * SinLat * SinLat);


    Matrix r(3, 1);

    r(1,1) = (N + h) * CosLat * cos(lon);

    r(2,1) = (N + h) * CosLat * sin(lon);

    r(3,1) = ((1.0 - e2) * N + h) * SinLat;

    return r;
}