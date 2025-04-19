#include "../include/Position.h"
#include <cmath>

Matrix Position(double lon, double lat, double h, double R_equ, double f)
{
    // Calcular excentricidad al cuadrado
    const double e2 = f * (2.0 - f);

    // Calcular (co)seno de latitud geodésica
    const double CosLat = cos(lat);
    const double SinLat = sin(lat);

    // Calcular radio de curvatura en el primer vertical
    const double N = R_equ / sqrt(1.0 - e2 * SinLat * SinLat);

    // Calcular componentes del vector de posición
    Matrix r(3, 1); // Vector columna 3x1

    // Componente X (ECEF)
    r(1,1) = (N + h) * CosLat * cos(lon);

    // Componente Y (ECEF)
    r(2,1) = (N + h) * CosLat * sin(lon);

    // Componente Z (ECEF)
    r(3,1) = ((1.0 - e2) * N + h) * SinLat;

    return r;
}