#ifndef PROYECTOTALLER_JPL_EPH_DE430_H
#define PROYECTOTALLER_JPL_EPH_DE430_H

#include <vector>
#include "Matrix.h"
#include "global.h"
#include "Cheb3D.h"


void JPL_Eph_DE430(double Mjd_TDB,
                   Matrix& r_Mercury, Matrix& r_Venus, Matrix& r_Earth,
                   Matrix& r_Mars, Matrix& r_Jupiter, Matrix& r_Saturn,
                   Matrix& r_Uranus, Matrix& r_Neptune, Matrix& r_Pluto,
                   Matrix& r_Moon, Matrix& r_Sun);


#endif //PROYECTOTALLER_JPL_EPH_DE430_H
