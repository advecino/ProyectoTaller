//
// Created by adria on 08/05/2025.
//

#ifndef PROYECTOTALLER_JPL_EPH_DE430_H
#define PROYECTOTALLER_JPL_EPH_DE430_H

#include <vector>


void JPL_Eph_DE430(
        double Mjd_TDB,
        std::vector<double>& r_Mercury,
        std::vector<double>& r_Venus,
        std::vector<double>& r_Earth,
        std::vector<double>& r_Mars,
        std::vector<double>& r_Jupiter,
        std::vector<double>& r_Saturn,
        std::vector<double>& r_Uranus,
        std::vector<double>& r_Neptune,
        std::vector<double>& r_Pluto,
        std::vector<double>& r_Moon,
        std::vector<double>& r_Sun
);


#endif //PROYECTOTALLER_JPL_EPH_DE430_H
