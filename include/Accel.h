//
// Created by advecino on 03/04/2025.
//

#ifndef PROYECTOTALLER_ACCEL_H
#define PROYECTOTALLER_ACCEL_H


#include <vector>
#include <tuple>

typedef struct {
    double MJD_J2000;
    double GM_Sun, GM_Moon, GM_Mercury, GM_Venus, GM_Mars, GM_Jupiter;
    double GM_Saturn, GM_Uranus, GM_Neptune, GM_Pluto;
} Constants;

typedef struct {
    double Mjd_UTC;
    int n, m;
    bool sun, moon, planets;
} AuxParams;

// Variables globales externas
extern Constants consts;
extern AuxParams AuxParam;
extern std::vector<double> eopdata;

// Función principal para calcular la aceleración
std::vector<double> Accel(double x, const std::vector<double>& Y);

// Funciones auxiliares
std::tuple<double, double, double, double, double, double, double, double, double> IERS(const std::vector<double>& eopdata, double mjd, char option);
std::tuple<double, double, double, double, double> timediff(double UT1_UTC, double TAI_UTC);
double Mjday_TDB(double Mjd_TT);
std::vector<double> AccelHarmonic(const std::vector<double>& r, const std::vector<std::vector<double>>& E, int n, int m);
std::vector<double> AccelPointMass(const std::vector<double>& r, const std::vector<double>& r_body, double GM);


#endif //PROYECTOTALLER_ACCEL_H
