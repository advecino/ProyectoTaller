#ifndef PROYECTOTALLER_TIMEDIFF_H
#define PROYECTOTALLER_TIMEDIFF_H


/**
 * @brief Estructura que contiene diferencias entre escalas temporales.
 */
struct TimeDiffs {
    double UT1_TAI;   ///< Diferencia UT1 - TAI [s]
    double UTC_GPS;   ///< Diferencia UTC - GPS [s]
    double UT1_GPS;   ///< Diferencia UT1 - GPS [s]
    double TT_UTC;    ///< Diferencia TT - UTC [s]
    double GPS_UTC;   ///< Diferencia GPS - UTC [s]
};


/**
 * @brief Calcula diferencias entre escalas de tiempo UT1, UTC, TAI, TT y GPS.
 * @param UT1_UTC Diferencia UT1 - UTC [s].
 * @param TAI_UTC Diferencia TAI - UTC [s].
 * @return Estructura con todas las diferencias de tiempo calculadas.
 */
TimeDiffs timediff(double UT1_UTC, double TAI_UTC);

#endif //PROYECTOTALLER_TIMEDIFF_H
