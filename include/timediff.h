//
// Created by adria on 11/05/2025.
//

#ifndef PROYECTOTALLER_TIMEDIFF_H
#define PROYECTOTALLER_TIMEDIFF_H


/**
 * @brief Calculates various time differences between time standards
 *
 * @param UT1_UTC Difference between UT1 and UTC [s]
 * @param TAI_UTC Difference between TAI and UTC [s]
 * @return struct TimeDiffs containing all computed time differences
 */
struct TimeDiffs {
    double UT1_TAI;
    double UTC_GPS;
    double UT1_GPS;
    double TT_UTC;
    double GPS_UTC;
};

TimeDiffs timediff(double UT1_UTC, double TAI_UTC);

#endif //PROYECTOTALLER_TIMEDIFF_H
