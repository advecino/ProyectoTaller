#include "../include/TimeDiff.h"
/*
%--------------------------------------------------------------------------
%
% Time differences [s]
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------*/


TimeDiffs timediff(double UT1_UTC, double TAI_UTC) {
    const double TT_TAI = 32.184;   // TT-TAI time difference [s]

    const double GPS_TAI = -19.0;   // GPS-TAI time difference [s]

    TimeDiffs diffs;

    double TT_GPS  =  TT_TAI-GPS_TAI;  // TT-GPS time difference [s]

    double TAI_GPS = -GPS_TAI;         // TAI-GPS time difference [s]

    diffs.UT1_TAI = UT1_UTC-TAI_UTC;  // UT1-TAI time difference [s]

    double UTC_TAI = -TAI_UTC;         // UTC-TAI time difference [s]

    diffs.UTC_GPS = UTC_TAI-GPS_TAI;  // UTC_GPS time difference [s]

    diffs.UT1_GPS = diffs.UT1_TAI-GPS_TAI;  // UT1-GPS time difference [s]

    diffs.TT_UTC  = TT_TAI-UTC_TAI;   //  TT-UTC time difference [s]

    diffs.GPS_UTC = GPS_TAI-UTC_TAI;  // GPS-UTC time difference [s]

    return diffs;
}
