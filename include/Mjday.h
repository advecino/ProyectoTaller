#ifndef PROYECTOTALLER_MJDAY_H
#define PROYECTOTALLER_MJDAY_H

/**
 * @brief Convierte fecha y hora UTC a fecha juliana modificada (MJD).
 * @param yr Año.
 * @param mon Mes.
 * @param day Día.
 * @param hr Hora (opcional, por defecto 0).
 * @param min Minutos (opcional, por defecto 0).
 * @param sec Segundos (opcional, por defecto 0.0).
 * @return Fecha juliana modificada (MJD).
 */
double Mjday(int yr, int mon, int day, int hr = 0, int min = 0, double sec = 0.0);


#endif //PROYECTOTALLER_MJDAY_H
