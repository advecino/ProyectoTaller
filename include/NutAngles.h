#ifndef PROYECTOTALLER_NUTANGLES_H
#define PROYECTOTALLER_NUTANGLES_H

/**
 * @brief Ángulos de nutación según el modelo IAU 1980.
 * @param Mjd_TT Fecha juliana modificada (Tiempo Terrestre).
 * @param dpsi Nutación en longitud [rad] (salida).
 * @param deps Nutación en oblicuidad [rad] (salida).
 */
void NutAngles(double Mjd_TT, double& dpsi, double& deps);


#endif //PROYECTOTALLER_NUTANGLES_H
