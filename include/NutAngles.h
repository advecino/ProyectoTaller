//
// Created by adria on 20/04/2025.
//

#ifndef PROYECTOTALLER_NUTANGLES_H
#define PROYECTOTALLER_NUTANGLES_H

/**
 * @brief Calcula los ángulos de nutación (IAU 1980)
 * @param Mjd_TT Fecha Juliana Modificada (Tiempo Terrestre)
 * @param dpsi Nutación en longitud [rad] (salida)
 * @param deps Nutación en oblicuidad [rad] (salida)
 */
void NutAngles(double Mjd_TT, double& dpsi, double& deps);


#endif //PROYECTOTALLER_NUTANGLES_H
