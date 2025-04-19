//
// Created by adria on 20/04/2025.
//

#ifndef PROYECTOTALLER_MEANOBLIQUITY_H
#define PROYECTOTALLER_MEANOBLIQUITY_H


/**
 * @brief Calcula la oblicuidad media de la eclíptica
 * @param Mjd_TT Fecha Juliana Modificada (Tiempo Terrestre)
 * @param MJD_J2000 Valor de MJD para la época J2000
 * @param Rad Factor de conversión de grados a radianes
 * @return Oblicuidad media de la eclíptica en radianes
 */
double MeanObliquity(double Mjd_TT, double MJD_J2000, double Rad);


#endif //PROYECTOTALLER_MEANOBLIQUITY_H
