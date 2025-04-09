//
// Created by advecino on 09/04/2025.
//
/**
 * @file GMST.h
 * @brief Declaración del cálculo del Tiempo Sidéreo Medio de Greenwich
 * @author M. Mahooti (adaptado por [tu nombre])
 * @date 20/03/2025
 */
#ifndef PROYECTOTALLER_GMST_H
#define PROYECTOTALLER_GMST_H


/**
 * @brief Calcula el Tiempo Sidéreo Medio de Greenwich
 * @param Mjd_UT1 Fecha Juliana Modificada UT1
 * @return TSG en radianes [0..2pi]
 */
double gmst(double Mjd_UT1);


#endif //PROYECTOTALLER_GMST_H
