/**
 * @file Mjday_TDB.h
 * @brief Calcula la Fecha Juliana Modificada para tiempo dinámico baricéntrico (TDB)
 */


#ifndef PROYECTOTALLER_MJDAY_TDB_H
#define PROYECTOTALLER_MJDAY_TDB_H


/**
 * @brief Calcula el Modified Julian Date para TDB
 * @param Mjd_TT Fecha Juliana Modificada (Tiempo Terrestre)
 * @return Mjd_TDB Fecha Juliana Modificada (Tiempo Dinámico Baricéntrico)
 */
double Mjday_TDB(double Mjd_TT);


#endif //PROYECTOTALLER_MJDAY_TDB_H
