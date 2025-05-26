#ifndef PROYECTOTALLER_CARGARDESDEFICHERO_H
#define PROYECTOTALLER_CARGARDESDEFICHERO_H

#include <string>


/**
 * @brief Carga la matriz PC desde un archivo .txt
 * @param filepath Ruta al archivo DE430Coeff.txt
 */
void cargarPC(const std::string& filepath);
/**
 * @brief Carga la matriz eopdata desde un archivo .txt
 * @param filepath Ruta al archivo eop19620101.txt
 */
void cargarEOP(const std::string& filepath);

#endif //PROYECTOTALLER_CARGARDESDEFICHERO_H
