#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include "../include/global.h"


void cargarPC(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("No se pudo abrir el fichero PC");
    }

    std::vector<std::vector<double>> rows;
    std::string line;

    while (std::getline(file, line)) {
        std::replace(line.begin(), line.end(), ',', '.');

        line.erase(std::remove_if(line.begin(), line.end(), [](unsigned char c) {
            return !std::isprint(c) && !std::isspace(c);
        }), line.end());

        std::stringstream ss(line);
        std::vector<double> row;
        std::string token;

        while (ss >> token) {
            try {
                row.push_back(std::stod(token));
            } catch (...) {
                std::cerr << "Error convirtiendo: [" << token << "]" << std::endl;
                throw;
            }
        }

        if (!row.empty()) {
            rows.push_back(row);
        }
    }

    file.close();

    int filas = rows.size();
    int cols = rows[0].size();
    PC = Matrix(filas, cols);

    for (int i = 1; i <= filas; ++i) {
        for (int j = 1; j <= cols; ++j) {
            PC(i, j) = rows[i - 1][j - 1];
        }
    }
}

/**
 * @brief Carga la matriz eopdata desde un archivo .txt con formato europeo.
 * @param filepath Ruta al archivo eop19620101.txt
 */
void cargarEOP(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("No se pudo abrir el fichero EOP");
    }

    std::vector<std::vector<double>> columnas(13); // 13 filas → cada vector representa una fila

    std::string line;
    int lineCount = 0;

    while (std::getline(file, line)) {
        ++lineCount;
        std::replace(line.begin(), line.end(), ',', '.');

        line.erase(std::remove_if(line.begin(), line.end(), [](unsigned char c) {
            return !std::isprint(c) && !std::isspace(c);
        }), line.end());

        std::stringstream ss(line);
        std::string token;
        std::vector<double> values;

        while (ss >> token) {
            try {
                values.push_back(std::stod(token));
            } catch (...) {
                std::cerr << "Error en línea " << lineCount << " al convertir: [" << token << "]" << std::endl;
                throw;
            }
        }

        if (values.size() != 13) {
            throw std::runtime_error("Línea " + std::to_string(lineCount) + ": se esperaban 13 valores.");
        }

        for (int i = 0; i < 13; ++i) {
            columnas[i].push_back(values[i]);
        }
    }

    file.close();

    int filas = 13;
    int cols = columnas[0].size();
    eopdata = Matrix(filas, cols);

    for (int i = 1; i <= filas; ++i) {
        for (int j = 1; j <= cols; ++j) {
            eopdata(i, j) = columnas[i - 1][j - 1];
        }
    }

}
