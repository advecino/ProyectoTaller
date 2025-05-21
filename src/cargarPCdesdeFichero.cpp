#include <fstream>
#include <sstream>
#include <iostream>
#include "../include/global.h"

/**
 * @brief Carga la matriz PC desde un archivo .txt con formato europeo
 * @param filepath Ruta al archivo DE430Coeff.txt
 */
void cargarPC(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("No se pudo abrir el fichero PC");
    }

    std::vector<std::vector<double>> rows;
    std::string line;

    while (std::getline(file, line)) {
        // Reemplazar comas por puntos (para convertir decimales europeos)
        std::replace(line.begin(), line.end(), ',', '.');

        // Eliminar caracteres no imprimibles (incluido BOM o espacios raros)
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

    std::cout << "Matriz PC cargada: " << filas << "x" << cols << std::endl;
}
