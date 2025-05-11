//
// Created by adria on 11/05/2025.
//

#include <stdexcept>
#include "IERS.h"

/*
%--------------------------------------------------------------------------
%
% IERS: Management of IERS time and polar motion data
%
% Last modified:   2018/02/01   M. Mahooti
%
%--------------------------------------------------------------------------*/




IERSResult IERS(const Matrix& eop, double Mjd_UTC, const std::string& interp) {
    int n = eop.getColumnas();
    if (n < 2) {
        throw std::runtime_error("Matriz EOP debe tener al menos 2 columnas");
    }

    IERSResult res;
    std::string s = "l";
    std::string ss = "n";
    if (interp == s) {
        // Interpolación lineal
        double mjd = floor(Mjd_UTC);
        int i = -1;

        // Encontrar el índice donde mjd coincide con la columna 4 de eop
        for (int j = 0; j < eop.getColumnas(); ++j) {
            if (eop(3, j) == mjd) {  // Columna 4 en MATLAB es índice 3 en C++
                i = j;
                break;
            }
        }

        if (i == -1) {
            throw std::runtime_error("No se encontró el valor MJD en los datos EOP.");
        }

        // Recuperar los valores de la fila i y i+1
        std::vector<double> preeop(12);
        std::vector<double> nexteop(12);

        for (int j = 1; j <= 12; ++j) {
            preeop[j] = eop(i, j);
            nexteop[j] = eop(i + 1, j);
        }

        double mfme = 1440 * (Mjd_UTC - floor(Mjd_UTC));
        double fixf = mfme / 1440;

        // Establecimiento de los parámetros de rotación de la Tierra IERS
        res.x_pole = preeop[4] + (nexteop[4] - preeop[4]) * fixf;
        res.y_pole = preeop[5] + (nexteop[5] - preeop[5]) * fixf;
        res.UT1_UTC = preeop[6] + (nexteop[6] - preeop[6]) * fixf;
        res.LOD = preeop[7] + (nexteop[7] - preeop[7]) * fixf;
        res.dpsi = preeop[8] + (nexteop[8] - preeop[8]) * fixf;
        res.deps = preeop[9] + (nexteop[9] - preeop[9]) * fixf;
        res.dx_pole = preeop[10] + (nexteop[10] - preeop[10]) * fixf;
        res.dy_pole = preeop[11] + (nexteop[11] - preeop[11]) * fixf;
        res.TAI_UTC = preeop[12];

        // Convertir las coordenadas del polo y otros valores
        res.x_pole /= Arcs;  // Coordenada del polo [rad]
        res.y_pole /= Arcs;  // Coordenada del polo [rad]
        res.dpsi /= Arcs;
        res.deps /= Arcs;
        res.dx_pole /= Arcs; // Coordenada del polo [rad]
        res.dy_pole /= Arcs; // Coordenada del polo [rad]
    } else if (interp == ss) {
        double mjd = floor(Mjd_UTC);
        int i = -1;

        // Encontrar el índice donde mjd coincide con la columna 4 de eop
        for (int j = 0; j < eop.getColumnas(); ++j) {
            if (eop(3, j) == mjd) {
                i = j;
                break;
            }
        }

        if (i == -1) {
            throw std::runtime_error("No se encontró el valor MJD en los datos EOP.");
        }


        std::vector<double> eop_row(12);


        for (int j = 1; j <= 12; ++j) {
            eop_row[j] = eop(i, j);
        }

        // Establecimiento de los parámetros de rotación de la Tierra IERS
        res.x_pole = eop_row[4] / Arcs;  // Coordenada del polo [rad]
        res.y_pole = eop_row[5] / Arcs;  // Coordenada del polo [rad]
        res.UT1_UTC = eop_row[6];                 // Diferencia de tiempo UT1-UTC [s]
        res.LOD = eop_row[7];                 // Longitud del día [s]
        res.dpsi = eop_row[8] / Arcs;
        res.deps = eop_row[9] / Arcs;
        res.dx_pole = eop_row[10] / Arcs; // Coordenada del polo [rad]
        res.dy_pole = eop_row[11] / Arcs; // Coordenada del polo [rad]
        res.TAI_UTC = eop_row[12];               // Diferencia de tiempo TAI-UTC [s]
    }
    return res;
}