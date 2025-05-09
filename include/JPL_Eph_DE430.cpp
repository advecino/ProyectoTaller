#include "JPL_Eph_DE430.h" // sin hacer

/*
%--------------------------------------------------------------------------
%
% JPL_Eph_DE430: Computes the sun, moon, and nine major planets' equatorial
%                position using JPL Ephemerides
%
% Inputs:
%   Mjd_TDB         Modified julian date of TDB
%
% Output:
%   r_Earth(solar system barycenter (SSB)),r_Mars,r_Mercury,r_Venus,
%   r_Jupiter,r_Saturn,r_Uranus,r_Neptune,r_Pluto,r_Moon,
%   r_Sun(geocentric equatorial position ([m]) referred to the
%   International Celestial Reference Frame (ICRF))
%
% Notes: Light-time is already taken into account
%
% Last modified:   2018/01/11   M. Mahooti
%
%--------------------------------------------------------------------------*/



/*
void JPL_Eph_DE430(double Mjd_TDB,
                   std::vector<double>& r_Earth,
                   std::vector<double>& r_Moon,
                   std::vector<double>& r_Sun,
                   std::vector<double>& r_Mercury,
                   std::vector<double>& r_Venus,
                   std::vector<double>& r_Mars,
                   std::vector<double>& r_Jupiter,
                   std::vector<double>& r_Saturn,
                   std::vector<double>& r_Uranus,
                   std::vector<double>& r_Neptune,
                   std::vector<double>& r_Pluto)
{
    double JD = Mjd_TDB + 2400000.5;

    // Buscar fila i en PC tal que PC[i][0] <= JD <= PC[i][1]
    int i = -1;
    for (int k = 0; k < PC.size(); ++k) {
        if (PC[k][0] <= JD && JD <= PC[k][1]) {
            i = k;
            break;
        }
    }
    if (i == -1) throw std::runtime_error("Fecha fuera del rango de efemérides");

    const auto& PCtemp = PC[i];
    double t1 = PCtemp[0] - 2400000.5;
    double dt = Mjd_TDB - t1;

    // Extraer coeficientes de la Tierra
    std::vector<double> Cx_Earth(13 * 2), Cy_Earth(13 * 2), Cz_Earth(13 * 2);
    int base = 230; // MATLAB: 231, en C++ es 230
    for (int j = 0; j < 13; ++j) {
        Cx_Earth[j] = PCtemp[base + j];
        Cy_Earth[j] = PCtemp[base + 13 + j];
        Cz_Earth[j] = PCtemp[base + 26 + j];
    }
    // Siguientes bloques +39
    for (int j = 0; j < 13; ++j) {
        Cx_Earth[13 + j] = PCtemp[base + 39 + j];
        Cy_Earth[13 + j] = PCtemp[base + 39 + 13 + j];
        Cz_Earth[13 + j] = PCtemp[base + 39 + 26 + j];
    }

    int j_block = (dt <= 16) ? 0 : 1;
    double Mjd0 = t1 + 16 * j_block;


    // Extraer coeficientes del bloque correspondiente y convertir a Matrix 13x1
    std::vector<double> Cx_block(Cx_Earth.begin() + 13 * j_block, Cx_Earth.begin() + 13 * (j_block + 1));
    std::vector<double> Cy_block(Cy_Earth.begin() + 13 * j_block, Cy_Earth.begin() + 13 * (j_block + 1));
    std::vector<double> Cz_block(Cz_Earth.begin() + 13 * j_block, Cz_Earth.begin() + 13 * (j_block + 1));

    Matrix Cx_mat(Cx_block);  // Constructor Matrix(const std::vector<double>&)
    Matrix Cy_mat(Cy_block);
    Matrix Cz_mat(Cz_block);

    Matrix r_Earth_mat = Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0 + 16, Cx_mat, Cy_mat, Cz_mat);

// Convertir Matrix a std::vector<double>
    r_Earth.clear();
    for (int i = 0; i < r_Earth_mat.getFilas(); ++i)
        r_Earth.push_back(r_Earth_mat(i, 0) * 1e3);  // Multiplica aquí si quieres los valores ya en metros

}*/