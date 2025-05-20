#ifndef PROYECTOTALLER_SAT_CONST_H
#define PROYECTOTALLER_SAT_CONST_H

/// @name Constantes matemáticas
/// @{
const double pi        = 3.14159265358979324;
const double pi2       = 2.0 * pi;               ///< 2π
const double Rad       = pi / 180.0;             ///< Radianes por grado
const double Deg       = 180.0 / pi;             ///< Grados por radián
const double Arcs      = 3600.0 * 180.0 / pi;     ///< Arcosegundos por radián
/// @}

/// @name Constantes generales
/// @{
const double MJD_J2000 = 51544.5;                ///< Fecha juliana modificada para la época J2000.0
const double T_B1950   = -0.500002108;           ///< Época B1950
const double c_light   = 299792457.999999984;    ///< Velocidad de la luz [m/s]
const double AU        = 149597870659.999996;    ///< Unidad astronómica [m]
/// @}

/// @name Parámetros físicos de la Tierra, Sol y Luna
/// @{
const double R_Earth   = 6378.137e3;             ///< Radio ecuatorial de la Tierra [m]
const double f_Earth   = 1.0 / 298.257223563;    ///< Achatamiento terrestre
const double R_Sun     = 696000.0e3;             ///< Radio del Sol [m]
const double R_Moon    = 1738.0e3;               ///< Radio de la Luna [m]
/// @}

/// @name Rotación terrestre
/// @{
const double omega_Earth = 7.2921158553e-5;      ///< Velocidad angular de la Tierra [rad/s]
/// @}

/// @name Coeficientes gravitacionales [m³/s²]
/// @{
const double GM_Earth    = 398600.4418e9;                     ///< Tierra
const double GM_Sun      = 1.327124399354841e20;              ///< Sol
const double GM_Moon     = GM_Earth / 81.3005869999999931;    ///< Luna
const double GM_Mercury  = 22032.08047272131e9;
const double GM_Venus    = 324858.7656168717e9;
const double GM_Mars     = 42828.28658876890e9;
const double GM_Jupiter  = 126712597.0817946e9;
const double GM_Saturn   = 37939519.70882996e9;
const double GM_Uranus   = 5780158.533597719e9;
const double GM_Neptune  = 6871307.771479524e9;
const double GM_Pluto    = 1020.864920706286e9;
/// @}

/// @name Presión de radiación solar
/// @{
const double P_Sol       = 4.560e-6;             ///< A 1 AU [N/m²]
/// @}

#endif //PROYECTOTALLER_SAT_CONST_H
