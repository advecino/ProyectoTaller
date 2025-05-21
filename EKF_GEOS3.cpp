#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "./include/Sat_const.h"
#include "./include/Mjday.h"
#include "./include/DEInteg.h"
#include "./include/anglesg.h"
#include "./include/Accel.h"
#include "./include/VarEqn.h"
#include "./include/IERS.h"
#include "./include/TimeDiff.h"
#include "./include/gmst.h"
#include "./include/R_z.h"
#include "./include/LTC.h"
#include "./include/AzElPa.h"
#include "./include/TimeUpdate.h"
#include "./include/MeasUpdate.h"
#include "./include/Position.h"
#include "./include/global.h"
#include "./include/Sat_const.h"
#include "./src/cargarPCdesdeFichero.cpp"

using namespace std;



AuxParam AuxParam;


int main() {

    try {
        cargarPC("../data/DE430Coeff.txt");
    } catch (const std::exception& e) {
        std::cerr << "Error al cargar PC: " << e.what() << std::endl;
    }
    try{
        cargarEOP("../data/eop19620101.txt");
    }catch(const std::exception& e){
        std::cerr << "Error al cargar eopdata: " << e.what() << std::endl;
    }

    // Leer coeficientes del modelo gravitacional
    ifstream fin("GGM03S.txt");
    if(!fin) {
        cerr << "Error al abrir GGM03S.txt" << endl;
        return 1;
    }

    for(int n=0; n<=180; n++) {
        for(int m=0; m<=n; m++) {
            int n_read, m_read;
            fin >> n_read >> m_read >> Cnm[n+1][m+1] >> Snm[n+1][m+1];
            // Leer y descartar los últimos 2 valores de la línea
            double temp;
            fin >> temp >> temp;
        }
    }
    fin.close();

    // Inicializar parámetros auxiliares
    AuxParam = {0, 0, 0, 0, false, false, false};


    ifstream eopfin("eop19620101.txt");
    if(!eopfin) {
        cerr << "Error al abrir eop19620101.txt" << endl;
        return 1;
    }

    int eop_count = 0;
    while(eopfin) {
        for(int i=0; i<13; i++) {
            eopfin >> eopdata[i][eop_count];
        }
        eop_count++;
    }
    eopfin.close();

    // Leer observaciones
    const int nobs = 46;
    double obs[nobs][4];

    ifstream obsfin("GEOS3.txt");
    if(!obsfin) {
        cerr << "Error al abrir GEOS3.txt" << endl;
        return 1;
    }

    for(int i=0; i<nobs; i++) {
        string line;
        getline(obsfin, line);

        int Y = stoi(line.substr(0,4));
        int M = stoi(line.substr(5,2));
        int D = stoi(line.substr(8,2));
        int h = stoi(line.substr(12,2));
        int m = stoi(line.substr(15,2));
        double s = stod(line.substr(18,6));
        double az = stod(line.substr(25,8));
        double el = stod(line.substr(35,7));
        double Dist = stod(line.substr(44,10));

        obs[i][0] = Mjday(Y,M,D,h,m,s);
        obs[i][1] = Rad*az;
        obs[i][2] = Rad*el;
        obs[i][3] = 1000.0*Dist;
    }
    obsfin.close();

    // Parámetros de medición
    double sigma_range = 92.5;          // [m]
    double sigma_az = 0.0224*Rad;      // [rad]
    double sigma_el = 0.0139*Rad;      // [rad]

    // Estación Kaena Point
    double lat = Rad*21.5748;          // [rad]
    double lon = Rad*(-158.2706);      // [rad]
    double alt = 300.20;               // [m]

    double Rs[3];
    Position(lon, lat, alt, Rs);

    // Fechas clave para determinación de órbita inicial
    double Mjd1 = obs[0][0];
    double Mjd2 = obs[8][0];
    double Mjd3 = obs[17][0];

    // Determinar órbita inicial
    double r2[3], v2[3];
    anglesg(obs[0][1], obs[8][1], obs[17][1], obs[0][2], obs[8][2], obs[17][2],
            Mjd1, Mjd2, Mjd3, Rs, Rs, Rs, r2, v2);

    double Y0_apr[6] = {r2[0], r2[1], r2[2], v2[0], v2[1], v2[2]};

    double Mjd0 = Mjday(1995,1,29,2,38,0);
    double Mjd_UTC = obs[8][0];

    // Configurar parámetros auxiliares
    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.n = 20;
    AuxParam.m = 20;
    AuxParam.sun = true;
    AuxParam.moon = true;
    AuxParam.planets = true;

    int n_eqn = 6;

    // Integrar hacia atrás hasta la época inicial
    double Y[6];
    memcpy(Y, Y0_apr, 6*sizeof(double));
    DEInteg(Accel, 0, -(obs[8][0]-Mjd0)*86400.0, 1e-13, 1e-6, 6, Y);

    // Matriz de covarianza inicial
    double P[6][6] = {{0}};
    for(int i=0; i<3; i++) P[i][i] = 1e8;
    for(int i=3; i<6; i++) P[i][i] = 1e3;

    // Matriz de transformación topocéntrica

    Matrix LT = LTC(lon, lat) ;

    // Variables para la ecuación de variación
    double yPhi[42] = {0};
    double Phi[6][6] = {{0}};

    // Bucle principal de mediciones
    double t = 0;
    for(int i=0; i<nobs; i++) {
        // Paso anterior
        double t_old = t;
        double Y_old[6];
        memcpy(Y_old, Y, 6*sizeof(double));

        // Incremento de tiempo y propagación
        Mjd_UTC = obs[i][0];
        t = (Mjd_UTC-Mjd0)*86400.0;

        // Obtener parámetros de orientación terrestre
        double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;
        IERSResult iers = IERS(eopdata, Mjd_UTC, 'l');

        double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;
        timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);

        double Mjd_TT = Mjd_UTC + TT_UTC/86400.0;
        double Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;

        AuxParam.Mjd_UTC = Mjd_UTC;
        AuxParam.Mjd_TT = Mjd_TT;

        // Configurar ecuación de variación
        for(int ii=0; ii<6; ii++) {
            yPhi[ii] = Y_old[ii];
            for(int j=0; j<6; j++) {
                if(ii==j) yPhi[6*j+ii+6] = 1.0;
                else yPhi[6*j+ii+6] = 0.0;
            }
        }

        // Integrar ecuación de variación
        DEInteg(VarEqn, 0, t-t_old, 1e-13, 1e-6, 42, yPhi);

        // Extraer matrices de transición de estado
        for(int j=0; j<6; j++) {
            for(int k=0; k<6; k++) {
                Phi[k][j] = yPhi[6*j+k+6];
            }
        }

        // Integrar estado
        DEInteg(Accel, 0, t-t_old, 1e-13, 1e-6, 6, Y);

        // Coordenadas topocéntricas
        double theta = gmst(Mjd_UT1);
        double U[3][3];
        R_z(theta, U);

        double r[3] = {Y[0], Y[1], Y[2]};
        double s[3];
        double Ur[3];
        // Multiplicar U*r
        for(int k=0; k<3; k++) {
            Ur[k] = 0;
            for(int l=0; l<3; l++) {
                Ur[k] += U[k][l] * r[l];
            }
        }
        // Restar Rs y multiplicar por LT
        for(int k=0; k<3; k++) {
            s[k] = 0;
            for(int l=0; l<3; l++) {
                s[k] += LT[k][l] * (Ur[l] - Rs[l]);
            }
        }

        // Actualización de tiempo
        TimeUpdate(P, Phi);

        // Azimut y derivadas parciales
        double Azim, Elev, dAds[3], dEds[3];
        AzElPa(s, Azim, Elev, dAds, dEds);

        double dAdY[6] = {0};
        for(int k=0; k<3; k++) {
            for(int l=0; l<3; l++) {
                dAdY[k] += dAds[l] * LT[l][k];
            }
        }

        // Actualización de medición (azimut)
        double K[6];
        MeasUpdate(Y, obs[i][1], Azim, sigma_az, dAdY, P, 6, K);

        // Elevación y derivadas parciales
        r[0] = Y[0]; r[1] = Y[1]; r[2] = Y[2];
        // Recalcular s
        for(int k=0; k<3; k++) {
            Ur[k] = 0;
            for(int l=0; l<3; l++) {
                Ur[k] += U[k][l] * r[l];
            }
        }
        for(int k=0; k<3; k++) {
            s[k] = 0;
            for(int l=0; l<3; l++) {
                s[k] += LT[k][l] * (Ur[l] - Rs[l]);
            }
        }

        AzElPa(s, Azim, Elev, dAds, dEds);

        double dEdY[6] = {0};
        for(int k=0; k<3; k++) {
            for(int l=0; l<3; l++) {
                dEdY[k] += dEds[l] * LT[l][k];
            }
        }

        // Actualización de medición (elevación)
        MeasUpdate(Y, obs[i][2], Elev, sigma_el, dEdY, P, 6, K);

        // Rango y derivadas parciales
        r[0] = Y[0]; r[1] = Y[1]; r[2] = Y[2];
        // Recalcular s
        for(int k=0; k<3; k++) {
            Ur[k] = 0;
            for(int l=0; l<3; l++) {
                Ur[k] += U[k][l] * r[l];
            }
        }
        for(int k=0; k<3; k++) {
            s[k] = 0;
            for(int l=0; l<3; l++) {
                s[k] += LT[k][l] * (Ur[l] - Rs[l]);
            }
        }

        double Dist = sqrt(s[0]*s[0] + s[1]*s[1] + s[2]*s[2]);
        double dDds[3] = {s[0]/Dist, s[1]/Dist, s[2]/Dist};

        double dDdY[6] = {0};
        for(int k=0; k<3; k++) {
            for(int l=0; l<3; l++) {
                dDdY[k] += dDds[l] * LT[l][k];
            }
        }

        // Actualización de medición (rango)
        MeasUpdate(Y, obs[i][3], Dist, sigma_range, dDdY, P, 6, K);
    }

    // Resultados finales
    double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;
    IERS(eopdata, obs[45][0], 'l', x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC);

    double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;
    timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);

    Mjd_TT = Mjd_UTC + TT_UTC/86400.0;
    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.Mjd_TT = Mjd_TT;

    // Integrar hacia atrás hasta la primera observación
    double Y0[6];
    memcpy(Y0, Y, 6*sizeof(double));
    DEInteg(Accel, 0, -(obs[45][0]-obs[0][0])*86400.0, 1e-13, 1e-6, 6, Y0);

    // Valores verdaderos para comparación
    double Y_true[6] = {
            5753.173e3, 2673.361e3, 3440.304e3,
            4.324207e3, -1.924299e3, -5.728216e3
    };

    // Imprimir resultados
    cout << "\nError de Estimación de Posición" << endl;
    printf("dX%10.1f [m]\n", Y0[0]-Y_true[0]);
    printf("dY%10.1f [m]\n", Y0[1]-Y_true[1]);
    printf("dZ%10.1f [m]\n", Y0[2]-Y_true[2]);

    cout << "\nError de Estimación de Velocidad" << endl;
    printf("dVx%8.1f [m/s]\n", Y0[3]-Y_true[3]);
    printf("dVy%8.1f [m/s]\n", Y0[4]-Y_true[4]);
    printf("dVz%8.1f [m/s]\n", Y0[5]-Y_true[5]);

    // Liberar memoria
    for(int i=0; i<181; i++) {
        free(Cnm[i]);
        free(Snm[i]);
    }
    free(Cnm);
    free(Snm);

    for(int i=0; i<13; i++) {
        free(eopdata[i]);
    }
    free(eopdata);

    return 0;
}