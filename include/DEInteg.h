#ifndef PROYECTOTALLER_DEINTEG_H
#define PROYECTOTALLER_DEINTEG_H


#include "Matrix.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <functional>

/**
 * @brief Integración numérica de EDOs con método variable de Shampine and Gordon.
 */
class DEInteg {
private:
    // Constantes de estado
    enum DE_STATE {
        DE_INIT = 1,    ///< Reinicio de integración
        DE_DONE = 2,    ///< Paso exitoso
        DE_BADACC = 3,  ///< Precisión no alcanzada
        DE_NUMSTEPS = 4,///< Exceso de pasos permitidos
        DE_STIFF = 5,   ///< Problema rígido sospechado
        DE_INVPARAM = 6 ///< Parámetros inválidos
    };

    // Variables miembro
    DE_STATE State_;
    bool PermitTOUT;
    double told;
    bool OldPermit;
    double delsgn;
    bool start;
    double x;
    Matrix yy;
    double h;
    double hold;
    double hnew;
    int k;
    int kold;
    bool phase1;
    bool nornd;
    bool stiff;
    int nostep;
    int kle4;
    int ifail;
    bool crash;

    // Vectores de trabajo
    Matrix wt;
    Matrix p;
    Matrix yp;
    Matrix phi;
    std::vector<double> g;
    std::vector<double> sig;
    std::vector<double> rho;
    std::vector<double> w;
    std::vector<double> alpha;
    std::vector<double> beta;
    std::vector<double> v;
    std::vector<double> psi_;

    // Constantes
    static const std::vector<double> two;
    static const std::vector<double> gstr;

public:
    /**
     * @brief Constructor por defecto.
     */
    DEInteg();

    /**
     * @brief Ejecuta la integración del sistema de EDOs.
     * @param func Función del sistema dy/dt = f(t, y).
     * @param t Tiempo inicial.
     * @param tout Tiempo final.
     * @param relerr Error relativo permitido.
     * @param abserr Error absoluto permitido.
     * @param y Vector estado inicial (modificado con el resultado final).
     * @return Vector estado final en tout.
     */
    Matrix integrate(std::function<Matrix(double, const Matrix&)> func,
                     double t, double tout, double relerr, double abserr,
                     Matrix& y);
};

#endif //PROYECTOTALLER_DEINTEG_H
