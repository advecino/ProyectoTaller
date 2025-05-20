#ifndef PROYECTOTALLER_RPOLY_H
#define PROYECTOTALLER_RPOLY_H


/**
 * @brief Estructura opaca para el estado interno del cálculo de raíces.
 */
// Opaque state handle
struct RPoly_State;


/**
 * @brief Reserva memoria para el estado necesario para hallar raíces de un polinomio real.
 * @param max_degree Grado máximo del polinomio.
 * @return Puntero al estado reservado.
 */
// Allocate state for finding roots of a real polynomial.
// Degrees up to max_degree are supported by returned state.
struct RPoly_State *real_poly_alloc(int max_degree);


/**
 * @brief Libera el estado previamente reservado.
 * @param s Puntero al estado a liberar.
 */
// Release state
void real_poly_release(struct RPoly_State *s);


/**
 * @brief Calcula las raíces de un polinomio con coeficientes reales.
 *
 * El array `p[]` contiene los coeficientes del polinomio:
 *   p(x) = p[0]*x^degree + ... + p[degree]
 *
 * @param p Coeficientes del polinomio.
 * @param degree Grado del polinomio.
 * @param state Estado interno reservado previamente.
 * @param zeror Array para las partes reales de las raíces.
 * @param zeroi Array para las partes imaginarias de las raíces.
 * @return Número de raíces encontradas (pares en zeror y zeroi).
 */
// Find roots of a polynomial with real coefficients.
// p[] holds coefficients of the polynomial:
//   p(x) = p[0]*x^degree + ... + p[degree]
//
// Caller must have already allocated 'state'.
// Returns number of roots stored in zeror[], zeroi[].
int real_poly_roots_compute(const double p[], int degree,
                            struct RPoly_State *state,
                            double zeror[], double zeroi[]);



/**
 * @brief Función de conveniencia para hallar raíces de un polinomio real.
 *
 * Internamente reserva y libera el estado automáticamente.
 *
 * @param p Coeficientes del polinomio.
 * @param degree Grado del polinomio.
 * @param zeror Array para las partes reales de las raíces.
 * @param zeroi Array para las partes imaginarias de las raíces.
 * @return Número de raíces encontradas.
 */
// Convenience function.
// Find roots of a polynomial with real coefficients.
// p[] holds coefficients of the polynomial:
//   p(x) = p[0]*x^degree + ... + p[degree]
//
// Returns number of roots stored in zeror[], zeroi[]
static inline int real_poly_roots(const double p[], int degree,
                                  double zeror[], double zeroi[])
{
    struct RPoly_State *state = real_poly_alloc(degree);
    int nr = real_poly_roots_compute(p, degree, state, zeror, zeroi);
    real_poly_release(state);

    return nr;
}



#endif //PROYECTOTALLER_RPOLY_H
