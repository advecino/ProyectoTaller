/**
 * @file Matrix.h
 * @brief Declaración de la clase Matrix para operaciones matriciales.
 * @author advecino
 * @date 20/03/2025
 */

#ifndef PROYECTOTALLER_MATRIX_H
#define PROYECTOTALLER_MATRIX_H

/**
 * @class Matrix
 * @brief Clase para representar y manipular matrices.
 */

class Matrix {
public:
    /**
    * @brief Constructor de la matriz con dimensiones especificadas.
    * @param fil Número de filas.
    * @param col Número de columnas.
    */
    Matrix(int fil, int col);

    /**
     * @brief Constructor que inicializa la matriz con valores dados.
     * @param fil Número de filas.
     * @param col Número de columnas.
     * @param v Array de valores para inicializar la matriz.
     * @param n Número de elementos en el array.
     */
    Matrix(int fil, int col, double v[],int n);

    /**
    * @brief Constructor de copia.
    * @param m Objeto Matrix a copiar.
    */
    Matrix(const Matrix& m);

    /**
    * @brief Destructor de la matriz.
    */
    ~Matrix();


    /**
    * @brief Sobrecarga del operador de asignación.
    * @param matrix2 Matriz a asignar.
    * @return Referencia a la matriz resultante.
    */
    Matrix& operator= (const Matrix& matrix2);
    /**
    * @brief Sobrecarga del operador de suma.
    * @param matrix2 Matriz a sumar.
    * @return Matriz resultante de la suma.
    */
    Matrix operator+ (const Matrix& matrix2);
    /**
     * @brief Sobrecarga del operador de resta.
     * @param matrix2 Matriz a restar.
     * @return Matriz resultante de la resta.
     */
    Matrix operator- (const Matrix& matrix2);
    /**
    * @brief Sobrecarga del operador de multiplicación.
    * @param matrix2 Matriz a multiplicar.
    * @return Matriz resultante de la multiplicación.
    */
    Matrix operator* (const Matrix& matrix2);
    /**
     * @brief Sobrecarga del operador de acceso a elementos de la matriz.
     * @param i Índice de la fila.
     * @param j Índice de la columna.
     * @return Referencia al valor en la posición (i, j).
     */
    double& operator()(const int i, const int j) const;

    /**
    * @brief Imprime la matriz en la consola.
    */
    void print();

private:
    /**
     * @brief Inicializa la matriz con valores en cero.
     */
    void initMatrix();

private:
    int fil;  ///< Número de filas de la matriz.
    int col; ///< Número de columnas de la matriz.
    double** matrix; ///< Puntero a la memoria donde se almacenan los valores de la matriz.

};


#endif //PROYECTOTALLER_MATRIX_H
