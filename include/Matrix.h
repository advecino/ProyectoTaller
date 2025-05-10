/**
 * @file Matrix.h
 * @brief Declaración de la clase Matrix para operaciones matriciales.
 * @author advecino
 * @date 20/03/2025
 */

#ifndef PROYECTOTALLER_MATRIX_H
#define PROYECTOTALLER_MATRIX_H

#include <vector>

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

    explicit Matrix(const std::vector<double>& values);

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
    Matrix operator* (const Matrix& matrix2) const;
    Matrix& operator/=(double scalar);
    /**
     * @brief Sobrecarga del operador de acceso a elementos de la matriz.
     * @param i Índice de la fila.
     * @param j Índice de la columna.
     * @return Referencia al valor en la posición (i, j).
     */
    double& operator()( int i,  int j) const;













    /**
     * @brief Suma un escalar a la matriz (modifica la matriz actual)
     * @param scalar El valor escalar a sumar
     * @return Referencia a la matriz resultante
     */
    Matrix& operator+=(double scalar);

    /**
     * @brief Resta un escalar a la matriz (modifica la matriz actual)
     * @param scalar El valor escalar a restar
     * @return Referencia a la matriz resultante
     */
    Matrix& operator-=(double scalar);
    /**
    * @brief Multiplica la matriz por un escalar
    * @param scalar El valor escalar a multiplicar
    * @return Referencia a la matriz resultante (para encadenar operaciones)
    */
    Matrix& operator*=(double scalar);
    /**
     * @brief Versiones no-miembro para suma/resta con escalar
     */



    friend Matrix operator+(const Matrix& m, double scalar);
    friend Matrix operator+(double scalar, const Matrix& m);
    friend Matrix operator-(const Matrix& m, double scalar);
    friend Matrix operator-(double scalar, const Matrix& m);
    friend Matrix operator*(const Matrix& m, double scalar);
    friend Matrix operator*(double scalar, const Matrix& m);
    friend Matrix operator/(const Matrix& m, double scalar);
    friend Matrix operator/(double scalar, const Matrix& m);





    /**
    * @brief Imprime la matriz en la consola.
    */
    void print() const;
    double norm() const;




    int getFilas() const;
    int getColumnas() const;
    Matrix transpuesta() const;
    Matrix inversa() const;


    /**
     * @brief Calcula el producto punto (dot product) entre dos vectores
     * @param a Primer vector (debe ser matriz nx1 o 1xn)
     * @param b Segundo vector (debe ser matriz nx1 o 1xn)
     * @return Resultado del producto punto
     * @throw std::invalid_argument Si las dimensiones no son compatibles
     */
    static double dot(const Matrix& a, const Matrix& b);

    /**
     * @brief Calcula el producto cruz (cross product) entre dos vectores 3D
     * @param a Primer vector (debe ser matriz 3x1 o 1x3)
     * @param b Segundo vector (debe ser matriz 3x1 o 1x3)
     * @return Matriz resultante del producto cruz (3x1)
     * @throw std::invalid_argument Si los vectores no son 3D
     */
    static Matrix cross(const Matrix& a, const Matrix& b);


private:
    /**
     * @brief Inicializa la matriz con valores en cero.
     */
    void initMatrix();
    int fil;  ///< Número de filas de la matriz.
    int col; ///< Número de columnas de la matriz.
    double** matrix; ///< Puntero a la memoria donde se almacenan los valores de la matriz.

};


#endif //PROYECTOTALLER_MATRIX_H
