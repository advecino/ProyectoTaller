/**
 * @file Matrix.cpp
 * @brief Implementación de la clase Matrix para operaciones matriciales.
 * @author advecino
 * @date 20/03/2025
 */
#include "../include/Matrix.h"
#include <iostream>
#include <iomanip>
#include <cmath>

/**
 * @brief Constructor de la matriz con dimensiones especificadas.
 * @param fil Número de filas.
 * @param col Número de columnas.
 */
Matrix::Matrix(int fil, int col) : fil(fil), col(col)
{
    initMatrix();
}

/**
 * @brief Constructor que inicializa la matriz con valores dados.
 * @param fil Número de filas.
 * @param col Número de columnas.
 * @param v Array de valores para inicializar la matriz.
 * @param n Número de elementos en el array.
 */
Matrix::Matrix(int fil, int col, double v[], int n): fil(fil), col(col)
{
    initMatrix();

    int k = 0;

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++){
            if (k < n)
                matrix[i][j] = v[k++];
            else
                matrix[i][j] = 0;
        }
}

/**
 * @brief Constructor de copia.
 * @param m Objeto Matrix a copiar.
 */
Matrix::Matrix(const Matrix& m)
{
    *this = m;
}


/**
 * @brief Inicializa la matriz con valores en cero.
 */
void Matrix::initMatrix()
{
    matrix = new double*[fil];
    for (int i = 0; i < fil; i++)
        matrix[i] = new double[col];

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            matrix[i][j] = 0.0;
}

/**
 * @brief Destructor de la matriz.
 */
Matrix::~Matrix()
{
    for (int i = 0; i < fil; i++)
        delete[] matrix[i];

    delete[] matrix;
}

/**
 * @brief Sobrecarga del operador de asignación.
 * @param matrix2 Matriz a asignar.
 * @return Referencia a la matriz resultante.
 */
Matrix& Matrix::operator=(const Matrix& matrix2)
{
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            this->matrix[i][j] = matrix2.matrix[i][j];

    return *this;
}

/**
 * @brief Sobrecarga del operador de suma.
 * @param matrix2 Matriz a sumar.
 * @return Matriz resultante de la suma.
 */
Matrix Matrix::operator+(const Matrix& matrix2)
{
    Matrix result(fil, col);

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] + matrix2.matrix[i][j];

    return result;
}

/**
 * @brief Sobrecarga del operador de resta.
 * @param matrix2 Matriz a restar.
 * @return Matriz resultante de la resta.
 */
Matrix Matrix::operator-(const Matrix& matrix2)
{
    Matrix result(fil, col);

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] - matrix2.matrix[i][j];

    return result;
}

/**
 * @brief Sobrecarga del operador de multiplicación.
 * @param matrix2 Matriz a multiplicar.
 * @return Matriz resultante de la multiplicación.
 */
Matrix Matrix::operator*(const Matrix& matrix2)
{
    Matrix result(fil, col);

    for (int i = 0; i < this->fil ; i++){
        for (int j = 0; j < matrix2.col; j++){
            result.matrix[i][j] = 0;
            for (int k = 0; k < this->col; k++){
                result.matrix[i][j] = result.matrix[i][j] + this->matrix[i][k] * matrix2.matrix[k][j];
            }
        }
    }

    return result;
}

/**
 * @brief Sobrecarga del operador de acceso a elementos de la matriz.
 * @param i Índice de la fila.
 * @param j Índice de la columna.
 * @return Referencia al valor en la posición (i, j).
 */
double& Matrix::operator()(const int i, const int j) const
{
    return matrix[i-1][j-1];
}

/**
 * @brief Imprime la matriz en la consola.
 */
void Matrix::print()
{
    for (int i = 0; i < fil; i++){
        for (int j = 0; j < col; j++){
            std::cout << std::fixed << std::setprecision(14) << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}



/**
 * @brief Obtiene el número de filas de la matriz
 * @return Número de filas
 */
int Matrix::getFilas() const {
    return fil;
}

/**
 * @brief Obtiene el número de columnas de la matriz
 * @return Número de columnas
 */
int Matrix::getColumnas() const {
    return col;
}

/**
 * @brief Calcula la matriz transpuesta
 * @return Nueva matriz que es la transpuesta de la original
 */
Matrix Matrix::transpuesta() const {
    Matrix result(fil, col);
    for(int i = 1; i <= fil; ++i) {
        for(int j = 1; j <= col; ++j) {
            result(j,i) = (*this)(i,j);
        }
    }
    return result;
}


/**
 * @brief Calcula la matriz inversa (solo para matrices 1x1 y 2x2)
 * @return Matriz inversa
 * @throw std::runtime_error Si la matriz no es cuadrada o es singular
 */
Matrix Matrix::inversa() const {
    if(fil != col) {
        throw std::runtime_error("Matrix no es cuadrada");
    }

    // Caso 1x1
    if(fil == 1) {
        if(fabs((*this)(1,1)) < 1e-15) {
            throw std::runtime_error("Matrix es singular (determinante cero)");
        }
        Matrix inv(1,1);
        inv(1,1) = 1.0 / (*this)(1,1);
        return inv;
    }
        // Caso 2x2
    else if(fil == 2) {
        double det = (*this)(1,1)*(*this)(2,2) - (*this)(1,2)*(*this)(2,1);
        if(fabs(det) < 1e-15) {
            throw std::runtime_error("Matrix es singular (determinante cero)");
        }

        Matrix inv(2,2);
        inv(1,1) = (*this)(2,2)/det;
        inv(1,2) = -(*this)(1,2)/det;
        inv(2,1) = -(*this)(2,1)/det;
        inv(2,2) = (*this)(1,1)/det;
        return inv;
    }
        // Matrices más grandes no soportadas
    else {
        throw std::runtime_error("Inversa no implementada para matrices >2x2");
    }
}

double Matrix::norm() const {
    if (col != 1) {
        throw std::invalid_argument("norm() solo aplicable a vectores columna");
    }

    double sum = 0.0;
    for (int i = 0; i < fil; ++i) {
        sum += matrix[i][0] * matrix[i][0];
    }
    return sqrt(sum);
}