/**
 * @file Matrix.cpp
 * @brief Implementación de la clase Matrix para operaciones matriciales.
 * @author advecino
 * @date 20/03/2025
 */
#include "../include/Matrix.h"
#include <iostream>
#include <iomanip>

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
