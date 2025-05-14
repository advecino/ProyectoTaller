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
Matrix::Matrix(const Matrix& m) : fil(m.fil), col(m.col) {
    matrix = new double*[fil];
    for (int i = 0; i < fil; i++) {
        matrix[i] = new double[col];
        for (int j = 0; j < col; j++) {
            matrix[i][j] = m.matrix[i][j];
        }
    }
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
 * @brief Constructor que crea una matriz de 1 fila a partir de un vector de valores.
 * @param values Vector de valores a usar como fila.
 */
Matrix::Matrix(const std::vector<double>& values) : fil(1), col(values.size())
{
    initMatrix();
    for (int j = 0; j < col; ++j) {
        matrix[0][j] = values[j];
    }
}


Matrix::Matrix(int fil, int col, double v[]) : fil(fil), col(col) {
    initMatrix();

    int k = 0;
    for (int i = 0; i < fil; i++) {
        for (int j = 0; j < col; j++) {
            matrix[i][j] = v[k++];
        }
    }
}

Matrix::Matrix(const std::initializer_list<double>& values, int cols) {
    int size = values.size();
    if (size % cols != 0) {
        throw std::invalid_argument("Number of values must be divisible by columns");
    }
    fil = size / cols;
    col = cols;
    matrix = new double*[fil];
    for (int i = 0; i < fil; i++) {
        matrix[i] = new double[col];
    }

    auto it = values.begin();
    for (int i = 0; i < fil; i++) {
        for (int j = 0; j < col; j++) {
            matrix[i][j] = *it++;
        }
    }
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
    if (this == &matrix2) return *this; // Evitar auto-asignación

    // Liberar memoria actual
    for (int i = 0; i < fil; i++) {
        delete[] matrix[i];
    }
    delete[] matrix;

    // Copiar datos de other
    fil = matrix2.fil;
    col = matrix2.col;
    matrix = new double*[fil];
    for (int i = 0; i < fil; i++) {
        matrix[i] = new double[col];
        for (int j = 0; j < col; j++) {
            matrix[i][j] = matrix2.matrix[i][j];
        }
    }
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
Matrix Matrix::operator*(const Matrix& matrix2) const
{
    Matrix result(fil, matrix2.col);

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
double& Matrix::operator()(int i, int j) const {
    if (i < 1 || i > fil || j < 1 || j > col) {
        throw std::out_of_range("Índice fuera de rango (1-based)");
    }
    return matrix[i - 1][j - 1];
}























// Suma por escalar (modifica la matriz actual)
Matrix& Matrix::operator+=(double scalar) {
    for (int i = 0; i < fil; i++) {
        for (int j = 0; j < col; j++) {
            matrix[i][j] += scalar;
        }
    }
    return *this;
}

Matrix& Matrix::operator-=(double scalar) {
    for (int i = 0; i < fil; ++i) {
        for (int j = 0; j < col; ++j) {
            matrix[i][j] -= scalar;
        }
    }
    return *this;
}


// Versiones no-miembro para suma
Matrix operator+(const Matrix& m, double scalar) {
    Matrix result(m.fil, m.col);
    for (int i = 0; i < m.fil; i++) {
        for (int j = 0; j < m.col; j++) {
            result.matrix[i][j] = m.matrix[i][j] + scalar;
        }
    }
    return result;
}

Matrix operator+(double scalar, const Matrix& m) {
    return m + scalar; // Conmutativa
}

Matrix operator-(const Matrix& m, double scalar) {
    Matrix result(m);
    return result -= scalar;
}


Matrix operator-(double scalar, const Matrix& m) {
    Matrix result(m.fil, m.col);
    for (int i = 0; i < m.fil; ++i) {
        for (int j = 0; j < m.col; ++j) {
            result.matrix[i][j] = scalar - m.matrix[i][j];
        }
    }
    return result;
}

Matrix& Matrix::operator*=(double scalar) {
    for (int i = 0; i < fil; i++) {
        for (int j = 0; j < col; j++) {
            matrix[i][j] *= scalar;
        }
    }
    return *this;
}

// Versión no-miembro (scalar * matrix)
Matrix operator*(double scalar, const Matrix& m) {
    Matrix result(m.fil, m.col);
    for (int i = 0; i < m.fil; i++) {
        for (int j = 0; j < m.col; j++) {
            result.matrix[i][j] = scalar * m.matrix[i][j];
        }
    }
    return result;
}

// Versión no-miembro (matrix * scalar)
Matrix operator*(const Matrix& m, double scalar) {
    return scalar * m; // Reutiliza la implementación anterior
}

/**
 * @brief Sobrecarga del operador de división por escalar (modifica la matriz actual)
 * @param scalar Escalar por el que dividir
 * @return Referencia a la matriz resultante
 */
Matrix& Matrix::operator/=(double scalar) {
    if (scalar == 0.0) {
        throw std::invalid_argument("División por cero no permitida");
    }

    for (int i = 0; i < fil; i++) {
        for (int j = 0; j < col; j++) {
            matrix[i][j] /= scalar;
        }
    }
    return *this;
}

/**
 * @brief Sobrecarga del operador de división matriz/escalar
 * @param m Matriz a dividir
 * @param scalar Escalar divisor
 * @return Nueva matriz resultante
 */
Matrix operator/(const Matrix& m, double scalar) {
    if (scalar == 0.0) {
        throw std::invalid_argument("División por cero no permitida");
    }

    Matrix result(m.fil, m.col);
    for (int i = 0; i < m.fil; i++) {
        for (int j = 0; j < m.col; j++) {
            result.matrix[i][j] = m.matrix[i][j] / scalar;
        }
    }
    return result;
}

/**
 * @brief Sobrecarga del operador de división escalar/matriz
 * @param scalar Escalar a dividir
 * @param m Matriz divisor
 * @return Nueva matriz resultante
 */
Matrix operator/(double scalar, const Matrix& m) {
    Matrix result(m.fil, m.col);
    for (int i = 0; i < m.fil; i++) {
        for (int j = 0; j < m.col; j++) {
            if (m.matrix[i][j] == 0.0) {
                throw std::invalid_argument("División por cero en elemento de matriz");
            }
            result.matrix[i][j] = scalar / m.matrix[i][j];
        }
    }
    return result;
}













/**
 * @brief Imprime la matriz en la consola.
 */
void Matrix::print() const
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
    if (fil != col) {
        throw std::runtime_error("La matriz debe ser cuadrada para calcular su inversa");
    }

    const int n = fil;
    Matrix result(n, n);
    Matrix temp(n, 2*n); // Matriz aumentada

    // Construir matriz aumentada [A|I]
    for (int i = 0; i < n; i++) {
        // Copiar la matriz original
        for (int j = 0; j < n; j++) {
            temp(i+1, j+1) = (*this)(i+1, j+1);
        }
        // Añadir matriz identidad
        for (int j = n; j < 2*n; j++) {
            temp(i+1, j+1) = (j - n == i) ? 1.0 : 0.0;
        }
    }

    // Eliminación hacia adelante
    for (int i = 1; i <= n; i++) {
        // Pivoteo parcial
        int maxRow = i;
        for (int k = i+1; k <= n; k++) {
            if (fabs(temp(k, i)) > fabs(temp(maxRow, i))) {
                maxRow = k;
            }
        }

        // Intercambiar filas si es necesario
        if (maxRow != i) {
            for (int j = 1; j <= 2*n; j++) {
                std::swap(temp(i, j), temp(maxRow, j));
            }
        }

        // Verificar si la matriz es singular
        if (fabs(temp(i, i)) < 1e-12) {
            throw std::runtime_error("Matriz singular, no se puede calcular la inversa");
        }

        // Hacer ceros debajo del pivote
        for (int k = i+1; k <= n; k++) {
            double factor = temp(k, i) / temp(i, i);
            for (int j = i; j <= 2*n; j++) {
                temp(k, j) -= factor * temp(i, j);
            }
        }
    }

    // Eliminación hacia atrás
    for (int i = n; i >= 1; i--) {
        // Normalizar la fila del pivote
        double pivot = temp(i, i);
        for (int j = i; j <= 2*n; j++) {
            temp(i, j) /= pivot;
        }

        // Hacer ceros arriba del pivote
        for (int k = i-1; k >= 1; k--) {
            double factor = temp(k, i);
            for (int j = i; j <= 2*n; j++) {
                temp(k, j) -= factor * temp(i, j);
            }
        }
    }

    // Extraer la matriz inversa (parte derecha de la matriz aumentada)
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            result(i, j) = temp(i, j + n);
        }
    }

    return result;
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

/**
 * @brief Calcula el producto punto entre dos vectores
 */
double Matrix::dot(const Matrix& a, const Matrix& b) {
    // Verificar que ambos son vectores (columna o fila)
    bool a_is_vector = (a.fil == 1 || a.col == 1);
    bool b_is_vector = (b.fil == 1 || b.col == 1);

    if (!a_is_vector || !b_is_vector) {
        throw std::invalid_argument("dot() solo aplicable a vectores");
    }

    // Obtener la longitud (usar la mayor dimensión)
    int n = (a.fil > a.col) ? a.fil : a.col;

    // Verificar que tienen la misma longitud
    if ((b.fil != n && b.col != n) || (b.fil != 1 && b.col != 1)) {
        throw std::invalid_argument("Vectores de dimensiones incompatibles para dot product");
    }

    double result = 0.0;

    // Producto punto para vectores columna
    if (a.col == 1 && b.col == 1) {
        for (int i = 0; i < n; ++i) {
            result += a.matrix[i][0] * b.matrix[i][0];
        }
    }
        // Producto punto para vectores fila
    else if (a.fil == 1 && b.fil == 1) {
        for (int j = 0; j < n; ++j) {
            result += a.matrix[0][j] * b.matrix[0][j];
        }
    }
        // Producto punto mixto (fila x columna)
    else {
        for (int i = 0; i < n; ++i) {
            double a_val = (a.col == 1) ? a.matrix[i][0] : a.matrix[0][i];
            double b_val = (b.col == 1) ? b.matrix[i][0] : b.matrix[0][i];
            result += a_val * b_val;
        }
    }

    return result;
}

/**
 * @brief Calcula el producto cruz entre dos vectores 3D
 */
Matrix Matrix::cross(const Matrix& a, const Matrix& b) {
    // Verificar que ambos son vectores 3D
    bool a_is_3d = (a.fil == 3 && a.col == 1) || (a.fil == 1 && a.col == 3);
    bool b_is_3d = (b.fil == 3 && b.col == 1) || (b.fil == 1 && b.col == 3);

    if (!a_is_3d || !b_is_3d) {
        throw std::invalid_argument("cross() solo aplicable a vectores 3D");
    }

    // Obtener componentes de los vectores
    double a1, a2, a3, b1, b2, b3;

    if (a.col == 1) { // Vector columna
        a1 = a.matrix[0][0];
        a2 = a.matrix[1][0];
        a3 = a.matrix[2][0];
    } else { // Vector fila
        a1 = a.matrix[0][0];
        a2 = a.matrix[0][1];
        a3 = a.matrix[0][2];
    }

    if (b.col == 1) { // Vector columna
        b1 = b.matrix[0][0];
        b2 = b.matrix[1][0];
        b3 = b.matrix[2][0];
    } else { // Vector fila
        b1 = b.matrix[0][0];
        b2 = b.matrix[0][1];
        b3 = b.matrix[0][2];
    }

    // Calcular producto cruz (siempre devuelve vector columna 3x1)
    Matrix result(3, 1);
    result.matrix[0][0] = a2 * b3 - a3 * b2;
    result.matrix[1][0] = a3 * b1 - a1 * b3;
    result.matrix[2][0] = a1 * b2 - a2 * b1;

    return result;
}


Matrix Matrix::getSubMatrix(int startRow, int endRow, int startCol, int endCol) const {
    // Validar índices
    if (startRow < 1 || endRow > fil || startCol < 1 || endCol > col ||
        startRow > endRow || startCol > endCol) {
        throw std::out_of_range("Índices de submatriz inválidos");
    }

    int newFil = endRow - startRow + 1;
    int newCol = endCol - startCol + 1;
    Matrix sub(newFil, newCol);

    for (int i = 1; i <= newFil; i++) {
        for (int j = 1; j <= newCol; j++) {
            sub(i, j) = (*this)(startRow + i - 1, startCol + j - 1);
        }
    }

    return sub;
}


Matrix Matrix::getColumn(int col) const {
    if (col < 0 || col > this->getColumnas()) {
        throw std::out_of_range("Índice de columna inválido");
    }
    Matrix result(this->getFilas(), 1);
    for (int i = 1; i < this->getFilas(); i++) {
        result(i, 1) = this->matrix[i-1][col-1];
    }
    return result;
}

Matrix Matrix::concatenate(const Matrix& other, int axis) const {
    if (axis == 0) { // Concatenación vertical
        if (col != other.col) {
            throw std::invalid_argument("El número de columnas debe coincidir para concatenar verticalmente");
        }

        Matrix result(fil + other.fil, col);
        for (int i = 0; i < fil; ++i) {
            for (int j = 0; j < col; ++j) {
                result.matrix[i][j] = matrix[i][j];
            }
        }
        for (int i = 0; i < other.fil; ++i) {
            for (int j = 0; j < other.col; ++j) {
                result.matrix[fil + i][j] = other.matrix[i][j];
            }
        }
        return result;
    }
    else if (axis == 1) { // Concatenación horizontal
        if (fil != other.fil) {
            throw std::invalid_argument("El número de filas debe coincidir para concatenar horizontalmente");
        }

        Matrix result(fil, col + other.col);
        for (int i = 0; i < fil; ++i) {
            for (int j = 0; j < col; ++j) {
                result.matrix[i][j] = matrix[i][j];
            }
            for (int j = 0; j < other.col; ++j) {
                result.matrix[i][col + j] = other.matrix[i][j];
            }
        }
        return result;
    }
    else {
        throw std::invalid_argument("El eje debe ser 0 (vertical) o 1 (horizontal)");
    }
}

void Matrix::setColumn(int col, const Matrix& column) {
    if (col < 1 || col > this->col) {
        throw std::invalid_argument("Column index out of bounds");
    }
    if (column.getFilas() != fil || column.getColumnas() != 1) {
        throw std::invalid_argument("Column vector must match matrix rows and be single column");
    }

    for (int i = 1; i <= fil; i++) {
        (*this)(i, col) = column(i, 1);
    }
}
