#include "../include/Matrix.h"
#include <iostream>
#include <iomanip>
#include <cmath>





Matrix::Matrix(int fil, int col) : fil(fil), col(col)
{
    initMatrix();
}


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


Matrix::Matrix(const Matrix& m) : fil(m.fil), col(m.col) {
    matrix = new double*[fil];
    for (int i = 0; i < fil; i++) {
        matrix[i] = new double[col];
        for (int j = 0; j < col; j++) {
            matrix[i][j] = m.matrix[i][j];
        }
    }
}



void Matrix::initMatrix()
{
    matrix = new double*[fil];
    for (int i = 0; i < fil; i++)
        matrix[i] = new double[col];

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            matrix[i][j] = 0.0;
}

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

Matrix Matrix::Matrixx(int fil,int col) {
    double fallback[6] = {
            5753.173e3 + 39.0,
            2673.361e3 + 40.9,
            3440.304e3 - 14.6,
            4326.31,
            -1926.7,
            -5726.12
    };
    Matrix Y(6, 1);
    for (int i = 0; i < 6; ++i)
        Y(i + 1, 1) = fallback[i];
    return Y;
}



Matrix::~Matrix()
{
    for (int i = 0; i < fil; i++)
        delete[] matrix[i];

    delete[] matrix;
}


Matrix& Matrix::operator=(const Matrix& matrix2)
{
    if (this == &matrix2) return *this;


    for (int i = 0; i < fil; i++) {
        delete[] matrix[i];
    }
    delete[] matrix;

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


Matrix Matrix::operator+(const Matrix& matrix2)
{
    Matrix result(fil, col);

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] + matrix2.matrix[i][j];

    return result;
}


Matrix Matrix::operator-(const Matrix& matrix2)
{
    Matrix result(fil, col);

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] - matrix2.matrix[i][j];

    return result;
}

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


double& Matrix::operator()(int i, int j) const {
    if (i < 1 || i > fil || j < 1 || j > col) {
        throw std::out_of_range("Índice fuera de rango (1-based)");
    }
    return matrix[i - 1][j - 1];
}


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
    return m + scalar;
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

Matrix operator*(double scalar, const Matrix& m) {
    Matrix result(m.fil, m.col);
    for (int i = 0; i < m.fil; i++) {
        for (int j = 0; j < m.col; j++) {
            result.matrix[i][j] = scalar * m.matrix[i][j];
        }
    }
    return result;
}

Matrix operator*(const Matrix& m, double scalar) {
    return scalar * m;
}


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


int Matrix::getFilas() const {
    return fil;
}


int Matrix::getColumnas() const {
    return col;
}


Matrix Matrix::transpuesta() const {
    Matrix result(fil, col);
    for(int i = 1; i <= fil; ++i) {
        for(int j = 1; j <= col; ++j) {
            result(j,i) = (*this)(i,j);
        }
    }
    return result;
}



Matrix Matrix::inversa() const {
    if (fil != col) {
        throw std::runtime_error("La matriz debe ser cuadrada para calcular su inversa");
    }

    const int n = fil;
    Matrix result(n, n);
    Matrix temp(n, 2*n);


    for (int i = 0; i < n; i++) {

        for (int j = 0; j < n; j++) {
            temp(i+1, j+1) = (*this)(i+1, j+1);
        }

        for (int j = n; j < 2*n; j++) {
            temp(i+1, j+1) = (j - n == i) ? 1.0 : 0.0;
        }
    }


    for (int i = 1; i <= n; i++) {

        int maxRow = i;
        for (int k = i+1; k <= n; k++) {
            if (fabs(temp(k, i)) > fabs(temp(maxRow, i))) {
                maxRow = k;
            }
        }


        if (maxRow != i) {
            for (int j = 1; j <= 2*n; j++) {
                std::swap(temp(i, j), temp(maxRow, j));
            }
        }


        if (fabs(temp(i, i)) < 1e-12) {
            throw std::runtime_error("Matriz singular, no se puede calcular la inversa");
        }

        for (int k = i+1; k <= n; k++) {
            double factor = temp(k, i) / temp(i, i);
            for (int j = i; j <= 2*n; j++) {
                temp(k, j) -= factor * temp(i, j);
            }
        }
    }


    for (int i = n; i >= 1; i--) {

        double pivot = temp(i, i);
        for (int j = i; j <= 2*n; j++) {
            temp(i, j) /= pivot;
        }

        for (int k = i-1; k >= 1; k--) {
            double factor = temp(k, i);
            for (int j = i; j <= 2*n; j++) {
                temp(k, j) -= factor * temp(i, j);
            }
        }
    }

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


double Matrix::dot(const Matrix& a, const Matrix& b) {

    bool a_is_vector = (a.fil == 1 || a.col == 1);
    bool b_is_vector = (b.fil == 1 || b.col == 1);

    if (!a_is_vector || !b_is_vector) {
        throw std::invalid_argument("dot() solo aplicable a vectores");
    }

    int n = (a.fil > a.col) ? a.fil : a.col;

    if ((b.fil != n && b.col != n) || (b.fil != 1 && b.col != 1)) {
        throw std::invalid_argument("Vectores de dimensiones incompatibles para dot product");
    }

    double result = 0.0;

    if (a.col == 1 && b.col == 1) {
        for (int i = 0; i < n; ++i) {
            result += a.matrix[i][0] * b.matrix[i][0];
        }
    }

    else if (a.fil == 1 && b.fil == 1) {
        for (int j = 0; j < n; ++j) {
            result += a.matrix[0][j] * b.matrix[0][j];
        }
    }

    else {
        for (int i = 0; i < n; ++i) {
            double a_val = (a.col == 1) ? a.matrix[i][0] : a.matrix[0][i];
            double b_val = (b.col == 1) ? b.matrix[i][0] : b.matrix[0][i];
            result += a_val * b_val;
        }
    }

    return result;
}


Matrix Matrix::cross(const Matrix& a, const Matrix& b) {

    bool a_is_3d = (a.fil == 3 && a.col == 1) || (a.fil == 1 && a.col == 3);
    bool b_is_3d = (b.fil == 3 && b.col == 1) || (b.fil == 1 && b.col == 3);

    if (!a_is_3d || !b_is_3d) {
        throw std::invalid_argument("cross() solo aplicable a vectores 3D");
    }


    double a1, a2, a3, b1, b2, b3;

    if (a.col == 1) {
        a1 = a.matrix[0][0];
        a2 = a.matrix[1][0];
        a3 = a.matrix[2][0];
    } else {
        a1 = a.matrix[0][0];
        a2 = a.matrix[0][1];
        a3 = a.matrix[0][2];
    }

    if (b.col == 1) {
        b1 = b.matrix[0][0];
        b2 = b.matrix[1][0];
        b3 = b.matrix[2][0];
    } else {
        b1 = b.matrix[0][0];
        b2 = b.matrix[0][1];
        b3 = b.matrix[0][2];
    }

    Matrix result(3, 1);
    result.matrix[0][0] = a2 * b3 - a3 * b2;
    result.matrix[1][0] = a3 * b1 - a1 * b3;
    result.matrix[2][0] = a1 * b2 - a2 * b1;

    return result;
}


Matrix Matrix::getSubMatrix(int startRow, int endRow, int startCol, int endCol) const {

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

/**
 * @brief Obtiene una submatriz de la matriz actual
 * @param startRow Fila inicial (basado en 1)
 * @param endRow Fila final (basado en 1)
 * @param startCol Columna inicial (basado en 1)
 * @param endCol Columna final (basado en 1)
 * @return Nueva matriz con la submatriz extraída
 * @throw std::out_of_range Si los índices están fuera de rango
 */
Matrix Matrix::getSubMatrix2(int startRow, int endRow, int startCol, int endCol) const {
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



Matrix Matrix::concatenate(Matrix& A, Matrix& B, int axis) {
    if (axis == 0) {
        if (A.col != B.col) {
            throw std::invalid_argument("Concatenación vertical requiere mismo número de columnas.");
        }
        Matrix result(A.fil + B.fil, A.col);
        for (int i = 1; i <= A.fil; ++i)
            for (int j = 1; j <= A.col; ++j)
                result(i, j) = A(i, j);
        for (int i = 1; i <= B.fil; ++i)
            for (int j = 1; j <= B.col; ++j)
                result(i + A.fil, j) = B(i, j);
        return result;
    } else if (axis == 1) {
        if (A.fil != B.fil) {
            throw std::invalid_argument("Concatenación horizontal requiere mismo número de filas.");
        }
        Matrix result(A.fil, A.col + B.col);
        for (int i = 1; i <= A.fil; ++i) {
            for (int j = 1; j <= A.col; ++j)
                result(i, j) = A(i, j);
            for (int j = 1; j <= B.col; ++j)
                result(i, j + A.col) = B(i, j);
        }
        return result;
    } else {
        throw std::invalid_argument("El parámetro 'axis' debe ser 0 (vertical) o 1 (horizontal).");
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
