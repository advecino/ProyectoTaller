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
    Matrix(int fil, int col, double v[], int n);

    /**
    * @brief Constructor de copia.
    * @param m Objeto Matrix a copiar.
    */
    Matrix(const Matrix &m);

    /**
    * @brief Constructor que inicializa una matriz columna a partir de un vector de valores.
    * @param values Vector con los valores de la matriz.
    */
    explicit Matrix(const std::vector<double> &values);

    /**
     * @brief Constructor que inicializa una matriz con una lista de inicialización.
     * @param values Lista de valores para inicializar la matriz.
     * @param cols Número de columnas (por defecto 1).
     */
    Matrix(const std::initializer_list<double> &values, int cols = 1);


    /**
    * @brief Constructor que inicializa la matriz con un array de datos.
    * @param fil Número de filas.
    * @param col Número de columnas.
    * @param v Array con los datos a cargar.
    */
    Matrix(int fil, int col, double v[]);
    /**
    * @brief Constructor de la matriz con dimensiones especificadas.
    * @param fil Número de filas.
    * @param col Número de columnas.
    */
    static Matrix Matrixx(int fil, int col);

    /**
     * @brief Devuelve una columna específica de la matriz.
     * @param col Índice de la columna a obtener.
     * @return Matriz columna correspondiente.
     */
    Matrix getColumn(int col) const;


    /**
     * @brief Asigna una matriz columna a una columna específica.
     * @param col Índice de la columna a modificar.
     * @param column Matriz columna a asignar.
     */
    void setColumn(int col, const Matrix &column);


/**
 * @brief Concatena dos matrices a lo largo de un eje.
 * @param A Primera matriz.
 * @param B Segunda matriz.
 * @param axis 0 para concatenación vertical (filas), 1 para horizontal (columnas).
 * @return Nueva matriz resultado de la concatenación.
 * @throw std::invalid_argument Si las dimensiones no son compatibles.
 */
    static Matrix concatenate(Matrix &A, Matrix &B, int axis = 0);


    /**
  * @brief Destructor de la matriz.
  */
    ~Matrix();

    /** @name Operadores */
    ///@{

    Matrix &operator=(const Matrix &matrix2);

    Matrix operator+(const Matrix &matrix2);

    Matrix operator-(const Matrix &matrix2);

    Matrix operator*(const Matrix &matrix2) const;

    Matrix &operator/=(double scalar);

    double &operator()(int i, int j) const;

    Matrix &operator+=(double scalar);

    Matrix &operator-=(double scalar);

    Matrix &operator*=(double scalar);

    friend Matrix operator+(const Matrix &m, double scalar);

    friend Matrix operator+(double scalar, const Matrix &m);

    friend Matrix operator-(const Matrix &m, double scalar);

    friend Matrix operator-(double scalar, const Matrix &m);

    friend Matrix operator*(const Matrix &m, double scalar);

    friend Matrix operator*(double scalar, const Matrix &m);

    friend Matrix operator/(const Matrix &m, double scalar);

    friend Matrix operator/(double scalar, const Matrix &m);

    ///@}

    /**
     * @brief Imprime la matriz por consola.
     */
    void print() const;

    /**
     * @brief Calcula la norma euclídea de la matriz.
     * @return Valor de la norma.
     */
    double norm() const;

    /**
     * @brief Devuelve el número de filas.
     * @return Número de filas de la matriz.
     */
    int getFilas() const;

    /**
     * @brief Devuelve el número de columnas.
     * @return Número de columnas de la matriz.
     */
    int getColumnas() const;

    /**
     * @brief Devuelve la transpuesta de la matriz.
     * @return Matriz transpuesta.
     */
    Matrix transpuesta() const;

    /**
     * @brief Devuelve la inversa de la matriz (si existe).
     * @return Matriz inversa.
     */
    Matrix inversa() const;

    /**
     * @brief Calcula el producto punto entre dos vectores.
     * @param a Primer vector (nx1 o 1xn).
     * @param b Segundo vector (nx1 o 1xn).
     * @return Resultado del producto escalar.
     * @throw std::invalid_argument Si las dimensiones no son compatibles.
     */
    static double dot(const Matrix &a, const Matrix &b);

    /**
     * @brief Calcula el producto cruz entre dos vectores tridimensionales.
     * @param a Primer vector (3x1 o 1x3).
     * @param b Segundo vector (3x1 o 1x3).
     * @return Vector resultante del producto cruz.
     * @throw std::invalid_argument Si los vectores no son tridimensionales.
     */
    static Matrix cross(const Matrix &a, const Matrix &b);

    /**
     * @brief Obtiene una submatriz.
     * @param startRow Fila inicial (1-based).
     * @param endRow Fila final (1-based).
     * @param startCol Columna inicial (1-based).
     * @param endCol Columna final (1-based).
     * @return Submatriz extraída.
     */
    Matrix getSubMatrix(int startRow, int endRow, int startCol, int endCol) const;

    /**
 * @brief Obtiene una submatriz de la matriz actual
 * @param startRow Fila inicial (basado en 1)
 * @param endRow Fila final (basado en 1)
 * @param startCol Columna inicial (basado en 1)
 * @param endCol Columna final (basado en 1)
 * @return Nueva matriz con la submatriz extraída
 * @throw std::out_of_range Si los índices están fuera de rango
 */
    Matrix getSubMatrix2(int startRow, int endRow, int startCol, int endCol) const;


private:
    /**
     * @brief Inicializa la matriz con ceros.
     */
    void initMatrix();

    int fil; ///< Número de filas.
    int col; ///< Número de columnas.
    double **matrix; ///< Puntero a la memoria de los elementos de la matriz.
};


#endif //PROYECTOTALLER_MATRIX_H
