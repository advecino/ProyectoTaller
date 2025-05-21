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

    /**
    * @brief Constructor que inicializa una matriz columna a partir de un vector de valores.
    * @param values Vector con los valores de la matriz.
    */
    explicit Matrix(const std::vector<double>& values);

    /**
     * @brief Constructor que inicializa una matriz con una lista de inicialización.
     * @param values Lista de valores para inicializar la matriz.
     * @param cols Número de columnas (por defecto 1).
     */
    Matrix(const std::initializer_list<double>& values, int cols = 1);


    /**
    * @brief Constructor que inicializa la matriz con un array de datos.
    * @param fil Número de filas.
    * @param col Número de columnas.
    * @param v Array con los datos a cargar.
    */
    Matrix(int fil, int col, double v[]);


    /**
     * @brief Devuelve una fila específica de la matriz.
     * @param fila Índice de la fila a obtener.
     * @return Matriz fila correspondiente.
     */
    Matrix getFila(int fila) const;
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
    void setColumn(int col, const Matrix& column);



    /**
     * @brief Conactenar 2 matrices verticalmente
     * @param m1 matriz1
     * @param m2 matriz2
     * @return matrizconcatenada
     */
    static Matrix concatenar( Matrix& m1,  Matrix& m2);



    /**
  * @brief Destructor de la matriz.
  */
    ~Matrix();

    /** @name Operadores */
    ///@{

    Matrix& operator= (const Matrix& matrix2);
    Matrix operator+ (const Matrix& matrix2);
    Matrix operator- (const Matrix& matrix2);
    Matrix operator* (const Matrix& matrix2) const;
    Matrix& operator/=(double scalar);
    double& operator()(int i, int j) const;

    Matrix& operator+=(double scalar);
    Matrix& operator-=(double scalar);
    Matrix& operator*=(double scalar);

    friend Matrix operator+(const Matrix& m, double scalar);
    friend Matrix operator+(double scalar, const Matrix& m);
    friend Matrix operator-(const Matrix& m, double scalar);
    friend Matrix operator-(double scalar, const Matrix& m);
    friend Matrix operator*(const Matrix& m, double scalar);
    friend Matrix operator*(double scalar, const Matrix& m);
    friend Matrix operator/(const Matrix& m, double scalar);
    friend Matrix operator/(double scalar, const Matrix& m);

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
    static double dot(const Matrix& a, const Matrix& b);

    /**
     * @brief Calcula el producto cruz entre dos vectores tridimensionales.
     * @param a Primer vector (3x1 o 1x3).
     * @param b Segundo vector (3x1 o 1x3).
     * @return Vector resultante del producto cruz.
     * @throw std::invalid_argument Si los vectores no son tridimensionales.
     */
    static Matrix cross(const Matrix& a, const Matrix& b);

    /**
     * @brief Obtiene una submatriz.
     * @param startRow Fila inicial (1-based).
     * @param endRow Fila final (1-based).
     * @param startCol Columna inicial (1-based).
     * @param endCol Columna final (1-based).
     * @return Submatriz extraída.
     */
    Matrix getSubMatrix(int startRow, int endRow, int startCol, int endCol) const;

private:
    /**
     * @brief Inicializa la matriz con ceros.
     */
    void initMatrix();

    int fil; ///< Número de filas.
    int col; ///< Número de columnas.
    double** matrix; ///< Puntero a la memoria de los elementos de la matriz.
};


#endif //PROYECTOTALLER_MATRIX_H
