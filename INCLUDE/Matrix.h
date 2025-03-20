//
// Created by advecino on 20/03/2025.
//

#ifndef PROYECTOTALLER_MATRIX_H
#define PROYECTOTALLER_MATRIX_H


class Matrix {
public:
    Matrix(int fil, int col);
    Matrix(int fil, int col, double v[],int n);
    Matrix(const Matrix& m);
    ~Matrix();

    Matrix& operator= (const Matrix& matrix2);
    Matrix operator+ (const Matrix& matrix2);
    Matrix operator- (const Matrix& matrix2);
    Matrix operator* (const Matrix& matrix2);
    double& operator()(const int i, const int j) const;

    void print();

private:
    void initMatrix();

};


#endif //PROYECTOTALLER_MATRIX_H
