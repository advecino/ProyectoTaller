//
// Created by advecino on 20/03/2025.
//
//
#include "Matrix.h"

void Matrix::initMatrix(){
    matrix = new double*[fil];
    for(int i=0;i<fil;i++){
        matrix[i] = new double[col];
    }
}