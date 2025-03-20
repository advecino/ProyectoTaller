//
// Created by advecino on 20/03/2025.
//
//
#include"../include/Matrix.h"
#include "../include/Matrix.h"

void Matrix::initMatrix(){
    matrix = new double*[fil];
    for(int i=0;i<fil;i++){
        matrix[i] = new double[col];
    }
    for (int i=0; i<fil;i++){
        for(int j = 0; j<col;j++){
            matrix[i][j]= 0.0;
        }
    }
}