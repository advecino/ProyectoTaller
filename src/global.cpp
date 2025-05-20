#include "../include/global.h"


double Cnm[300][300] = {0};
double Snm[300][300] = {0};

Matrix PC = [](){
        Matrix I(3,3);
        for(int i=1;i<=3;++i) I(i,i) = 1.0;
        return I;
    }();
Matrix eopdata(13, 3);