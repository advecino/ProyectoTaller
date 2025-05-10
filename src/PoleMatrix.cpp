//
// Created by adria on 11/05/2025.
//

#include "../include/PoleMatrix.h"


Matrix PoleMatrix(double xp, double yp) {
    return R_y(-xp) * R_x(-yp);
}