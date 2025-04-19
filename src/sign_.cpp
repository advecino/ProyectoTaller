//
// Created by adria on 20/04/2025.
//

#include "../include/sign_.h"


#include <cmath>

double sign_(double a, double b) {
    if (b >= 0.0) {
        return std::abs(a);
    } else {
        return -std::abs(a);
    }
}
