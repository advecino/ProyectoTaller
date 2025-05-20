#include "../include/sign_.h"


#include <cmath>

double sign_(double a, double b) {
    if (b >= 0.0) {
        return std::abs(a);
    } else {
        return -std::abs(a);
    }
}
