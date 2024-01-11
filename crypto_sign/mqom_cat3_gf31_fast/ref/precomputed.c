#include "precomputed.h"
#include "parameters.h"

const uint8_t interpolation_matrix_default[36] = {
    1, 22, 29, 29, 4, 8, 0, 5, 22, 12, 2, 22, 0, 26, 27, 8, 14, 18, 0, 24, 9, 17, 30, 13, 0, 22, 9, 28, 25, 9, 0, 25, 28, 30, 18, 23
};
const uint8_t interpolation_matrix_last[25] = {
    1, 16, 26, 28, 22, 0, 4, 6, 17, 5, 0, 28, 28, 29, 8, 0, 22, 8, 27, 5, 0, 23, 25, 23, 22
};

const uint8_t vanishing_polynomial[6] = {
    0, 4, 26, 23, 23, 16
};

const uint8_t inv_n1 = 26;
