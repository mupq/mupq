#include "precomputed.h"
#include "parameters.h"

const uint8_t interpolation_matrix_default[25] = {
    1, 16, 26, 28, 22, 0, 4, 6, 17, 5, 0, 28, 28, 29, 8, 0, 22, 8, 27, 5, 0, 23, 25, 23, 22
};
const uint8_t interpolation_matrix_last[16] = {
    1, 24, 1, 5, 0, 3, 13, 16, 0, 14, 2, 15, 0, 21, 15, 26
};

const uint8_t vanishing_polynomial[5] = {
    0, 24, 12, 4, 21
};

const uint8_t inv_n1 = 25;
