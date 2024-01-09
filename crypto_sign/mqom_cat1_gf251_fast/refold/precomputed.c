#include "precomputed.h"
#include "parameters.h"

const uint8_t interpolation_matrix_default[16] = {
    1, 40, 1, 209, 0, 3, 123, 126, 0, 124, 2, 125, 0, 84, 125, 42
};
const uint8_t interpolation_matrix_last[9] = {
    1, 124, 126, 0, 2, 250, 0, 125, 126
};

const uint8_t vanishing_polynomial[4] = {
    0, 245, 11, 245
};

const uint8_t inv_n1 = 63;
