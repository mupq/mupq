#include "precomputed.h"
#include "parameters.h"

const uint8_t interpolation_matrix_default[25] = {
    1, 228, 242, 146, 136, 0, 4, 163, 127, 209, 0, 248, 193, 249, 63, 0, 85, 165, 43, 209, 0, 188, 241, 188, 136
};
const uint8_t interpolation_matrix_last[9] = {
    1, 124, 126, 0, 2, 250, 0, 125, 126
};

const uint8_t vanishing_polynomial[5] = {
    0, 24, 201, 35, 241
};

const uint8_t inv_n1 = 201;
