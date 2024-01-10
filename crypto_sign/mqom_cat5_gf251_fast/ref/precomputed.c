#include "precomputed.h"
#include "parameters.h"

const uint8_t interpolation_matrix_default[36] = {
    1, 27, 96, 198, 157, 23, 0, 5, 140, 118, 104, 136, 0, 246, 239, 16, 22, 230, 0, 87, 119, 25, 250, 21, 0, 187, 13, 197, 241, 115, 0, 201, 146, 199, 230, 228
};
const uint8_t interpolation_matrix_last[9] = {
    1, 124, 126, 0, 2, 250, 0, 125, 126
};

const uint8_t vanishing_polynomial[6] = {
    0, 131, 23, 26, 85, 236
};

const uint8_t inv_n1 = 42;
