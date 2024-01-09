#include <stdio.h>
#include <stdlib.h>

#include "witness.h"
#include "field.h"

void uncompress_instance(instance_t* inst) {
    if(inst->A == NULL) {
        // We assume here that
        //   inst->A == NULL iff inst->b == NULL
        inst->A = (void*) malloc(PARAM_m*PARAM_MATRIX_BYTESIZE);
        inst->b = (void*) malloc(PARAM_m*PARAM_n);
        prg_context entropy_ctx;
        samplable_t entropy = prg_to_samplable(&entropy_ctx);
        prg_init(&entropy_ctx, inst->seed, NULL);
        memset(*inst->A, 0, PARAM_m*PARAM_MATRIX_BYTESIZE);
        for(unsigned int i=0; i<PARAM_m; i++)
            random_points((*inst->A)[i], (PARAM_n*(PARAM_n+1))/2, &entropy);
        random_points(*inst->b, PARAM_m*PARAM_n, &entropy);
    }
}

static void compute_the_mq_equations_outputs(
        uint8_t y[PARAM_m], uint8_t x[PARAM_n],
        uint8_t A[PARAM_m][PARAM_MATRIX_BYTESIZE], uint8_t b[PARAM_m][PARAM_n]) {

    uint8_t tmp[PARAM_n];
    for(int i=0; i<PARAM_m; i++) {
        // tmp = b_i + A_i^T x
        memcpy(tmp, b[i], PARAM_n);
        // We use "matrows_muladd" instead of "matcols_muladd"
        //    since the matrix A_i is transposed
        int ind = 0;
        for(uint32_t j=0; j<PARAM_n; j++) {
            tmp[j] = add_points(tmp[j], innerproduct_points(x+j, A[i]+ind, PARAM_n-j));
            ind += (PARAM_n-j);
        }
        // y_i = tmp^T x
        y[i] = 0;
        matcols_muladd(&y[i], x, tmp, PARAM_n, 1);
    }
}

void generate_instance_with_solution(instance_t** inst, solution_t** sol, samplable_t* entropy) {
    // Allocate
    *sol = (solution_t*) malloc(sizeof(solution_t));
    *inst = (instance_t*) malloc(sizeof(instance_t));
    (*inst)->A = NULL;
    (*inst)->b = NULL;

    // Extended Witness
    random_points((*sol)->x, PARAM_n, entropy);

    // Sample a seed
    byte_sample(entropy, (*inst)->seed, PARAM_SEED_SIZE);

    // Build random matrices
    uncompress_instance(*inst);

    // Build y
    compute_the_mq_equations_outputs(
        (*inst)->y, (*sol)->x,
        (*(*inst)->A), (*(*inst)->b)
    );
}

void free_instance_solution(solution_t* sol) {
    free(sol);
}

void free_instance(instance_t* inst) {
    if(inst->A) {
        free(inst->A);
        inst->A = NULL;
    }
    if(inst->b) {
        free(inst->b);
        inst->b = NULL;
    }
    free(inst);
}

int are_same_instances(instance_t* inst1, instance_t* inst2) {
    int ret = 0;
    ret |= memcmp(inst1->seed, inst2->seed, PARAM_SEED_SIZE);
    ret |= memcmp(inst1->y, inst2->y, PARAM_m);
    return ((ret) ? 0 : 1);
}

int is_correct_solution(instance_t* inst, solution_t* sol) {
    int ret = 0;
    uncompress_instance(inst);

    uint8_t y_candidate[PARAM_m];
    compute_the_mq_equations_outputs(
        y_candidate, sol->x,
        (*inst->A), (*inst->b)
    );
    ret = memcmp(inst->y, y_candidate, PARAM_m);

    return (ret == 0);
}
