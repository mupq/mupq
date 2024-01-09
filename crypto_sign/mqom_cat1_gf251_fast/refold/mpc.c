#include "mpc.h"
#include "precomputed.h"

int is_valid_plain_broadcast(mpc_broadcast_t const* plain_br) {
    uint8_t ret = 0;
    // Just test if "plain_br->v" is zero
    for(unsigned int i=0; i<sizeof(plain_br->v); i++)
        ret |= ((uint8_t*) plain_br->v)[i];
    return (ret == 0);
}

void compute_hint(mpc_hint_t* hint, const mpc_wit_t* wit, const mpc_unif_t* unif, const instance_t* inst, mpc_challenge_1_t* mpc_challenge_1) {
    uint8_t tmp[PARAM_eta];

    // Compute w = linA x
    uint8_t w[PARAM_n][PARAM_eta] = {0};
    #ifdef PARAM_PRECOMPUTE_LIN_A
    (void) inst;
    matcols_muladd_triangular(w, wit->x, mpc_challenge_1->linA, PARAM_n, PARAM_eta);
    #else
    uint8_t intermediary[PARAM_n];
    for(unsigned int j=0; j<PARAM_m; j++) {
        memset(intermediary, 0, sizeof(intermediary));
        matcols_muladd_triangular(intermediary, wit->x, (*inst->A)[j], PARAM_n, 1);
        for(unsigned int i=0; i<PARAM_n; i++) {
            mul_points_mixed(tmp, intermediary[i], mpc_challenge_1->gamma[j]);
            add_points_ext(w[i], w[i], tmp);
        }
    }
    #endif

    // Interpolation of w -> \hat{W}
    uint8_t w_hat[PARAM_n][PARAM_eta] = {0};
    for(unsigned int i=0; i<PARAM_eta; i++) {
        uint8_t input[PARAM_n] = {0};
        uint8_t partial[PARAM_n] = {0};
        for(unsigned int j=0; j<PARAM_n; j++)
            input[j] = w[j][i];
        for(unsigned int offset=0; offset<(PARAM_n2-1)*PARAM_n1; offset+=PARAM_n1)
            matcols_muladd(partial+offset, input+offset, interpolation_matrix_default, PARAM_n1, PARAM_n1);
        matcols_muladd(partial+(PARAM_n-PARAM_last), input+(PARAM_n-PARAM_last), interpolation_matrix_last, PARAM_last, PARAM_last);        
        for(unsigned int j=0; j<PARAM_n; j++)
            w_hat[j][i] = partial[j];
    }

    // Interpolation of x
    uint8_t x_extended[PARAM_n2][PARAM_n1] = {0};
    memcpy((uint8_t*) x_extended, wit->x, PARAM_n);
    uint8_t x_poly[PARAM_n2][PARAM_n1] = {0};
    for(unsigned int i=0; i<PARAM_n2; i++)
        matcols_muladd(x_poly[i], x_extended[i], interpolation_matrix_default, PARAM_n1, PARAM_n1);

    // Compute Q
    memset(hint->q_poly, 0, sizeof(hint->q_poly));
    for(unsigned int j=0; j<PARAM_n2; j++) {
        unsigned int offset = j*PARAM_n1;
        unsigned int size = (j == PARAM_n2-1) ? PARAM_last :PARAM_n1;

        uint8_t w_tilde[PARAM_n1+1][PARAM_eta] = {0};            
        // we assume here that "vanishing(0)=0"
        for(unsigned int i=1; i<PARAM_n1; i++)
            mul_points_mixed(w_tilde[i], vanishing_polynomial[i], unif->a[j]);
        memcpy(w_tilde[PARAM_n1],unif->a[j],PARAM_eta);
        // \tilde{w}(X) += \hat{w}(X)
        add_tab_points(w_tilde, w_hat+offset, size*PARAM_eta);

        for(unsigned int i=0; i<PARAM_n1+1; i++)
            for(unsigned int k=0; k<PARAM_n1; k++) {
                if(i+k == 0) continue;
                mul_points_mixed(tmp, x_poly[j][k], w_tilde[i]);
                add_points_ext(hint->q_poly[i+k-1], hint->q_poly[i+k-1], tmp);
            }
    }
}

static void run_multiparty_computation(mpc_broadcast_t* broadcast,
                                mpc_challenge_1_t* mpc_challenge_1, mpc_challenge_2_t* mpc_challenge_2,
                                mpc_share_t const* share, mpc_broadcast_t const* plain_br,
                                const instance_t* inst, char has_sharing_offset, char entire_computation) {

    




    uint8_t tmp[PARAM_eta];
    uint8_t acc[PARAM_eta];

    // [[w]] = linA [[x]]
    uint8_t w[PARAM_n][PARAM_eta] = {0};
    #ifdef PARAM_PRECOMPUTE_LIN_A
    matcols_muladd_triangular(w, share->wit.x, mpc_challenge_1->linA, PARAM_n, PARAM_eta);
    #else
    uint8_t intermediary[PARAM_n];
    for(unsigned int j=0; j<PARAM_m; j++) {
        memset(intermediary, 0, sizeof(intermediary));
        matcols_muladd_triangular(intermediary, share->wit.x, (*inst->A)[j], PARAM_n, 1);
        for(unsigned int i=0; i<PARAM_n; i++) {
            mul_points_mixed(tmp, intermediary[i], mpc_challenge_1->gamma[j]);
            add_points_ext(w_hat[i], w_hat[i], tmp);
        }
    }
    #endif
    // Interpolation of [[w]]
    uint8_t w_hat[PARAM_n][PARAM_eta] = {0};
    for(unsigned int i=0; i<PARAM_eta; i++) {
        uint8_t input[PARAM_n] = {0};
        uint8_t partial[PARAM_n] = {0};
        for(unsigned int j=0; j<PARAM_n; j++)
            input[j] = w[j][i];
        for(unsigned int offset=0; offset<(PARAM_n2-1)*PARAM_n1; offset+=PARAM_n1)
            matcols_muladd(partial+offset, input+offset, interpolation_matrix_default, PARAM_n1, PARAM_n1);
        matcols_muladd(partial+(PARAM_n-PARAM_last), input+(PARAM_n-PARAM_last), interpolation_matrix_last, PARAM_last, PARAM_last);        
        for(unsigned int j=0; j<PARAM_n; j++)
            w_hat[j][i] = partial[j];
    }

    // Parse [[W^]]
    // [[W~]] = [[W^]] + [[a]]*vanishing
    uint8_t w_poly[PARAM_n2][PARAM_n1+1][PARAM_eta] = {0};   
    for(unsigned int j=0; j<PARAM_n2; j++) {
        unsigned int offset = j*PARAM_n1;
        unsigned int size = (j == PARAM_n2-1) ? PARAM_last : PARAM_n1;

        // we assume here that "vanishing(0)=0"
        for(unsigned int i=1; i<PARAM_n1; i++)
            mul_points_mixed(w_poly[j][i], vanishing_polynomial[i], share->unif.a[j]);
        memcpy(w_poly[j][PARAM_n1],share->unif.a[j],PARAM_eta);
        // \tilde{w}(X) += \hat{w}(X)
        add_tab_points(w_poly[j], w_hat+offset, size*PARAM_eta);
    }

    // [[\alpha]] = [[W~]](r)
    memset(broadcast->alpha, 0, sizeof(broadcast->alpha));
    for(unsigned int j=0; j<PARAM_n2; j++) {
        for(unsigned int i=0; i<PARAM_n1+1; i++) {
            mul_points_ext(tmp, mpc_challenge_2->r_powers[i], w_poly[j][i]);
            add_points_ext(broadcast->alpha[j], broadcast->alpha[j], tmp);
        }
    }

    if(entire_computation) {
        // [[z]] = sum_{j=1}^m gamma_j * (y_j - b_j^T [[x]])
        uint8_t z[PARAM_eta] = {0};
        for(int j=0; j<PARAM_m; j++) {
            uint8_t res = innerproduct_points((*inst->b)[j], share->wit.x, PARAM_n);
            if(has_sharing_offset)
                res = sub_points(inst->y[j], res);
            else
                res = neg_point(res);
            mul_points_mixed(tmp, res, mpc_challenge_1->gamma[j]);
            add_points_ext(z, z, tmp);
        }

        // Interpolation of [[x]]
        uint8_t x_extended[PARAM_n2][PARAM_n1] = {0};
        memcpy((uint8_t*) x_extended, share->wit.x, PARAM_n);
        uint8_t x_poly[PARAM_n2][PARAM_n1] = {0};
        for(unsigned int i=0; i<PARAM_n2; i++)
            matcols_muladd(x_poly[i], x_extended[i], interpolation_matrix_default, PARAM_n1, PARAM_n1);

        // [[beta]] = (1/n1) * ([[z]] - sum_i i*[[Q]](i))
        uint8_t beta[PARAM_eta];
        memset(beta, 0, sizeof(beta));
        for(uint8_t i=0; i<PARAM_n1; i++) {
            uint8_t i_power = i;
            for(unsigned int k=0; k<2*PARAM_n1-1; k++) {
                mul_points_mixed(tmp, i_power, share->hint.q_poly[k]);
                add_points_ext(beta, beta, tmp);
                i_power = mul_points(i_power, i);
            }
        }
        sub_points_ext(beta, z, beta);
        mul_points_mixed(beta, inv_n1, beta);

        // [[v]] = r*[[Q]](r) + beta - sum_j \alpha [[X^]](r)
        memset(broadcast->v, 0, sizeof(broadcast->v));
        for(unsigned int j=0; j<PARAM_n2; j++) {
            memset(acc, 0, sizeof(tmp));
            for(unsigned int i=0; i<PARAM_n1; i++) {
                mul_points_mixed(tmp, x_poly[j][i], mpc_challenge_2->r_powers[i]);
                add_points_ext(acc, acc, tmp);
            }
            mul_points_ext(acc, acc, plain_br->alpha[j]);
            add_points_ext(broadcast->v, broadcast->v, acc);
        }
        neg_tab_points(broadcast->v, PARAM_eta);
        add_points_ext(broadcast->v, broadcast->v, beta);
        for(unsigned int k=0; k<2*PARAM_n1-1; k++) {
            mul_points_ext(tmp, share->hint.q_poly[k], mpc_challenge_2->r_powers[k+1]);
            add_points_ext(broadcast->v, broadcast->v, tmp);
        }
    } else {
        memset(broadcast->v, 0, sizeof(broadcast->v));
    }

}

void mpc_compute_plain_broadcast(mpc_broadcast_t* plain_br,
                                mpc_challenge_1_t* mpc_challenge_1, mpc_challenge_2_t* mpc_challenge_2,
                                mpc_share_t const* share, const instance_t* inst) {
    run_multiparty_computation(plain_br, mpc_challenge_1, mpc_challenge_2, share, plain_br, inst, 1, 0);
}

void mpc_compute_communications(mpc_broadcast_t* broadcast,
                                mpc_challenge_1_t* mpc_challenge_1, mpc_challenge_2_t* mpc_challenge_2,
                                mpc_share_t const* share, mpc_broadcast_t const* plain_br,
                                const instance_t* inst, char has_sharing_offset) {
    run_multiparty_computation(broadcast, mpc_challenge_1, mpc_challenge_2, share, plain_br, inst, has_sharing_offset, 1);
}

void expand_mpc_challenge_hash_1(mpc_challenge_1_t** challenges, const uint8_t* digest, size_t nb, instance_t* inst) {

    xof_context entropy_ctx;
    xof_init(&entropy_ctx);
    xof_update(&entropy_ctx, digest, PARAM_DIGEST_SIZE);
    samplable_t entropy = xof_to_samplable(&entropy_ctx);

    for(size_t num=0; num<nb; num++) {
        vec_rnd((uint8_t*) challenges[num]->gamma, sizeof((*challenges)->gamma), &entropy);

        #ifdef PARAM_PRECOMPUTE_LIN_A
            uint8_t gamma_[PARAM_eta][PARAM_m];
            uint8_t linA_[PARAM_eta][PARAM_MATRIX_BYTESIZE] = {0};
            for(unsigned int i=0; i<PARAM_eta; i++) {
                for(unsigned int j=0; j<PARAM_m; j++)
                    gamma_[i][j] = challenges[num]->gamma[j][i];

                mat128cols_muladd(linA_[i], gamma_[i], (*inst->A), PARAM_m, PARAM_MATRIX_BYTESIZE);

                for(unsigned int j=0; j<PARAM_MATRIX_BYTESIZE; j++)
                        challenges[num]->linA[j][i] = linA_[i][j];
            }
        #else
        (void) inst;
        #endif
    }
}

void expand_mpc_challenge_hash_2(mpc_challenge_2_t** challenges, const uint8_t* digest, size_t nb, instance_t* inst) {

    xof_context entropy_ctx;
    xof_init(&entropy_ctx);
    xof_update(&entropy_ctx, digest, PARAM_DIGEST_SIZE);
    xof_final(&entropy_ctx);
    samplable_t entropy = xof_to_samplable(&entropy_ctx);

    (void) inst;
    for(size_t num=0; num<nb; num++) {
        unsigned int valid_challenge = 0;
        while(!valid_challenge) {
            vec_rnd((uint8_t*) challenges[num]->r, sizeof((*challenges)->r), &entropy);

            uint8_t val = 0;
            for(unsigned int i=1; i<PARAM_eta; i++)
                val |= challenges[num]->r[i];
            valid_challenge = !((val == 0) && (challenges[num]->r[0] < PARAM_n1));
        }
    }

    for(size_t num=0; num<nb; num++) {
        uint8_t* point = challenges[num]->r;
        // powers[0] = 1
        memset(&challenges[num]->r_powers[0], 0, PARAM_eta);
        challenges[num]->r_powers[0][0] = 1;
        // powers[1] = point
        memcpy(&challenges[num]->r_powers[1], point, PARAM_eta);
        // powers[i] = evals[i-1]*point
        for(unsigned int i=2; i<2*PARAM_n1; i++)
            mul_points_ext(
                challenges[num]->r_powers[i],
                challenges[num]->r_powers[i-1],
                point
            );
    }
}
