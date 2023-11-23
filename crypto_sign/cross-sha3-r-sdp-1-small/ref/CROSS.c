/**
 *
 * Reference ISO-C11 Implementation of CROSS.
 *
 * @version 1.1 (March 2023)
 *
 * @author Alessandro Barenghi <alessandro.barenghi@polimi.it>
 * @author Gerardo Pelosi <gerardo.pelosi@polimi.it>
 *
 * This code is hereby placed in the public domain.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHORS ''AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 * OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 **/

#include "CROSS.h"
#include "csprng_hash.h"
#include "fq_arith.h"
#include "seedtree.h"
#include "merkle_tree.h"
#include <assert.h>
#include <stdio.h>
#include "pack_unpack.h"

#if defined(RSDP)
static
void expand_public_seed(FQ_ELEM V_tr[N-K][K],
                        const uint8_t seed_pub[SEED_LENGTH_BYTES]){
  CSPRNG_STATE_T CSPRNG_state_mat;
  initialize_csprng(&CSPRNG_state_mat, seed_pub, SEED_LENGTH_BYTES);
  CSPRNG_fq_mat(V_tr,&CSPRNG_state_mat);
}
#elif defined(RSDPG)
static
void expand_public_seed(FQ_ELEM V_tr[N-K][K],
                        FZ_ELEM W_tr[N-M][M],
                        const uint8_t seed_pub[SEED_LENGTH_BYTES]){
  uint8_t seedV_seed_W[2][SEED_LENGTH_BYTES];
  CSPRNG_STATE_T CSPRNG_state_mat;
  initialize_csprng(&CSPRNG_state_mat, seed_pub, SEED_LENGTH_BYTES);
  csprng_randombytes((uint8_t *)seedV_seed_W,
                     2*SEED_LENGTH_BYTES,
                     &CSPRNG_state_mat);

  initialize_csprng(&CSPRNG_state_mat, seedV_seed_W[0], SEED_LENGTH_BYTES);
  CSPRNG_fq_mat(V_tr,&CSPRNG_state_mat);

  initialize_csprng(&CSPRNG_state_mat, seedV_seed_W[1], SEED_LENGTH_BYTES);
  CSPRNG_fz_mat(W_tr,&CSPRNG_state_mat);
}
#endif


void CROSS_keygen(prikey_t *SK,
                  pubkey_t *PK){

  /* generation of random material for public and private key */
  randombytes(SK->seed,SEED_LENGTH_BYTES);
  uint8_t seede_seed_pub[2][SEED_LENGTH_BYTES];

  CSPRNG_STATE_T CSPRNG_state;
  initialize_csprng(&CSPRNG_state,SK->seed,SEED_LENGTH_BYTES);
  csprng_randombytes((uint8_t *)seede_seed_pub,
                     2*SEED_LENGTH_BYTES,
                     &CSPRNG_state);
  memcpy(PK->seed_pub,seede_seed_pub[1],SEED_LENGTH_BYTES);

  /* expansion of matrix/matrices */

  FQ_ELEM V_tr[N-K][K];
#if defined(RSDP)
  expand_public_seed(V_tr,PK->seed_pub);
#elif defined(RSDPG)
  FZ_ELEM W_tr[N-M][M];
  expand_public_seed(V_tr,W_tr,PK->seed_pub);
#endif

  /* expansion of secret key material */
  FZ_ELEM eta[N];
  CSPRNG_STATE_T CSPRNG_state_eta;
  initialize_csprng(&CSPRNG_state_eta, seede_seed_pub[0], SEED_LENGTH_BYTES);

#if defined(RSDP)
  CSPRNG_zz_vec(eta,&CSPRNG_state_eta);
#elif defined(RSDPG)
  FZ_ELEM zeta[M];
  CSPRNG_zz_inf_w(zeta,&CSPRNG_state_eta);
  fz_inf_w_by_fz_matrix(eta,zeta,W_tr);
  fz_dz_norm_sigma(eta);
#endif

  /* compute public syndrome */
  FQ_ELEM pub_syn[N-K];
  restr_vec_by_fq_matrix(pub_syn,eta,V_tr);
  fq_dz_norm_synd(pub_syn);
  pack_fq_syn(PK->s,pub_syn);
}


#if defined(RSDP)
static
void expand_private_seed(FZ_ELEM eta[N],
                         FQ_ELEM V_tr[N-K][K],
                         const uint8_t seed[SEED_LENGTH_BYTES]){
  uint8_t seede_seed_pub[2][SEED_LENGTH_BYTES];
  CSPRNG_STATE_T CSPRNG_state;
  initialize_csprng(&CSPRNG_state,seed,SEED_LENGTH_BYTES);
  csprng_randombytes((uint8_t *)seede_seed_pub,
                     2*SEED_LENGTH_BYTES,
                     &CSPRNG_state);

  expand_public_seed(V_tr,seede_seed_pub[1]);

  CSPRNG_STATE_T CSPRNG_state_eta;
  initialize_csprng(&CSPRNG_state_eta, seede_seed_pub[0], SEED_LENGTH_BYTES);
  CSPRNG_zz_vec(eta,&CSPRNG_state_eta);
}
#elif defined(RSDPG)
static
void expand_private_seed(FZ_ELEM eta[N],
                         FZ_ELEM zeta[M],
                         FQ_ELEM V_tr[N-K][K],
                         FZ_ELEM W_tr[N-M][M],
                         const uint8_t seed[SEED_LENGTH_BYTES]){
  uint8_t seede_seed_pub[2][SEED_LENGTH_BYTES];
  CSPRNG_STATE_T CSPRNG_state;
  initialize_csprng(&CSPRNG_state,seed,SEED_LENGTH_BYTES);
  csprng_randombytes((uint8_t *)seede_seed_pub,
                     2*SEED_LENGTH_BYTES,
                     &CSPRNG_state);

  expand_public_seed(V_tr,W_tr,seede_seed_pub[1]);

  CSPRNG_STATE_T CSPRNG_state_eta;
  initialize_csprng(&CSPRNG_state_eta, seede_seed_pub[0], SEED_LENGTH_BYTES);
  CSPRNG_zz_inf_w(zeta,&CSPRNG_state_eta);
  fz_inf_w_by_fz_matrix(eta,zeta,W_tr);
  fz_dz_norm_sigma(eta);
}
#endif

/* sign cannot fail */
void CROSS_sign(const prikey_t *SK,
               const char *const m,
               const uint64_t mlen,
               sig_t *sig){
    /* Key material expansion */
    FQ_ELEM V_tr[N-K][K];
    FZ_ELEM eta[N];
#if defined(RSDP)
    expand_private_seed(eta,V_tr,SK->seed);
#elif defined(RSDPG)
    FZ_ELEM zeta[M];
    FZ_ELEM W_tr[N-M][M];
    expand_private_seed(eta,zeta,V_tr,W_tr,SK->seed);
#endif

    uint8_t root_seed[SEED_LENGTH_BYTES];
    randombytes(root_seed,SEED_LENGTH_BYTES);
    randombytes(sig->salt,SALT_LENGTH_BYTES);

    uint8_t seed_tree[SEED_LENGTH_BYTES*NUM_NODES_OF_SEED_TREE];
    generate_seed_tree_from_root(seed_tree,root_seed,sig->salt);
    uint8_t * seed_tree_leaves = seed_tree +
                                 SEED_LENGTH_BYTES*(NUM_LEAVES_OF_SEED_TREE-1);

    FZ_ELEM eta_tilde[T][N];
    FZ_ELEM sigma[T][N];
    FQ_ELEM u_tilde[T][N];
    FQ_ELEM s_tilde[N-K];

#if defined(RSDP)
    uint8_t cmt_0_i_input[sizeof(FQ_ELEM)*(N-K)+
                          sizeof(FZ_ELEM)*N+
                          SALT_LENGTH_BYTES+sizeof(uint16_t)];
    const int offset_salt = sizeof(FQ_ELEM)*(N-K)+sizeof(FZ_ELEM)*N;
    const int offset_round_idx = sizeof(FQ_ELEM)*(N-K)+sizeof(FZ_ELEM)*N+SALT_LENGTH_BYTES;
#elif defined(RSDPG)
    FZ_ELEM zeta_tilde[M];
    FZ_ELEM delta[T][M];
    uint8_t cmt_0_i_input[sizeof(FQ_ELEM)*(N-K)+
                          sizeof(FZ_ELEM)*M+
                          SALT_LENGTH_BYTES+sizeof(uint16_t)];
    const int offset_salt = sizeof(FQ_ELEM)*(N-K)+sizeof(FZ_ELEM)*M;
    const int offset_round_idx = sizeof(FQ_ELEM)*(N-K)+sizeof(FZ_ELEM)*M+SALT_LENGTH_BYTES;
#endif

    uint8_t cmt_1_i_input[SEED_LENGTH_BYTES+
                          SALT_LENGTH_BYTES+sizeof(uint16_t)];

    /* cmt_0_i_input is syndrome||sigma ||salt ; place salt at the end */
    memcpy(cmt_0_i_input+offset_salt, sig->salt, SALT_LENGTH_BYTES);
    /* cmt_1_i_input is seed[i]||salt */
    memcpy(cmt_1_i_input+SEED_LENGTH_BYTES, sig->salt, SALT_LENGTH_BYTES);

    uint8_t cmt_0[T][HASH_DIGEST_LENGTH] = {0};
    uint8_t cmt_1[T][HASH_DIGEST_LENGTH] = {0};

    CSPRNG_STATE_T CSPRNG_state;
    for(int i = 0; i<T; i++){
        uint8_t seed_u_t_seed_e_t[2*SEED_LENGTH_BYTES];
        /* expand seed[i] into seed_e and seed_u */
        initialize_csprng(&CSPRNG_state,
                          seed_tree_leaves+SEED_LENGTH_BYTES*i,
                          SEED_LENGTH_BYTES);
        csprng_randombytes(seed_u_t_seed_e_t,
                           2*SEED_LENGTH_BYTES,
                           &CSPRNG_state);

        /* expand eta_tilde */
#if defined(RSDP)
        initialize_csprng(&CSPRNG_state,
                          seed_u_t_seed_e_t+SEED_LENGTH_BYTES,
                          SEED_LENGTH_BYTES);
        CSPRNG_zz_vec(eta_tilde[i], &CSPRNG_state);
#elif defined(RSDPG)
        initialize_csprng(&CSPRNG_state,
                          seed_u_t_seed_e_t+SEED_LENGTH_BYTES,
                          SEED_LENGTH_BYTES);
        CSPRNG_zz_inf_w(zeta_tilde, &CSPRNG_state);
        restr_inf_w_sub(delta[i], zeta,zeta_tilde);
        fz_dz_norm_delta(delta[i]);

        fz_inf_w_by_fz_matrix(eta_tilde[i],zeta_tilde,W_tr);
        fz_dz_norm_sigma(eta_tilde[i]);
#endif
        restr_vec_sub(sigma[i], eta,eta_tilde[i]);

        FQ_ELEM v[N];
        convert_restr_vec_to_fq(v,sigma[i]);
        fz_dz_norm_sigma(sigma[i]);
        /* expand u_tilde */
        initialize_csprng(&CSPRNG_state,
                          seed_u_t_seed_e_t,
                          SEED_LENGTH_BYTES);
        CSPRNG_fq_vec(u_tilde[i], &CSPRNG_state);

        FQ_ELEM u[N];
        fq_vec_by_fq_vec_pointwise(u,v,u_tilde[i]);
        fq_vec_by_fq_matrix(s_tilde,u,V_tr);
        fq_dz_norm_synd(s_tilde);

        /* container for s-tilde, sigma_i and salt */
        memcpy(cmt_0_i_input, s_tilde, sizeof(FQ_ELEM)*(N-K));
        const int offset_sigma_delta = sizeof(FQ_ELEM)*(N-K);

#if defined(RSDP)
        memcpy(cmt_0_i_input + offset_sigma_delta, sigma[i], sizeof(FZ_ELEM)*N);
#elif defined(RSDPG)
        memcpy(cmt_0_i_input + offset_sigma_delta, delta[i], sizeof(FZ_ELEM)*M);
#endif
        /* Fixed endianness marshalling of round counter */
        cmt_0_i_input[offset_round_idx] = (i >> 8) &0xFF;
        cmt_0_i_input[offset_round_idx+1] = i & 0xFF;
        
        hash(cmt_0[i],cmt_0_i_input,sizeof(cmt_0_i_input));
        memcpy(cmt_1_i_input,
               seed_tree_leaves+SEED_LENGTH_BYTES*i,
               SEED_LENGTH_BYTES);
        cmt_1_i_input[SEED_LENGTH_BYTES+SALT_LENGTH_BYTES] = (i >> 8) &0xFF;
        cmt_1_i_input[SEED_LENGTH_BYTES+SALT_LENGTH_BYTES+1] = i & 0xFF;        
        hash(cmt_1[i],cmt_1_i_input,sizeof(cmt_1_i_input));
    }

    /* vector containing d_0 and d_1 from spec */
    uint8_t commit_digests[2][HASH_DIGEST_LENGTH];
    uint8_t merkle_tree_0[NUM_NODES_OF_MERKLE_TREE * HASH_DIGEST_LENGTH];

    merkle_tree_root_compute(commit_digests[0], merkle_tree_0, cmt_0);
    hash(commit_digests[1], (unsigned char*)cmt_1, sizeof(cmt_1));
    hash(sig->digest_01,
              (unsigned char*) commit_digests,
              sizeof(commit_digests));

    /* first challenge extraction */
    uint8_t beta_buf[2*HASH_DIGEST_LENGTH+SALT_LENGTH_BYTES];
    /* place d_m at the beginning of the input of the hash generating d_beta*/
    hash(beta_buf,(uint8_t*) m,mlen);
    memcpy(beta_buf+HASH_DIGEST_LENGTH, sig->digest_01, HASH_DIGEST_LENGTH);
    memcpy(beta_buf+2*HASH_DIGEST_LENGTH, sig->salt, SALT_LENGTH_BYTES);

    uint8_t d_beta[HASH_DIGEST_LENGTH];
    hash(d_beta,beta_buf,2*HASH_DIGEST_LENGTH+SALT_LENGTH_BYTES);

    FQ_ELEM beta[T];
    initialize_csprng(&CSPRNG_state,d_beta,HASH_DIGEST_LENGTH);
    CSPRNG_fq_vec_beta(beta, &CSPRNG_state);

    /* Computation of the first round of responses */
    FQ_ELEM y[T][N];
    for(int i = 0; i < T; i++){
        fq_vec_by_restr_vec_scaled(y[i],
                                   eta_tilde[i],
                                   beta[i],
                                   u_tilde[i]);
        fq_dz_norm(y[i]);
    }

    /* Second challenge extraction */
    uint8_t digest_b_buf[T*N*sizeof(FQ_ELEM)+HASH_DIGEST_LENGTH];
    memcpy(digest_b_buf,y,T*N*sizeof(FQ_ELEM));
    memcpy(digest_b_buf+T*N*sizeof(FQ_ELEM),d_beta,HASH_DIGEST_LENGTH);
    hash(sig->digest_b, digest_b_buf, sizeof(digest_b_buf));

    uint8_t fixed_weight_b[T]={0};
    /*since w > t-w, we generate a weight t-w string and flip the contents  */
    expand_digest_to_fixed_weight(fixed_weight_b,sig->digest_b);
    for(int i = 0; i<T; i++){
        fixed_weight_b[i] = !fixed_weight_b[i];
    }

    /* Computation of the second round of responses */
    uint16_t mtp_len;
    merkle_tree_proof_compute(sig->mtp,&mtp_len,merkle_tree_0,cmt_0,fixed_weight_b);
    publish_seeds(sig->stp,seed_tree,fixed_weight_b);

    int published_rsps = 0;
    for(int i = 0; i<T; i++){
        if(fixed_weight_b[i] == 0){
            assert(published_rsps < T-W);
            pack_fq_vec(sig->rsp_0[published_rsps].y, y[i]);
#if defined(RSDP)
            pack_fz_vec(sig->rsp_0[published_rsps].sigma, sigma[i]);
#elif defined(RSDPG)
            pack_fz_rsdp_g_vec(sig->rsp_0[published_rsps].delta, delta[i]);
#endif
            memcpy(sig->rsp_1[published_rsps], cmt_1[i], HASH_DIGEST_LENGTH);
            published_rsps++;
        }
    }
}

/* verify returns 1 if signature is ok, 0 otherwise */
int CROSS_verify(const pubkey_t *const PK,
                 const char *const m,
                 const uint64_t mlen,
                 const sig_t *const sig){
    CSPRNG_STATE_T CSPRNG_state;

    FQ_ELEM V_tr[N-K][K];
#if defined(RSDP)
    expand_public_seed(V_tr,PK->seed_pub);
#elif defined(RSDPG)
    FZ_ELEM W_tr[N-M][M];
    expand_public_seed(V_tr,W_tr,PK->seed_pub);
#endif

    FQ_ELEM pub_syn[N-K];
    unpack_fq_syn(pub_syn,PK->s);

    uint8_t beta_buf[2*HASH_DIGEST_LENGTH+SALT_LENGTH_BYTES];
    hash(beta_buf,(uint8_t*) m,mlen);
    memcpy(beta_buf+HASH_DIGEST_LENGTH, sig->digest_01, HASH_DIGEST_LENGTH);
    memcpy(beta_buf+2*HASH_DIGEST_LENGTH, sig->salt, SALT_LENGTH_BYTES);

    uint8_t d_beta[HASH_DIGEST_LENGTH];
    hash(d_beta,beta_buf,sizeof(beta_buf));

    FQ_ELEM beta[T];
    initialize_csprng(&CSPRNG_state,d_beta,HASH_DIGEST_LENGTH);
    CSPRNG_fq_vec_beta(beta, &CSPRNG_state);

    uint8_t fixed_weight_b[T]={0};
    expand_digest_to_fixed_weight(fixed_weight_b,sig->digest_b);
    /* flip fixed weight string, obtaining T-W set bits, change after
     * parameter convention becomes flipped */
    for(int i = 0; i<T; i++){
        fixed_weight_b[i] = !fixed_weight_b[i];
    }

    uint8_t seed_tree[SEED_LENGTH_BYTES*NUM_NODES_OF_SEED_TREE];
    regenerate_leaves(seed_tree, fixed_weight_b, sig->stp, sig->salt);

    uint8_t * seed_tree_leaves = seed_tree +
                                 SEED_LENGTH_BYTES*(NUM_LEAVES_OF_SEED_TREE-1);

#if defined(RSDP)
    uint8_t cmt_0_i_input[sizeof(FQ_ELEM)*(N-K)+
                          sizeof(FZ_ELEM)*N+
                          SALT_LENGTH_BYTES+sizeof(uint16_t)];
    const int offset_salt = sizeof(FQ_ELEM)*(N-K)+sizeof(FZ_ELEM)*N;
    const int offset_round_idx = sizeof(FQ_ELEM)*(N-K)+sizeof(FZ_ELEM)*N+SALT_LENGTH_BYTES;    
#elif defined(RSDPG)
    uint8_t cmt_0_i_input[sizeof(FQ_ELEM)*(N-K)+
                          sizeof(FZ_ELEM)*M+
                          SALT_LENGTH_BYTES+sizeof(uint16_t)];
    const int offset_salt = sizeof(FQ_ELEM)*(N-K)+sizeof(FZ_ELEM)*M;
    const int offset_round_idx = sizeof(FQ_ELEM)*(N-K)+sizeof(FZ_ELEM)*M+SALT_LENGTH_BYTES;    
#endif
    /* cmt_0_i_input is syndrome||sigma ||salt */
    memcpy(cmt_0_i_input+offset_salt, sig->salt, SALT_LENGTH_BYTES);

    uint8_t cmt_1_i_input[SEED_LENGTH_BYTES+SALT_LENGTH_BYTES+sizeof(uint16_t)];
    memcpy(cmt_1_i_input+SEED_LENGTH_BYTES, sig->salt, SALT_LENGTH_BYTES);
    /* cmt_1_i_input is seed[i]||salt */

    uint8_t cmt_0[T][HASH_DIGEST_LENGTH] = {0};
    uint8_t cmt_1[T][HASH_DIGEST_LENGTH] = {0};

    uint8_t seed_u_t_seed_e_t[2*SEED_LENGTH_BYTES];
    FZ_ELEM eta_tilde[N];
    FQ_ELEM u_tilde[N];

    FQ_ELEM y_tilde[N] = {0};
    FQ_ELEM s_tilde[N-K] = {0};

    FQ_ELEM y[T][N];

    int used_rsps = 0;
    int is_signature_ok = 1;
    for(int i = 0; i< T; i++){
        if(fixed_weight_b[i] == 1){
            memcpy(cmt_1_i_input,
                   seed_tree_leaves+SEED_LENGTH_BYTES*i,
                   SEED_LENGTH_BYTES);
            cmt_1_i_input[SEED_LENGTH_BYTES+SALT_LENGTH_BYTES] = (i >> 8) &0xFF;
            cmt_1_i_input[SEED_LENGTH_BYTES+SALT_LENGTH_BYTES+1] = i & 0xFF; 
            hash(cmt_1[i],cmt_1_i_input,sizeof(cmt_1_i_input));
            /* expand seed[i] into seed_e and seed_u */
            initialize_csprng(&CSPRNG_state,
                              seed_tree_leaves+SEED_LENGTH_BYTES*i,
                              SEED_LENGTH_BYTES);
            csprng_randombytes(seed_u_t_seed_e_t,
                               2*SEED_LENGTH_BYTES,
                               &CSPRNG_state);
#if defined(RSDP)
            /* expand eta_tilde */
            initialize_csprng(&CSPRNG_state,
                              seed_u_t_seed_e_t+SEED_LENGTH_BYTES,
                              SEED_LENGTH_BYTES);
            CSPRNG_zz_vec(eta_tilde, &CSPRNG_state);
#elif defined(RSDPG)
            initialize_csprng(&CSPRNG_state,
                              seed_u_t_seed_e_t+SEED_LENGTH_BYTES,
                              SEED_LENGTH_BYTES);
            FZ_ELEM zeta_tilde[M];
            CSPRNG_zz_inf_w(zeta_tilde, &CSPRNG_state);
            fz_inf_w_by_fz_matrix(eta_tilde,zeta_tilde,W_tr);
            fz_dz_norm_sigma(eta_tilde);
#endif
            /* expand u_tilde */
            initialize_csprng(&CSPRNG_state,
                              seed_u_t_seed_e_t,
                              SEED_LENGTH_BYTES);
            CSPRNG_fq_vec(u_tilde, &CSPRNG_state);
            fq_vec_by_restr_vec_scaled(y[i],
                                       eta_tilde,
                                       beta[i],
                                       u_tilde);
            fq_dz_norm(y[i]);
        } else {
            /* place y[i] in the buffer for later on hashing */
            unpack_fq_vec(y[i], sig->rsp_0[used_rsps].y);

            FZ_ELEM sigma_local[N];
#if defined(RSDP)
            /*sigma is memcpy'ed directly into cmt_0 input buffer */
            FZ_ELEM* sigma_ptr = cmt_0_i_input+sizeof(FQ_ELEM)*(N-K);
	    unpack_fz_vec(sigma_local, sig->rsp_0[used_rsps].sigma);
            memcpy(sigma_ptr, sigma_local, sizeof(FZ_ELEM)*N);
            is_signature_ok = is_signature_ok &&
                              is_zz_vec_in_restr_group(sigma_local);
#elif defined(RSDPG)
            /*delta is memcpy'ed directly into cmt_0 input buffer */
            FZ_ELEM* sigma_ptr = cmt_0_i_input+sizeof(FQ_ELEM)*(N-K);
	    unpack_fz_rsdp_g_vec(sigma_ptr, sig->rsp_0[used_rsps].delta);
            is_signature_ok = is_signature_ok &&
                              is_zz_inf_w_valid(sigma_ptr);
            fz_inf_w_by_fz_matrix(sigma_local,sigma_ptr,W_tr);

#endif
            memcpy(cmt_1[i], sig->rsp_1[used_rsps], HASH_DIGEST_LENGTH);
            used_rsps++;

            FQ_ELEM v[N];
            convert_restr_vec_to_fq(v,sigma_local);
            fq_vec_by_fq_vec_pointwise(y_tilde,v,y[i]);

            fq_vec_by_fq_matrix(s_tilde,y_tilde,V_tr);
            fq_dz_norm_synd(s_tilde);
            fq_synd_minus_fq_vec_scaled((FQ_ELEM*) cmt_0_i_input,
                                        s_tilde,
                                        beta[i],
                                        pub_syn);
            fq_dz_norm_synd((FQ_ELEM*) cmt_0_i_input);
            cmt_0_i_input[offset_round_idx] = (i >> 8) &0xFF;
            cmt_0_i_input[offset_round_idx+1] = i & 0xFF;            
            
            hash(cmt_0[i], cmt_0_i_input, sizeof(cmt_0_i_input));
        }
    } /* end for iterating on ZKID iterations */


    uint8_t commit_digests[2][HASH_DIGEST_LENGTH];
    merkle_tree_root_recompute(commit_digests[0],
                               cmt_0,
                               sig->mtp,
                               fixed_weight_b);

    uint8_t digest_01_recomputed[HASH_DIGEST_LENGTH];
    hash(commit_digests[1], (unsigned char*)cmt_1, sizeof(cmt_1));
    hash(digest_01_recomputed,
              (unsigned char*) commit_digests,
              sizeof(commit_digests));

    uint8_t digest_b_recomputed[HASH_DIGEST_LENGTH];

    uint8_t digest_b_buf[T*N*sizeof(FQ_ELEM)+HASH_DIGEST_LENGTH];
    memcpy(digest_b_buf,y,T*N*sizeof(FQ_ELEM));
    memcpy(digest_b_buf+T*N*sizeof(FQ_ELEM),d_beta,HASH_DIGEST_LENGTH);
    hash(digest_b_recomputed, digest_b_buf, sizeof(digest_b_buf));

    int does_digest_01_match = ( memcmp(digest_01_recomputed,
                                        sig->digest_01,
                                        HASH_DIGEST_LENGTH) == 0);
    int does_digest_b_match = ( memcmp(digest_b_recomputed,
                                        sig->digest_b,
                                        HASH_DIGEST_LENGTH) == 0);
    is_signature_ok = is_signature_ok &&
                      does_digest_01_match &&
                      does_digest_b_match;
    return is_signature_ok;
}
