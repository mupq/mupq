#include "sign-mpcith.h"
#include "sample.h"
#include "tree.h"
#include "commit.h"
#include "hash.h"
#include "keygen.h"
#include "mpc.h"
#include "views.h"
#include "parameters-all.h"
#include "serialization.h"
#include <string.h>
#include <stdio.h>

/***********************************************
 *            Signature Structure              *
 ***********************************************/

typedef struct proof_hypercube_7r_t {
    uint8_t* seed_info; // Information required to compute the tree with seeds of of all opened parties
    uint8_t* unopened_digest;
    uint8_t* hint_digest;

    // Plaintext broadcast messages
    mpc_broadcast_t* plain_broadcast;
    
    // Last party's share
    mpc_wit_t* wit;
    mpc_hint_t* hint;
} proof_hypercube_7r_t;

typedef struct signature_hypercube_7r_t {
    uint8_t* salt;
    uint8_t* mpc_challenge_1_hash;
    uint8_t* mpc_challenge_2_hash;
    uint8_t* view_challenge_hash;
    proof_hypercube_7r_t proofs[PARAM_NB_EXECUTIONS];
    uint8_t* allocated_memory; // Just to manage the memory
} signature_hypercube_7r_t;

// For parsing
typedef struct const_proof_hypercube_7r_t {
    const uint8_t* seed_info; // Information required to compute the tree with seeds of of all opened parties
    const uint8_t* unopened_digest;
    const uint8_t* hint_digest;

    // Plaintext broadcast messages
    mpc_broadcast_t* plain_broadcast;
    
    // Last party's share
    const mpc_wit_t* wit;
    const mpc_hint_t* hint;
} const_proof_hypercube_7r_t;

typedef struct const_signature_hypercube_7r_t {
    const uint8_t* salt;
    const uint8_t* mpc_challenge_1_hash;
    const uint8_t* mpc_challenge_2_hash;
    const uint8_t* view_challenge_hash;
    const_proof_hypercube_7r_t proofs[PARAM_NB_EXECUTIONS];
    uint8_t* allocated_memory; // Just to manage the memory
} const_signature_hypercube_7r_t;

// Free signature structure
void free_signature(signature_hypercube_7r_t* sig);
void free_const_signature(const_signature_hypercube_7r_t* sig);

// For signing
signature_hypercube_7r_t* init_signature_structure(const uint8_t* salt, uint8_t* buf, size_t buflen);
int build_signature(const signature_hypercube_7r_t* sig, uint8_t* buf, size_t buflen, uint16_t hidden_views[PARAM_NB_EXECUTIONS]);

// For verification
const_signature_hypercube_7r_t* parse_signature(const uint8_t* buf, size_t buflen, uint16_t hidden_views[PARAM_NB_EXECUTIONS]);

/***********************************************
 *             Hash for challenge              *
 ***********************************************/

static void hash_for_mpc_challenge_1(uint8_t challenge_hash[PARAM_DIGEST_SIZE], uint8_t seed_commitments[PARAM_NB_EXECUTIONS][PARAM_NB_PARTIES][PARAM_DIGEST_SIZE],
    const instance_t* inst, const uint8_t* salt, const uint8_t* message, size_t mlen)
{
    hash_context ctx;
    hash_init_prefix(&ctx, HASH_PREFIX_FIRST_CHALLENGE);
    if(inst != NULL)
        hash_update_instance(&ctx, inst);
    if(mlen > 0)
        hash_update(&ctx, message, mlen);
    hash_update(&ctx, salt, PARAM_SALT_SIZE);
    hash_update(&ctx, (uint8_t*) seed_commitments, PARAM_NB_EXECUTIONS*PARAM_NB_PARTIES*PARAM_DIGEST_SIZE);
    hash_final(&ctx, challenge_hash);
}

static void hash_for_mpc_challenge_2(uint8_t challenge_hash[PARAM_DIGEST_SIZE], const uint8_t mpc_challenge_hash_1[PARAM_DIGEST_SIZE], uint8_t hint_commitments[PARAM_NB_EXECUTIONS][PARAM_DIGEST_SIZE],
    const uint8_t* salt, const uint8_t* message, size_t mlen)
{
    hash_context ctx;
    hash_init_prefix(&ctx, HASH_PREFIX_SECOND_CHALLENGE);
    if(mlen > 0)
        hash_update(&ctx, message, mlen);
    hash_update(&ctx, salt, PARAM_SALT_SIZE);
    hash_update(&ctx, mpc_challenge_hash_1, PARAM_DIGEST_SIZE);
    hash_update(&ctx, (uint8_t*) hint_commitments, PARAM_NB_EXECUTIONS*PARAM_DIGEST_SIZE);
    hash_final(&ctx, challenge_hash);
}

static void hash_for_view_challenge(uint8_t challenge_hash[PARAM_DIGEST_SIZE], const uint8_t mpc_challenge_hash_2[PARAM_DIGEST_SIZE],
    mpc_broadcast_t* broadcast[PARAM_NB_EXECUTIONS][PARAM_HYPERCUBE_DIMENSION], mpc_broadcast_t* plain_broadcast[], 
    const uint8_t* salt, const uint8_t* message, size_t mlen)
{
    hash_context ctx;
    hash_init_prefix(&ctx, HASH_PREFIX_THIRD_CHALLENGE);
    if(mlen > 0)
        hash_update(&ctx, message, mlen);
    hash_update(&ctx, salt, PARAM_SALT_SIZE);
    hash_update(&ctx, mpc_challenge_hash_2, PARAM_DIGEST_SIZE);
    for (size_t e=0; e<PARAM_NB_EXECUTIONS; e++)
        for(size_t i=0; i<PARAM_HYPERCUBE_DIMENSION; i++)
            hash_update(&ctx, (uint8_t*) broadcast[e][i], PARAM_BR_SIZE);
    for (size_t e=0; e<PARAM_NB_EXECUTIONS; e++)
        hash_update(&ctx, (uint8_t*) plain_broadcast[e], PARAM_UNIF_SIZE); // We do not commit the v's zero values
    hash_final(&ctx, challenge_hash);
}

/***********************************************
 *             Signing Algorithm               *
 ***********************************************/

int mpcith_hypercube_7r_sign(uint8_t* sig, size_t* siglen,
                const uint8_t* m, size_t mlen,
                const uint8_t* sk,
                const uint8_t* salt, const uint8_t* seed
                ) {
    uint32_t i, j, e, p;
    int ret;

#ifndef NDEBUG
    // In debug mode, let us check if the key generation
    //    produces valid key pair. 
    ret = mqom_validate_keys(NULL, sk);
    if(ret)
        printf("Error: SK invalid for signing.\n");
#endif

    // Deserialize the private key
    mqom_secret_key_t ssk;
    ret = deserialize_secret_key(&ssk, sk, PARAM_SECRETKEYBYTES);
    if(ret < 0)
        return ret;
    uncompress_instance(ssk.inst);
    mpc_wit_t* plain_wit = ssk.wit;

    // Signature Structure
    signature_hypercube_7r_t* ssig = init_signature_structure(salt, sig, PARAM_SIGNATURE_SIZEBYTES);

    /********************************************/
    /**********     INITIALIZATION     **********/
    /********************************************/

    prg_context entropy_ctx;
    samplable_t entropy = prg_to_samplable(&entropy_ctx);
    #ifdef PARAM_RND_EXPANSION_X4
    prg_context_x4 entropy_ctx_x4;
    samplable_x4_t entropy_x4 = prg_to_samplable_x4(&entropy_ctx_x4);
    #else /* PARAM_RND_EXPANSION_X4 */
    prg_context entropy_ctx_x4[4];
    samplable_t entropy_x4[4];
    for(j=0; j<4; j++)
        entropy_x4[j]= prg_to_samplable(&entropy_ctx_x4[j]);
    #endif /* PARAM_RND_EXPANSION_X4 */
    seed_tree_t* seeds_tree[PARAM_NB_EXECUTIONS];
    for(e=0; e<PARAM_NB_EXECUTIONS; e++)
        seeds_tree[e] = malloc_seed_tree(PARAM_HYPERCUBE_DIMENSION);
    uint8_t seed_commitments[PARAM_NB_EXECUTIONS][PARAM_NB_PARTIES][PARAM_DIGEST_SIZE];
    uint8_t hint_commitments[PARAM_NB_EXECUTIONS][PARAM_DIGEST_SIZE];

    // Derive the root seeds for all parallel executions
    uint8_t rseed[PARAM_NB_EXECUTIONS][PARAM_SEED_SIZE];
    prg_init(&entropy_ctx, seed, ssig->salt);
    byte_sample(&entropy, (uint8_t*) rseed, PARAM_NB_EXECUTIONS*PARAM_SEED_SIZE);

    // We manipulate here pointers
    //    to avoid & when using
    //    to be more consistant with plain_wit
    mpc_share_t* mshares[PARAM_NB_EXECUTIONS][PARAM_HYPERCUBE_DIMENSION];
    for(e=0; e<PARAM_NB_EXECUTIONS; e++)
        for(p=0; p<PARAM_HYPERCUBE_DIMENSION; p++) {
            mshares[e][p] = new_share();
            vec_set_zero(mshares[e][p], PARAM_SHARE_SIZE);
        }
    mpc_unif_t* plain_unif[PARAM_NB_EXECUTIONS];
    mpc_hint_t* plain_hint[PARAM_NB_EXECUTIONS];

    /********************************************/
    /********   COMMIT PARTIES' INPUTS   ********/
    /********************************************/

    mpc_share_t* shares[4];
    for(j=0; j<4; j++)
        shares[j] = new_share();
    mpc_share_t* share = shares[0];
    mpc_share_t* acc = new_share();
    for(e=0; e<PARAM_NB_EXECUTIONS; e++) {
        // Build the seed tree of the current execution
        expand_seed_tree(seeds_tree[e], rseed[e], ssig->salt);
        uint8_t** seeds = get_leaves(seeds_tree[e]);
        vec_set_zero(acc, PARAM_SHARE_SIZE);

        i=0;
        // Let us treat the parties four by four...
        for(; (i+3)+1<PARAM_NB_PARTIES; i+=4) {
            const uint8_t* ptr_seeds[4] = {seeds[i], seeds[i+1], seeds[i+2], seeds[i+3]};
            #ifdef PARAM_RND_EXPANSION_X4
            prg_init_x4(&entropy_ctx_x4, ptr_seeds, NULL);
            vec_rnd_x4((void**) shares, PARAM_SHARE_SIZE, &entropy_x4);
            #else /* PARAM_RND_EXPANSION_X4 */
            // No need of the salt, we do not care about
            //    collisions of the leave seeds
            prg_init_x4_array(entropy_ctx_x4, ptr_seeds, NULL);
            for(j=0; j<4; j++)
                vec_rnd((void*) shares[j], PARAM_SHARE_SIZE, &entropy_x4[j]);
            #endif /* PARAM_RND_EXPANSION_X4 */
            for(j=0; j<4; j++) {
                vec_add(acc, shares[j], PARAM_SHARE_SIZE);
                for(p=0; p<PARAM_HYPERCUBE_DIMENSION; p++)
                    if((((i+j)>>p) & 1) == 0) // If pth bit of i is zero
                        vec_add(mshares[e][p], shares[j], PARAM_SHARE_SIZE);
            }
            uint8_t* digests[4] = {
                seed_commitments[e][i+0], seed_commitments[e][i+1],
                seed_commitments[e][i+2], seed_commitments[e][i+3], 
            };
            const uint16_t is[4] = {(uint16_t)i, (uint16_t)(i+1), (uint16_t)(i+2), (uint16_t)(i+3)};
            commit_seed_x4(digests, (uint8_t const*const*) &seeds[i], ssig->salt, (uint16_t)e, is);
        }
        // Now we treat the last parties
        for(;i<PARAM_NB_PARTIES; i++) {
            // No need of the salt, we do not care about
            //    collisions of the leave seeds
            prg_init(&entropy_ctx, seeds[i], NULL);
            if(i != PARAM_NB_PARTIES-1) {
                // Expand the input share from seed
                vec_rnd(share, PARAM_SHARE_SIZE, &entropy);
                vec_add(acc, share, PARAM_SHARE_SIZE);

                // Aggregate to get the shares of the "main parties"
                for(p=0; p<PARAM_HYPERCUBE_DIMENSION; p++)
                    if(((i>>p) & 1) == 0) // If pth bit of i is zero
                        vec_add(mshares[e][p], share, PARAM_SHARE_SIZE);
            
                // Commit to the party's input (by committing to its seed)
                commit_seed(seed_commitments[e][i], seeds[i], ssig->salt, (uint16_t)e, (uint16_t)i);
            } else {
                // Compute plain unif
                plain_unif[e] = new_unif();
                vec_rnd(plain_unif[e], PARAM_UNIF_SIZE, &entropy);
                vec_add(plain_unif[e], get_unif(acc), PARAM_UNIF_SIZE);
                vec_normalize(plain_unif[e], PARAM_UNIF_SIZE);

                // For the moment, ssig->proofs[e].hint just contains the hint accumulator
                vec_set(ssig->proofs[e].hint, get_hint(acc), PARAM_HINT_SIZE);

                // Compute aux
                vec_set(ssig->proofs[e].wit, plain_wit, PARAM_WIT_SIZE);
                vec_sub(ssig->proofs[e].wit, get_wit(acc), PARAM_WIT_SIZE);
                vec_normalize(ssig->proofs[e].wit, PARAM_WIT_SIZE);

                // Commit to the party's input
                commit_seed_and_wit(seed_commitments[e][i], seeds[i], ssig->proofs[e].wit, ssig->salt, (uint16_t)e, (uint16_t)i);
            }
        }
    }
    free(acc);
    for(j=0; j<4; j++)
        free(shares[j]);


    // Expand the MPC challenge N°1
    mpc_challenge_1_t* mpc_challenges_1[PARAM_NB_EXECUTIONS];
    for(e=0; e<PARAM_NB_EXECUTIONS; e++)
        mpc_challenges_1[e] = new_challenge_1();
    hash_for_mpc_challenge_1(ssig->mpc_challenge_1_hash, seed_commitments, ssk.inst, ssig->salt, NULL, 0);
    expand_mpc_challenge_hash_1(mpc_challenges_1, ssig->mpc_challenge_1_hash, PARAM_NB_EXECUTIONS, ssk.inst);

    /********************************************/
    /*********    COMPUTE THE HINTS    **********/
    /********************************************/

    for(e=0; e<PARAM_NB_EXECUTIONS; e++) {
        plain_hint[e] = new_hint();
        compute_hint(plain_hint[e], plain_wit, plain_unif[e], ssk.inst, mpc_challenges_1[e]);
        vec_neg(ssig->proofs[e].hint, PARAM_HINT_SIZE);
        vec_add(ssig->proofs[e].hint, plain_hint[e], PARAM_HINT_SIZE);
        commit_hint(hint_commitments[e], ssig->proofs[e].hint, ssig->salt, (uint16_t)e, PARAM_NB_PARTIES);
        memcpy(ssig->proofs[e].hint_digest, hint_commitments[e], PARAM_DIGEST_SIZE);
    }

    // Expand the MPC challenge N°2
    mpc_challenge_2_t* mpc_challenges_2[PARAM_NB_EXECUTIONS];
    for(e=0; e<PARAM_NB_EXECUTIONS; e++)
        mpc_challenges_2[e] = new_challenge_2();
    hash_for_mpc_challenge_2(ssig->mpc_challenge_2_hash, ssig->mpc_challenge_1_hash, hint_commitments, ssig->salt, NULL, 0);
    expand_mpc_challenge_hash_2(mpc_challenges_2, ssig->mpc_challenge_2_hash, PARAM_NB_EXECUTIONS, ssk.inst);

    /********************************************/
    /*********  SIMULATE MPC PROTOCOL  **********/
    /********************************************/

    mpc_broadcast_t* plain_broadcast[PARAM_NB_EXECUTIONS];
    mpc_broadcast_t* broadcast[PARAM_NB_EXECUTIONS][PARAM_HYPERCUBE_DIMENSION];
    for(e=0; e<PARAM_NB_EXECUTIONS; e++) {
        mpc_share_t* plain = new_share();
        vec_set(get_wit(plain), plain_wit, PARAM_WIT_SIZE);
        vec_set(get_unif(plain), plain_unif[e], PARAM_UNIF_SIZE);
        vec_set(get_hint(plain), plain_hint[e], PARAM_HINT_SIZE);

        // Compute the values open by the parties
        plain_broadcast[e] = new_br();
        mpc_compute_plain_broadcast(plain_broadcast[e], mpc_challenges_1[e], mpc_challenges_2[e], plain, ssk.inst);
        vec_set(ssig->proofs[e].plain_broadcast, plain_broadcast[e], PARAM_BR_SIZE);
    
        // Compute the broadcast messages of the main parties
        for(p=0; p<PARAM_HYPERCUBE_DIMENSION; p++) {
            broadcast[e][p] = new_br();
            vec_normalize(mshares[e][p], PARAM_SHARE_SIZE);
            mpc_compute_communications(broadcast[e][p], mpc_challenges_1[e], mpc_challenges_2[e], mshares[e][p], ssig->proofs[e].plain_broadcast, ssk.inst, 0);
        }

        free(plain);
    }

    // Expand the view challenge
    uint16_t hidden_views[PARAM_NB_EXECUTIONS];
    hash_for_view_challenge(ssig->view_challenge_hash, ssig->mpc_challenge_2_hash, broadcast, plain_broadcast, ssig->salt, m, mlen);
    expand_view_challenge_hash((uint16_t*) hidden_views, ssig->view_challenge_hash, PARAM_NB_EXECUTIONS, 1);

    /********************************************/
    /**********   FINALIZE SIGNATURE   **********/
    /********************************************/

    for(e=0; e<PARAM_NB_EXECUTIONS; e++) {
        uint16_t num_unopened_party = hidden_views[e];
        get_seed_path(ssig->proofs[e].seed_info, seeds_tree[e], num_unopened_party);
        memcpy(ssig->proofs[e].unopened_digest, seed_commitments[e][num_unopened_party], PARAM_DIGEST_SIZE);
    }

    // Serialize the signature
    ret = build_signature(ssig, sig, PARAM_SIGNATURE_SIZEBYTES, hidden_views);
    free_signature(ssig);

    for(e=0; e<PARAM_NB_EXECUTIONS; e++)
        free_seed_tree(seeds_tree[e]);
    mqom_free_keys_internal(NULL, &ssk);
    for(e=0; e<PARAM_NB_EXECUTIONS; e++) {
        for(p=0; p<PARAM_HYPERCUBE_DIMENSION; p++) {
            free(mshares[e][p]);
            free(broadcast[e][p]);
        }
        free(plain_unif[e]);
        free(plain_hint[e]);
        free(plain_broadcast[e]);
        free(mpc_challenges_1[e]);
        free(mpc_challenges_2[e]);
    }

    if(ret < 0)
        return ret;
    *siglen = (size_t) ret;
    return 0;
}

/***********************************************
 *            Verication Algorithm             *
 ***********************************************/

int mpcith_hypercube_7r_sign_verify(const uint8_t* sig, size_t siglen,
                const uint8_t* m, size_t mlen,
                const uint8_t* pk
                ) {
    uint32_t i, j, e, p;
    int ret;

    // Deserialize the public key
    mqom_public_key_t ppk;
    ret = deserialize_public_key(&ppk, pk, PARAM_PUBLICKEYBYTES);
    if (ret < 0)
        return ret;
    uncompress_instance(ppk.inst);

    // Parse the signature
    // Note: while parsing, it expands the view challenge.
    uint16_t hidden_views[PARAM_NB_EXECUTIONS];
    const_signature_hypercube_7r_t* ssig = parse_signature(sig, siglen, hidden_views);
    if(ssig == NULL) {
        ret = -1;
        mqom_free_keys_internal(&ppk, NULL);
        return ret;
    }

    // Initialize
    uint8_t seed_commitments[PARAM_NB_EXECUTIONS][PARAM_NB_PARTIES][PARAM_DIGEST_SIZE];
    uint8_t hint_commitments[PARAM_NB_EXECUTIONS][PARAM_DIGEST_SIZE];
    prg_context entropy_ctx;
    samplable_t entropy = prg_to_samplable(&entropy_ctx);
    #ifdef PARAM_RND_EXPANSION_X4
    prg_context_x4 entropy_ctx_x4;
    samplable_x4_t entropy_x4 = prg_to_samplable_x4(&entropy_ctx_x4);
    #else /* PARAM_RND_EXPANSION_X4 */
    prg_context entropy_ctx_x4[4];
    samplable_t entropy_x4[4];
    for(j=0; j<4; j++)
        entropy_x4[j]= prg_to_samplable(&entropy_ctx_x4[j]);
    #endif /* PARAM_RND_EXPANSION_X4 */

    // Expand MPC Challenge
    mpc_challenge_1_t* mpc_challenges_1[PARAM_NB_EXECUTIONS];
    for(e=0; e<PARAM_NB_EXECUTIONS; e++)
        mpc_challenges_1[e] = new_challenge_1();
    expand_mpc_challenge_hash_1(mpc_challenges_1, ssig->mpc_challenge_1_hash, PARAM_NB_EXECUTIONS, ppk.inst);

    mpc_challenge_2_t* mpc_challenges_2[PARAM_NB_EXECUTIONS];
    for(e=0; e<PARAM_NB_EXECUTIONS; e++)
        mpc_challenges_2[e] = new_challenge_2();
    expand_mpc_challenge_hash_2(mpc_challenges_2, ssig->mpc_challenge_2_hash, PARAM_NB_EXECUTIONS, ppk.inst);

    mpc_share_t* shares[4];
    for(j=0; j<4; j++)
        shares[j] = new_share();
    mpc_share_t* share = shares[0];
    mpc_share_t* mshares[PARAM_HYPERCUBE_DIMENSION];
    for(p=0; p<PARAM_HYPERCUBE_DIMENSION; p++)
        mshares[p] = new_share();
    mpc_broadcast_t* plain_broadcast[PARAM_NB_EXECUTIONS];
    mpc_broadcast_t* broadcast[PARAM_NB_EXECUTIONS][PARAM_HYPERCUBE_DIMENSION];
    for(e=0; e<PARAM_NB_EXECUTIONS; e++) {
        plain_broadcast[e] = new_br();
        for(p=0; p<PARAM_HYPERCUBE_DIMENSION; p++)
            broadcast[e][p] = new_br();
    }

    for(e=0; e<PARAM_NB_EXECUTIONS; e++) {
        // Get the open leaf seeds
        seed_tree_t* seed_tree = malloc_seed_tree(PARAM_HYPERCUBE_DIMENSION);
        reconstruct_tree(seed_tree, hidden_views[e], ssig->proofs[e].seed_info, ssig->salt);
        uint8_t** seeds = get_leaves(seed_tree);

        // Get the commitment of the hint
        if(hidden_views[e] == PARAM_NB_PARTIES-1)
            memcpy(hint_commitments[e], ssig->proofs[e].hint_digest, PARAM_DIGEST_SIZE);
        else
            commit_hint(hint_commitments[e], ssig->proofs[e].hint, ssig->salt, (uint16_t)e, PARAM_NB_PARTIES);

        // Get the plain broadcast
        vec_set(plain_broadcast[e], ssig->proofs[e].plain_broadcast, PARAM_BR_SIZE);

        // Get the share of the main parties
        for(p=0; p<PARAM_HYPERCUBE_DIMENSION; p++)
            vec_set_zero(mshares[p], PARAM_SHARE_SIZE);

        i=0;
        // Let us treat the parties four by four...
        for(; (i+3)+1<PARAM_NB_PARTIES; i+=4) {
            uint8_t* digests[4] = {
                seed_commitments[e][i+0], seed_commitments[e][i+1],
                seed_commitments[e][i+2], seed_commitments[e][i+3], 
            };
            const uint16_t is[4] = {(uint16_t)i, (uint16_t)(i+1), (uint16_t)(i+2), (uint16_t)(i+3)};
            commit_seed_x4(digests, (uint8_t const*const*) &seeds[i], ssig->salt, (uint16_t)e, is);
            const uint8_t* ptr_seeds[4] = {seeds[i], seeds[i+1], seeds[i+2], seeds[i+3]};
            #ifdef PARAM_RND_EXPANSION_X4
            prg_init_x4(&entropy_ctx_x4, ptr_seeds, NULL);
            vec_rnd_x4((void**) shares, PARAM_SHARE_SIZE, &entropy_x4);
            #else /* PARAM_RND_EXPANSION_X4 */
            // No need of the salt, we do not care about
            //    collisions of the leave seeds
            prg_init_x4_array(entropy_ctx_x4, ptr_seeds, NULL);
            for(j=0; j<4; j++)
                vec_rnd((void*) shares[j], PARAM_SHARE_SIZE, &entropy_x4[j]);
            #endif /* PARAM_RND_EXPANSION_X4 */
            for(j=0; j<4; j++) {
                if(i+j != hidden_views[e]) {
                    // Aggregate to get the shares of the "main parties"
                    uint16_t xor = (((uint16_t)(i+j)) ^ hidden_views[e]);
                    for(p=0; p<PARAM_HYPERCUBE_DIMENSION; p++)
                        if(((xor>>p) & 1) == 1) // If the pth bit of i* is different of the pth bit of i
                            vec_add(mshares[p], shares[j], PARAM_SHARE_SIZE);
                } else {
                    memcpy(seed_commitments[e][i+j], ssig->proofs[e].unopened_digest, PARAM_DIGEST_SIZE);
                }
            }
        }
        // Now we treat the last parties
        for(; i<PARAM_NB_PARTIES; i++) {
            if(i == hidden_views[e]) {
                memcpy(seed_commitments[e][i], ssig->proofs[e].unopened_digest, PARAM_DIGEST_SIZE);
                continue;
            }
            prg_init(&entropy_ctx, seeds[i], NULL);
            if(i != PARAM_NB_PARTIES-1) {
                // Expand the input share from seed
                vec_rnd(share, PARAM_SHARE_SIZE, &entropy);

                // Recompute the party's commitment
                commit_seed(seed_commitments[e][i], seeds[i], ssig->salt, (uint16_t)e, (uint16_t)i);
            } else {
                // Get the party's input share
                vec_rnd(get_unif(share), PARAM_UNIF_SIZE, &entropy);
                vec_set(get_wit(share), ssig->proofs[e].wit, PARAM_WIT_SIZE);
                vec_set(get_hint(share), ssig->proofs[e].hint, PARAM_HINT_SIZE);

                // Recompute the party's commitment
                commit_seed_and_wit(seed_commitments[e][i], seeds[i], ssig->proofs[e].wit, ssig->salt, (uint16_t)e, (uint16_t)i);
            }

            // Aggregate to get the shares of the "main parties"
            uint16_t xor = (((uint16_t)i) ^ hidden_views[e]);
            for(p=0; p<PARAM_HYPERCUBE_DIMENSION; p++)
                if(((xor>>p) & 1) == 1) // If the pth bit of i* is different of the pth bit of i
                    vec_add(mshares[p], share, PARAM_SHARE_SIZE);
        }

        // Re-emulate the computation of the open main parties
        for(p=0; p<PARAM_HYPERCUBE_DIMENSION; p++) {
            vec_normalize(mshares[p], PARAM_SHARE_SIZE);
            if(((hidden_views[e]>>p) & 1) == 1) {
                // It the index of the hidden mparty is set,
                //    it implies that the current mshare is for the mparty with unset bit.
                mpc_compute_communications(broadcast[e][p], mpc_challenges_1[e], mpc_challenges_2[e], mshares[p], ssig->proofs[e].plain_broadcast, ppk.inst, 0);
            } else {
                mpc_compute_communications(broadcast[e][p], mpc_challenges_1[e], mpc_challenges_2[e], mshares[p], ssig->proofs[e].plain_broadcast, ppk.inst, 1);
                vec_neg(broadcast[e][p], PARAM_BR_SIZE);
                vec_add(broadcast[e][p], ssig->proofs[e].plain_broadcast, PARAM_BR_SIZE);
                vec_normalize(broadcast[e][p], PARAM_BR_SIZE);
            }
        }

        free_seed_tree(seed_tree);
    }
    for(j=0; j<4; j++)
        free(shares[j]);
    for(p=0; p<PARAM_HYPERCUBE_DIMENSION; p++)
        free(mshares[p]);

    // Recompute the hash digests of the challenges
    //    and check they are consistent with the ones in the signature
    uint8_t mpc_challenge_1_hash[PARAM_DIGEST_SIZE], mpc_challenge_2_hash[PARAM_DIGEST_SIZE], view_challenge_hash[PARAM_DIGEST_SIZE];
    hash_for_mpc_challenge_1(mpc_challenge_1_hash, seed_commitments, ppk.inst, ssig->salt, NULL, 0);
    hash_for_mpc_challenge_2(mpc_challenge_2_hash, ssig->mpc_challenge_1_hash, hint_commitments, ssig->salt, NULL, 0);
    hash_for_view_challenge(view_challenge_hash, ssig->mpc_challenge_2_hash, broadcast, plain_broadcast, ssig->salt, m, mlen);
    ret = (memcmp(mpc_challenge_1_hash, ssig->mpc_challenge_1_hash, PARAM_DIGEST_SIZE) != 0);
    ret |= (memcmp(mpc_challenge_2_hash, ssig->mpc_challenge_2_hash, PARAM_DIGEST_SIZE) != 0);
    ret |= (memcmp(view_challenge_hash, ssig->view_challenge_hash, PARAM_DIGEST_SIZE) != 0);

    mqom_free_keys_internal(&ppk, NULL);
    free_const_signature(ssig);
    for(e=0; e<PARAM_NB_EXECUTIONS; e++) {
        for(p=0; p<PARAM_HYPERCUBE_DIMENSION; p++)
            free(broadcast[e][p]);
        free(plain_broadcast[e]);
        free(mpc_challenges_1[e]);
        free(mpc_challenges_2[e]);
    }
    return ret;
}

/***********************************************
 *             About Serialization             *
 ***********************************************/

void free_signature(signature_hypercube_7r_t* sig) {
    if(sig->allocated_memory != NULL)
        free(sig->allocated_memory);
    free(sig);
}

void free_const_signature(const_signature_hypercube_7r_t* sig) {
    if(sig->allocated_memory != NULL)
        free(sig->allocated_memory);
    free(sig);
}

const_signature_hypercube_7r_t* parse_signature(const uint8_t* buf, size_t buflen, uint16_t hidden_views[PARAM_NB_EXECUTIONS]) {
    size_t bytes_required = 3*PARAM_DIGEST_SIZE + PARAM_SALT_SIZE;
    if(buflen < bytes_required) {
        return NULL;
    }

    const_signature_hypercube_7r_t* sig = malloc(sizeof(const_signature_hypercube_7r_t));
    sig->salt = buf;                  buf += PARAM_SALT_SIZE;
    sig->mpc_challenge_1_hash = buf;  buf += PARAM_DIGEST_SIZE;
    sig->mpc_challenge_2_hash = buf;  buf += PARAM_DIGEST_SIZE;
    sig->view_challenge_hash = buf;   buf += PARAM_DIGEST_SIZE;
    expand_view_challenge_hash((uint16_t*) hidden_views, sig->view_challenge_hash, PARAM_NB_EXECUTIONS, 1);

    for(size_t e=0; e<PARAM_NB_EXECUTIONS; e++) {
        bytes_required += PARAM_SEED_SIZE*PARAM_HYPERCUBE_DIMENSION + PARAM_DIGEST_SIZE + PARAM_COMPRESSED_BR_SIZE;
        if(hidden_views[e] != PARAM_NB_PARTIES-1)
            bytes_required += PARAM_WIT_SHORT_SIZE + PARAM_HINT_SHORT_SIZE;
        else
            bytes_required += PARAM_DIGEST_SIZE;
    }
    if(buflen != bytes_required) {
        free(sig);
        return NULL;
    }

    sig->allocated_memory = malloc(
        PARAM_NB_EXECUTIONS*(
            PARAM_BR_SIZE+PARAM_WIT_PARSED_SIZE+PARAM_HINT_PARSED_SIZE
        )
    );
    uint8_t* buf_ = sig->allocated_memory;
    for(size_t e=0; e<PARAM_NB_EXECUTIONS; e++) {
        sig->proofs[e].unopened_digest = buf;
        buf += PARAM_DIGEST_SIZE;
        sig->proofs[e].plain_broadcast = (void*) buf_;
        buf_ += PARAM_BR_SIZE;
        uncompress_plain_broadcast(sig->proofs[e].plain_broadcast, buf);
        buf += PARAM_COMPRESSED_BR_SIZE;
        sig->proofs[e].seed_info = buf;
        buf += PARAM_SEED_SIZE*PARAM_HYPERCUBE_DIMENSION;
        if(hidden_views[e] != PARAM_NB_PARTIES-1) {
            wit_deserialize(&sig->proofs[e].wit, buf, buf_);
            buf += PARAM_WIT_SHORT_SIZE;
            buf_ += PARAM_WIT_PARSED_SIZE;
            hint_deserialize(&sig->proofs[e].hint, buf, buf_);
            buf += PARAM_HINT_SHORT_SIZE;
            buf_ += PARAM_HINT_PARSED_SIZE;
        } else {
            sig->proofs[e].hint_digest = buf;
            buf += PARAM_DIGEST_SIZE;
        }
    }

    return sig;
}

signature_hypercube_7r_t* init_signature_structure(const uint8_t* salt, uint8_t* buf, size_t buflen) {
    size_t bytes_required = 3*PARAM_DIGEST_SIZE + PARAM_SALT_SIZE;
    if(buflen < bytes_required) {
        return NULL;
    }

    signature_hypercube_7r_t* sig = malloc(sizeof(signature_hypercube_7r_t));
    sig->salt = buf;               buf += PARAM_SALT_SIZE;
    sig->mpc_challenge_1_hash = buf;   buf += PARAM_DIGEST_SIZE;
    sig->mpc_challenge_2_hash = buf;   buf += PARAM_DIGEST_SIZE;
    sig->view_challenge_hash = buf;  buf += PARAM_DIGEST_SIZE;
    memcpy(sig->salt, salt, PARAM_SALT_SIZE);

    sig->allocated_memory = malloc(
        PARAM_NB_EXECUTIONS*(
            PARAM_HYPERCUBE_DIMENSION*PARAM_SEED_SIZE // seed tree
             + PARAM_DIGEST_SIZE // unopen commitment
             + PARAM_BR_SIZE // plain broadcast
             + PARAM_WIT_SIZE // wit of last party
             + PARAM_HINT_SIZE // hint of last party
             + PARAM_DIGEST_SIZE // hint commitment
        )
    );
    uint8_t* buf_ = sig->allocated_memory;
    for(size_t e=0; e<PARAM_NB_EXECUTIONS; e++) {
        sig->proofs[e].unopened_digest = buf_;
        buf_ += PARAM_DIGEST_SIZE;
        sig->proofs[e].plain_broadcast = (void*) buf_;
        buf_ += PARAM_BR_SIZE;
        sig->proofs[e].seed_info = buf_;
        buf_ += PARAM_HYPERCUBE_DIMENSION*PARAM_SEED_SIZE;
        sig->proofs[e].wit = (void*) buf_;   buf_ += PARAM_WIT_SIZE;
        sig->proofs[e].hint = (void*) buf_;  buf_ += PARAM_HINT_SIZE;
        sig->proofs[e].hint_digest = buf_;   buf_ += PARAM_DIGEST_SIZE;
    }

    return sig;
}

int build_signature(const signature_hypercube_7r_t* sig, uint8_t* buf, size_t buflen, uint16_t hidden_views[PARAM_NB_EXECUTIONS]) {
    size_t bytes_required = (3*PARAM_DIGEST_SIZE + PARAM_SALT_SIZE);
    bytes_required += PARAM_NB_EXECUTIONS*PARAM_COMPRESSED_BR_SIZE;
    bytes_required += PARAM_NB_EXECUTIONS*PARAM_DIGEST_SIZE;
    for(size_t e=0; e<PARAM_NB_EXECUTIONS; e++) {
        bytes_required += PARAM_SEED_SIZE*PARAM_HYPERCUBE_DIMENSION;
        if(hidden_views[e] != PARAM_NB_PARTIES-1)
            bytes_required += PARAM_WIT_SHORT_SIZE + PARAM_HINT_SHORT_SIZE;
        else
            bytes_required += PARAM_DIGEST_SIZE;
    }
    if(buflen < bytes_required)
        return -1;

    uint8_t* bufBase = buf;
    buf += (3*PARAM_DIGEST_SIZE + PARAM_SALT_SIZE); // salt and challenge hashes are already set (via init)
    for(size_t e=0; e<PARAM_NB_EXECUTIONS; e++) {
        memcpy(buf, sig->proofs[e].unopened_digest, PARAM_DIGEST_SIZE);
        buf += PARAM_DIGEST_SIZE;
        compress_plain_broadcast(buf, sig->proofs[e].plain_broadcast);
        buf += PARAM_COMPRESSED_BR_SIZE;
        memcpy(buf, sig->proofs[e].seed_info, PARAM_SEED_SIZE*PARAM_HYPERCUBE_DIMENSION);
        buf += PARAM_SEED_SIZE*PARAM_HYPERCUBE_DIMENSION;
        if(hidden_views[e] != PARAM_NB_PARTIES-1) {
            wit_serialize(buf, sig->proofs[e].wit);
            buf += PARAM_WIT_SHORT_SIZE;
            hint_serialize(buf, sig->proofs[e].hint);
            buf += PARAM_HINT_SHORT_SIZE;
        } else {
            memcpy(buf, sig->proofs[e].hint_digest, PARAM_DIGEST_SIZE);
            buf += PARAM_DIGEST_SIZE;
        }
    }

    return (int)(buf - bufBase);
}
