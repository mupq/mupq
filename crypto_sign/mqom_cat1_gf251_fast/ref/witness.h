#ifndef MQOM_WITNESS_H
#define MQOM_WITNESS_H

#include <stdint.h>
#include "hash.h"
#include "sample.h"

#include <stdint.h>
#include "parameters.h"

// Instance Definition
//   This structure represents an instance of
//   the problem on which the security of the
//   signature scheme relies. It corresponds
//   to the public key.
//   Some member can be pointers when they are
//   generated at each signing and verification
//   from the others members. 
typedef struct instance_t {
    uint8_t seed[PARAM_SEED_SIZE];
    uint8_t y[PARAM_m];
    uint8_t (*A)[PARAM_m][PARAM_MATRIX_BYTESIZE];
    uint8_t (*b)[PARAM_m][PARAM_n];
} instance_t;

// Solution Definition
//   This structure represents a solution for
//   an instance presented by "instance_t".
//   It is part of the secret key of the
//   signature scheme.
//   It corresponds to the extended solution,
//   meaning that it contains all the secret
//   values which can be deterministically
//   built from the solution itself and which
//   are inputs of the underlying MPC protocol.
typedef struct solution_t {
    uint8_t x[PARAM_n];
} solution_t;

// "PARAM_INSTANCE_SIZE" corresponds to the size
//   in bytes of the structure of type "instance_t"
//   when serialized.
#define PARAM_INSTANCE_SIZE (PARAM_SEED_SIZE+PARAM_m)

// "PARAM_SOL_SIZE" corresponds to the size
//   in bytes of the structure of type "solution_t"
//   when serialized.
#define PARAM_SOL_SIZE (PARAM_n)




/**
 * @brief Generate a random pair (inst, sol) where
 *    "inst" is a problem instance and "sol" is one of
 *    its extended solution.
 * 
 * @param inst a pointer to a pointer which will point where a memory block
 *     containing a problem instance after the call to this function.
 * @param sol a pointer to a pointer which will point where a memory block
 *     containing a solution after the call to this function.
 * @param entropy a entropy source of type "samplable_t"
 */
void generate_instance_with_solution(instance_t** inst, solution_t** sol, samplable_t* entropy);

/**
 * @brief Check if "sol" is a valid extended solution of
 *    the instance represented by "inst".
 * 
 * @param inst a problem instance
 * @param sol the tested solution
 * @return 1 if it is the case, 0 otherwise
 */
int is_correct_solution(instance_t* inst, solution_t* sol);

/**
 * @brief Check if "inst1" and "inst2" corresponds to
 *    the same problem instance.
 * 
 * @param inst1 the first tested instance
 * @param inst2 the second tested instance
 * @return 1 if it is the case, 0 otherwise
 */
int are_same_instances(instance_t* inst1, instance_t* inst2);

/**
 * @brief Serialize a instance structure in a byte buffer
 * 
 * @param buf the buffer (allocated array of PARAM_INSTANCE_SIZE bytes)
 * @param inst the instance to serialized
 */
void serialize_instance(uint8_t* buf, const instance_t* inst);

/**
 * @brief Deserialize the instance structure which has been
 *    serialized in a byte buffer by "serialize_instance"
 * 
 * @param buf the buffer containing the instance to deserialize
 * @return the deserialized instance
 */
instance_t* deserialize_instance(const uint8_t* buf);

/**
 * @brief Update a hash context by putting in the hash input
 *      the instance structure
 * 
 * @param ctx the hash structure
 * @param inst the instance to hash
 */
void hash_update_instance(hash_context* ctx, const instance_t* inst);

/**
 * @brief Serialize a solution structure in a byte buffer
 * 
 * @param buf the buffer (allocated array of PARAM_SOL_SIZE bytes)
 * @param sol the solution to serialized
 */
void serialize_instance_solution(uint8_t* buf, const solution_t* sol);

/**
 * @brief Deserialize the solution structure which has been
 *    serialized in a byte buffer by "serialize_instance_solution"
 * 
 * @param buf the buffer containing the solution to deserialize
 * @return the deserialized solution
 */
solution_t* deserialize_instance_solution(const uint8_t* buf);

/**
 * @brief Uncompress a instance structure. If some of the members
 *    is still undefined, compute them from the other members.
 * 
 * @param inst the instance to uncompress   
 */
void uncompress_instance(instance_t* inst);

/**
 * @brief Deallocates the memory previously allocated
 *    for a solution by a call to generate_instance_with_solution
 *    or deserialize_instance_solution.
 * 
 * @param sol the pointer to a memory block previously allocated for a solution.
 */
void free_instance_solution(solution_t* sol);

/**
 * @brief Deallocates the memory previously allocated
 *    for a instance by a call to generate_instance_with_solution
 *    or deserialize_instance.
 *
 * @param inst the pointer to a memory block previously allocated for a instance.
 */
void free_instance(instance_t* inst);

#endif /* MQOM_WITNESS_H */
