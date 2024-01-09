#ifndef MQOM_SERIALIZATION_H
#define MQOM_SERIALIZATION_H

#include "parameters-all.h"
#include "mpc-all.h"

// Witness
#ifndef PARAM_WIT_SHORT_SIZE
static inline void wit_serialize(uint8_t* buf, const mpc_wit_t* wit) {
    vec_set(buf, wit, PARAM_WIT_SIZE);
}
static inline void wit_deserialize(mpc_wit_t const** wit, const uint8_t* buf, const uint8_t* memory) {
    (void) memory; // No need to use additional memory
    (*wit) = (const mpc_wit_t*) buf;
}
#define PARAM_WIT_SHORT_SIZE PARAM_WIT_SIZE
#define PARAM_WIT_PARSED_SIZE (0)
#else
#define wit_serialize(buf, wit) wit_compress(buf, wit)
static inline void wit_deserialize(mpc_wit_t const** wit, const uint8_t* buf, uint8_t* memory) {
    wit_decompress((mpc_wit_t*) memory, buf);
    (*wit) = (const mpc_wit_t*) memory;
}
#define PARAM_WIT_PARSED_SIZE PARAM_WIT_SIZE
#endif

// Unif
#ifndef PARAM_UNIF_SHORT_SIZE
static inline void unif_serialize(uint8_t* buf, const mpc_unif_t* unif) {
    vec_set(buf, unif, PARAM_UNIF_SIZE);
}
static inline void unif_deserialize(mpc_unif_t const** unif, const uint8_t* buf, const uint8_t* memory) {
    (void) memory; // No need to use additional memory
    (*unif) = (const mpc_unif_t*) buf;
}
#define PARAM_UNIF_SHORT_SIZE PARAM_UNIF_SIZE
#define PARAM_UNIF_PARSED_SIZE (0)
#else
#define unif_serialize(buf, unif) unif_compress(buf, unif)
static inline void unif_deserialize(mpc_unif_t const** unif, const uint8_t* buf, uint8_t* memory) {
    unif_decompress((mpc_unif_t*) memory, buf);
    (*unif) = (const mpc_unif_t*) memory;
}
#define PARAM_UNIF_PARSED_SIZE PARAM_UNIF_SIZE
#endif

// Hint
#ifndef PARAM_HINT_SHORT_SIZE
static inline void hint_serialize(uint8_t* buf, const mpc_hint_t* hint) {
    vec_set(buf, hint, PARAM_HINT_SIZE);
}
static inline void hint_deserialize(mpc_hint_t const** hint, const uint8_t* buf, const uint8_t* memory) {
    (void) memory; // No need to use additional memory
    (*hint) = (const mpc_hint_t*) buf;
}
#define PARAM_HINT_SHORT_SIZE PARAM_HINT_SIZE
#define PARAM_HINT_PARSED_SIZE (0)
#else
#define hint_serialize(buf, hint) hint_compress(buf, hint)
static inline void hint_deserialize(mpc_hint_t const** hint, const uint8_t* buf, uint8_t* memory) {
    hint_decompress((mpc_hint_t*) memory, buf);
    (*hint) = (const mpc_hint_t*) memory;
}
#define PARAM_HINT_PARSED_SIZE PARAM_HINT_SIZE
#endif

// Broadcast
#ifndef PARAM_BR_SHORT_SIZE
static inline void br_serialize(uint8_t* buf, const mpc_broadcast_t* br) {
    vec_set(buf, br, PARAM_BR_SIZE);
}
static inline void br_deserialize(mpc_broadcast_t const** br, const uint8_t* buf, const uint8_t* memory) {
    (void) memory; // No need to use additional memory
    (*br) = (const mpc_broadcast_t*) buf;
}
#define PARAM_BR_SHORT_SIZE PARAM_BR_SIZE
#define PARAM_BR_PARSED_SIZE (0)
#else
#define br_serialize(buf, br) br_compress(buf, br)
static inline void br_deserialize(mpc_broadcast_t const** br, const uint8_t* buf, uint8_t* memory) {
    br_decompress((mpc_broadcast_t*) memory, buf);
    (*br) = (const mpc_broadcast_t*) memory;
}
#define PARAM_BR_PARSED_SIZE PARAM_BR_SIZE
#endif

#endif /* MQOM_SERIALIZATION_H */
