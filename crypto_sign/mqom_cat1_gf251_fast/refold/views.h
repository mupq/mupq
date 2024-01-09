#ifndef MQOM_VIEWS_H
#define MQOM_VIEWS_H

#include <stdint.h>
#include <stddef.h>
#include "parameters-all.h"

/**
 * @brief Derive a hash digest into
 *    "nb_sets" sets of "nb_parties_by_set" parties.
 *    There is no redundancy in a set of parties.
 *    A party is designed by its number between
 *         0 and PARAM_NB_PARTIES.
 *    The output is the flattened list of these sets.
 *    Each set is sorted according to the party numbers.
 * 
 * @param digest the hash digest to transform.
 *              Must be a string of PARAM_DIGEST_SIZE bytes.
 * @param nb_sets the number of sets.
 * @param nb_parties_by_set the number of parties in each set.
 * @param opened_views the flattened list of the party set.
 *              Must be a uint16_t array of length "nb_sets*nb_parties_by_set"
 */
void expand_view_challenge_hash(uint16_t* opened_views, const uint8_t* digest, unsigned int nb_sets, unsigned int nb_parties_by_set);

#endif /* MQOM_VIEWS_H */
