
#include "views.h"
#include "hash.h"
#include "mpc-all.h"

void expand_view_challenge_hash(uint16_t* opened_views, const uint8_t* digest, unsigned int nb_sets, unsigned int nb_parties_by_set) {

    xof_context entropy_ctx;
    xof_init(&entropy_ctx);
    xof_update(&entropy_ctx, digest, PARAM_DIGEST_SIZE);
    samplable_t entropy = xof_to_samplable(&entropy_ctx);

    uint16_t ind = 0;
    uint16_t mask = (1<<PARAM_LOG_NB_PARTIES) - 1;
    uint8_t tmp[2];
    uint16_t value;
    uint8_t is_unique;
    for(uint32_t num=0; num<nb_sets; num++) {
        for(uint32_t i=0; i<nb_parties_by_set; i++) {
            do {
                do {
                    byte_sample(&entropy, tmp, sizeof(uint16_t));
                    value = ((uint16_t) tmp[0]) | ((uint16_t) (tmp[1]<<8)); // Deal with endianness
                    value &= mask;
                } while(value >= PARAM_NB_PARTIES);

                // Check if unique
                is_unique = 1;
                for(size_t j=0; j<i && is_unique; j++)
                    if(opened_views[ind-i+j] == value)
                        is_unique = 0;
            } while(!is_unique);
            opened_views[ind++] = value;
        }
    }

    // Sort the open views
    //    Have no effect when nb_parties_by_set = 1
    for(unsigned int e=0; e<nb_sets; e++) {
        for(size_t p=0; p<nb_parties_by_set-1; p++) {
            for(size_t i=p+1; i<nb_parties_by_set; i++) {
                if(opened_views[nb_parties_by_set*e+p] > opened_views[nb_parties_by_set*e+i]) {
                    value = opened_views[nb_parties_by_set*e+i];
                    opened_views[nb_parties_by_set*e+i] = opened_views[nb_parties_by_set*e+p];
                    opened_views[nb_parties_by_set*e+p] = value;
                }
            }
        }
    }
}

