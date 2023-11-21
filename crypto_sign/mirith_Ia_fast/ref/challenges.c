/*
 * Copyright 2023 Carlo Sanna, Javier Verbel, and Floyd Zweydinger.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "challenges.h"

void get_second_challenges(uint32_t i_star[TAU], const hash_t hash2)
{
    size_t l;
    prng_t prng;
    
    prng_init(&prng, hash2, NULL);
    
    for (l = 0; l < TAU; l++)
    {
        uint32_t  r;
        
        prng_bytes(&prng, &r, sizeof(r));
        
        /* NOTE: It is well-known that this method to generate
         * a random integer in [0, N_PARTIES - 1] is slighly biased
         * toward small values. */
        i_star[l] = r % N_PARTIES;
    }
}
