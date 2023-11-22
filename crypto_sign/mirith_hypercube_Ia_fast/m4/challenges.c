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

void get_second_challenges(int i_star[TAU], const hash_t hash2)
{
    int l;
    prng_t prng;
    
    prng_init(&prng, hash2, NULL);
    
    for (l = 0; l < TAU; l++)
    {
        int r;
        
        prng_bytes(&prng, &r, sizeof(r));
        
        /* NOTE: It is well-known that this method to generate
         * a random integer in [0, N_PARTIES_ROUND - 1] is slighly biased
         * toward small values. */
        i_star[l] = r % N_PARTIES_ROUND;
        
        if (i_star[l] < 0)
        {
            i_star[l] += N_PARTIES_ROUND;
        }
        /* * */
    }
}
