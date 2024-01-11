#include "snova.h"

#include "gf16_matrix_inline.h"
#include "fips202.h"
#include "aes.h"

static gf16m_t S[l_SNOVA] = {0};

/**
 * Generate elements of F16[S]
 */
void gen_S_array() {
    be_aI(S[0], 1);
    be_the_S(S[1]);
    for (int index = 2; index < l_SNOVA; index++) {
        gf16m_mul(S[index - 1], S[1], S[index]);
    }
}

/**
 * SNOVA init
 */
void snova_init() {
    init_gf16_tables();
    gen_S_array();
}


/**
 * Using AES-CTR encryption as a hash function
 * AES ciphertext padded with zeros.
 * The iv is also padded with zeros.
 * Using input value as the AES key.
 * The ciphertext obtained from AES encryption serves as the output of the hash
 * function.
 * @param pt_seed_array - Pointer to the hash input. (Fixed length of 16)
 * @param pt_output_array - Pointer to the hash output. (Fixed length of
 * bytes_prng_public)
 */
void hash_aes128(uint8_t* pt_seed_array, uint8_t* pt_output_array) {
    const uint8_t iv[12] = {0};
    aes128ctx ctx;
    aes128_ctr_keyexp(&ctx, pt_seed_array);
    aes128_ctr(pt_output_array, bytes_prng_public, iv, &ctx);
}

/**
 * Convert one byte of data to GF16 representation (using only half of the
 * byte). Example: <bytes 12 34 56 78 9a bc> -> <bytes 02 01 04 03 05 ..... 0c
 * 0b>
 * @param byte_array - input (bytes)
 * @param gf16_array - output (GF16)
 * @param num_of_GF16s - GF16 amount
 */
void convert_bytes_to_GF16s(uint8_t* byte_array, gf16_t* gf16_array,
                            int num_of_GF16s) {
    uint16_t* GF16_array_2B = (uint16_t*)gf16_array;
    for (int index = 0; index < num_of_GF16s >> 1; index++) {
        *(GF16_array_2B++) =
            (byte_array[index] & 0xF) ^ ((byte_array[index] & 0xF0) << 4);
    }
    if (num_of_GF16s % 2 == 1) {
        gf16_array[num_of_GF16s - 1] = byte_array[num_of_GF16s >> 1] & 0xF;
    }
}

/**
 * Convert two GF16 values to one byte.
 * Example:
 *  <bytes 02 01 04 03 05 ..... 0c 0b> -> <bytes 12 34 56 78 9a bc>
 * @param byte_array - output (bytes)
 * @param gf16_array - input (GF16)
 * @param num_of_GF16s - GF16 amount
 */
void convert_GF16s_to_bytes(uint8_t* byte_array, gf16_t* gf16_array,
                            int num_of_GF16s) {
    for (int index = 0; index < num_of_GF16s / 2; index++) {
        byte_array[index] =
            gf16_array[index << 1] | (gf16_array[(index << 1) + 1] << 4);
    }
    if (num_of_GF16s % 2 == 1)
        byte_array[num_of_GF16s / 2] = gf16_array[num_of_GF16s - 1];
}

/**
 * Convert one byte of data to GF16 representation (using only half of the
 * byte). cut_in_half Example: <bytes 12 34 56 78 9a bc> -> <bytes 02 04 06 08
 * 0a 01 03 05 07 09 0b>
 * @param byte_array - input (bytes)
 * @param gf16_array - output (GF16)
 * @param num_of_GF16s - GF16 amount
 */
void convert_bytes_to_GF16s_cut_in_half(uint8_t* byte_array, gf16_t* gf16_array,
                                        int num_of_GF16s) {
    uint64_t* pt_GF16_8;
    uint64_t* pt_GF16half_8;
    uint64_t* pt_byte_8;
    pt_GF16_8 = (uint64_t*)(gf16_array + ((num_of_GF16s + 1) >> 1));
    pt_GF16half_8 = (uint64_t*)gf16_array;
    pt_byte_8 = (uint64_t*)byte_array;
    for (int index = (num_of_GF16s >> 4); index > 0; --index) {
        *(pt_GF16half_8++) = ((*(pt_byte_8))) & 0x0F0F0F0F0F0F0F0F;
        *(pt_GF16_8++) = ((*pt_byte_8++) >> 4) & 0x0F0F0F0F0F0F0F0F;
    }
    for (int index = (num_of_GF16s >> 4) << 3; index < (num_of_GF16s >> 1);
         ++index) {
        gf16_array[index + ((num_of_GF16s + 1) >> 1)] = byte_array[index] >> 4;
    }
    for (int index = (num_of_GF16s >> 4) << 3;
         index < ((num_of_GF16s + 1) >> 1); ++index) {
        gf16_array[index] = byte_array[index] & 0xF;
    }
}

/**
 * Convert two GF16 values to one byte.
 * Example:
 *  <bytes 02 04 06 08 0a 01 03 05 07 09 0b> -> <bytes 12 34 56 78 9a bc>
 * @param byte_array - output (bytes)
 * @param gf16_array - input (GF16)
 * @param num_of_GF16s - GF16 amount
 */
void convert_GF16s_to_bytes_merger_in_half(uint8_t* byte_array,
                                           gf16_t* gf16_array,
                                           int num_of_GF16s) {
    uint64_t* pt_GF16_8;
    uint64_t* pt_GF16half_8;
    uint64_t* pt_byte_8;
    pt_GF16_8 = (uint64_t*)(gf16_array + ((num_of_GF16s + 1) >> 1));
    pt_GF16half_8 = (uint64_t*)gf16_array;
    pt_byte_8 = (uint64_t*)byte_array;
    for (int index = (num_of_GF16s >> 4); index > 0; --index) {
        (*(pt_byte_8++)) = (*(pt_GF16_8++) << 4) | (*(pt_GF16half_8++));
    }
    for (int index = (num_of_GF16s >> 4) << 3; index < (num_of_GF16s >> 1);
         ++index) {
        byte_array[index] =
            gf16_array[index] |
            ((gf16_array[index + ((num_of_GF16s + 1) >> 1)]) << 4);
    }
    if (num_of_GF16s & 1 == 1) {
        byte_array[num_of_GF16s >> 1] = gf16_array[num_of_GF16s >> 1];
    }
}

/**
 * @param c - output
 * @param pt_matrix - input
 */
void gen_a_FqS(gf16_t* c, gf16m_t* pt_matrix) {
    gf16m_t temp;
    be_aI(*pt_matrix, c[0]);
    for (int i = 1; i < rank - 1; ++i) {
        gf16m_scale(S[i], c[i], temp);
        gf16m_add(*pt_matrix, temp, *pt_matrix);
    }
    gf16m_scale(S[rank - 1],
                (c[rank - 1] != 0) ? c[rank - 1] : 16 - (c[0] + (c[0] == 0)),
                temp);
    gf16m_add(*pt_matrix, temp, *pt_matrix);
}

/**
 * Generate the linear map T12
 * @param T12 - output
 * @param seed - input
 */
void gen_seeds_and_T12(T12_t T12, uint8_t* seed) {
    gf16_t* pt_array;
    uint8_t prng_output_private[bytes_prng_private];
    gf16_t GF16_prng_output_private[GF16s_prng_private];

    shake256(prng_output_private, bytes_prng_private, seed, seed_length_private);
    convert_bytes_to_GF16s(prng_output_private, GF16_prng_output_private,
                           GF16s_prng_private);

    // cout << "T12======" << endl;
    pt_array = GF16_prng_output_private;
    for (int j = 0; j < v_SNOVA; ++j) {
        for (int k = 0; k < o_SNOVA; ++k) {
            gen_a_FqS(pt_array, T12[j] + k);
            // --- T12[j][k].printout();
            // printf("T12[%d][%d] = \n", j, k);
            // print_gf16m(key->map2.T12[j][k]);
            pt_array += l_SNOVA;
        }
    }
}

/**
 * Generate the random part of public key
 * @param map - P11 P12 P21 Aalpha Balpha Qalpha1 Qalpha2
 * @param pt_public_key_seed - input
 */
void gen_A_B_Q_P(map_group1* map, uint8_t* pt_public_key_seed) {
    uint8_t prng_output_public[bytes_prng_public];
    uint8_t temp[lsq_SNOVA * l_SNOVA];
    // ----- pt temp -----
    hash_aes128(pt_public_key_seed, prng_output_public);
    convert_bytes_to_GF16s(prng_output_public, (uint8_t*)map,
                           GF16s_prng_public);
    memcpy(temp, ((uint8_t*)map->Qalpha1), lsq_SNOVA * l_SNOVA);
    // printf("Aalpha: \n");
    for (int index = 0; index < lsq_SNOVA; ++index) {
        be_invertible_by_add_aS(map->Aalpha[index]);
        // print_byte(map->Aalpha[index], sq_rank);
        // Aalpha[index].printout();
    }

    // cout << "B======" << endl;
    for (int index = 0; index < lsq_SNOVA; ++index) {
        be_invertible_by_add_aS(map->Balpha[index]);
        // Balpha[index].printout();
    }
    gf16_t* pt_array = (uint8_t*)(map->Qalpha1) + l_SNOVA * lsq_SNOVA;
    for (int index = 0; index < lsq_SNOVA; ++index) {
        gen_a_FqS(pt_array, map->Qalpha2 + index);
        be_invertible_by_add_aS(map->Qalpha2[index]);
        pt_array += l_SNOVA;
        // Qalpha2[index].printout();
    }

    pt_array = temp;
    for (int index = 0; index < lsq_SNOVA; ++index) {
        gen_a_FqS(pt_array, map->Qalpha1 + index);
        be_invertible_by_add_aS(map->Qalpha1[index]);
        // Qalpha1[index].printout();
        pt_array += l_SNOVA;
    }
}

/**
 * Generate private key (F part)
 * @param map2 - output: F11 F12 F21
 * @param map1 - input: P11 P12 P21 Aalpha Balpha Qalpha1 Qalpha2
 * @param T12 - input
 */
void gen_F(map_group2* map2, map_group1* map1, T12_t T12) {
    gf16m_t temp;
    for (int i = 0; i < m_SNOVA; ++i) {
        for (int j = 0; j < v_SNOVA; ++j) {
            for (int k = 0; k < v_SNOVA; ++k) {
                gf16m_clone(map2->F11[i][j][k], map1->P11[i][j][k]);
            }
        }
    }

    for (int i = 0; i < m_SNOVA; ++i) {
        for (int j = 0; j < v_SNOVA; ++j) {
            for (int k = 0; k < o_SNOVA; ++k) {
                gf16m_clone(map2->F12[i][j][k], map1->P12[i][j][k]);
                for (int index = 0; index < v_SNOVA; ++index) {
                    gf16m_mul(map1->P11[i][j][index], T12[index][k], temp);
                    gf16m_add(map2->F12[i][j][k], temp, map2->F12[i][j][k]);
                    // F12[i][j][k].printout();
                }
            }
        }
    }

    for (int i = 0; i < m_SNOVA; ++i) {
        for (int j = 0; j < o_SNOVA; ++j) {
            for (int k = 0; k < v_SNOVA; ++k) {
                gf16m_clone(map2->F21[i][j][k], map1->P21[i][j][k]);
                for (int index = 0; index < v_SNOVA; ++index) {
                    gf16m_mul(T12[index][j], map1->P11[i][index][k], temp);
                    gf16m_add(map2->F21[i][j][k], temp, map2->F21[i][j][k]);
                    // F21[i][j][k].printout();
                }
            }
        }
    }
}

/**
 * Generate public key (P22 part)
 * @param T12 - input
 * @param P21 - input
 * @param F12 - input
 * @param outP22 - output
 */
void gen_P22(T12_t T12, P21_t P21, F12_t F12, P22_byte_t outP22) {
    gf16m_t temp1, temp2;
    P22_t P22 = {0};
    for (int i = 0; i < m_SNOVA; ++i) {
        for (int j = 0; j < o_SNOVA; ++j) {
            for (int k = 0; k < o_SNOVA; ++k) {
                for (int index = 0; index < v_SNOVA; ++index) {
                    gf16m_mul(T12[index][j], F12[i][index][k], temp1);
                    gf16m_mul(P21[i][j][index], T12[index][k], temp2);
                    gf16m_add(temp1, temp2, temp1);
                    gf16m_add(P22[i][j][k], temp1, P22[i][j][k]);
                    // P22[i][j][k].printout();
                }
            }
        }
    }
    convert_GF16s_to_bytes(outP22, (uint8_t*)P22,
                           m_SNOVA * o_SNOVA * o_SNOVA * lsq_SNOVA);
}

/**
 * P22 byte to GF16
 * @param P22_gf16s - output
 * @param P22_bytes - input
 */
void input_P22(uint8_t* P22_gf16s, uint8_t* P22_bytes) {
    convert_bytes_to_GF16s(P22_bytes, P22_gf16s,
                           m_SNOVA * o_SNOVA * o_SNOVA * lsq_SNOVA);
}

/**
 * generate snova key elements.
 * @param key_elems - pointer to output snova key elements.
 * @param pk_seed - pointer to input public key seed.
 * @param sk_seed - pointer to input private key elements.
 */
void generate_keys_core(snova_key_elems* key_elems, uint8_t* pk_seed,
                        uint8_t* sk_seed) {
    gen_seeds_and_T12(key_elems->T12, sk_seed);

    memcpy(key_elems->pk.pt_public_key_seed, pk_seed, seed_length_public);

    gen_A_B_Q_P(&(key_elems->map1), pk_seed);

    gen_F(&(key_elems->map2), &(key_elems->map1), key_elems->T12);

    gen_P22(key_elems->T12, key_elems->map1.P21, key_elems->map2.F12,
            key_elems->pk.P22);
}

/**
 * Pack expanded private key. esk = (key_elems, pt_private_key_seed).
 * @param esk - pointer to output expanded private key.
 * @param key_elems - pointer to input snova key elements.
 * @param pt_private_key_seed - pointer to input private key seed.
 */
void sk_pack(uint8_t* esk, snova_key_elems* key_elems,
             uint8_t* pt_private_key_seed) {
    uint8_t* sk_gf16_ptr = (uint8_t*)(key_elems->map1.Aalpha);
    convert_GF16s_to_bytes_merger_in_half(
        esk, sk_gf16_ptr,
        (bytes_sk - (seed_length_public + seed_length_private)) * 2);
    memcpy(esk + (bytes_sk - (seed_length_public + seed_length_private)),
           key_elems->pk.pt_public_key_seed, seed_length_public);
    memcpy(esk + (bytes_sk - seed_length_private), pt_private_key_seed,
           seed_length_private);
}

/**
 * Unpack expanded secret key. skupk = (esk).
 * @param skupk - pointer to output private key (unpack).
 * @param esk - pointer to input expanded private key.
 */
void sk_unpack(sk_gf16* skupk, const uint8_t* esk) {
    convert_bytes_to_GF16s_cut_in_half(
        esk, (uint8_t*)skupk,
        (bytes_sk - (seed_length_public + seed_length_private)) * 2);
    memcpy(skupk->pt_public_key_seed,
           esk + (bytes_sk - (seed_length_public + seed_length_private)),
           seed_length_public + seed_length_private);
}

/**
 * Pack public key. pk = (key_elems).
 * @param pk - pointer to output public key.
 * @param key_elems - pointer to input snova key elements.
 */
void pk_pack(uint8_t* pk, snova_key_elems* key_elems) {
    memcpy(pk, key_elems->pk.pt_public_key_seed, bytes_pk);
}

/**
 * Generates public and private key. where private key is the seed of private
 * key.
 * @param pkseed - pointer to input public key seed.
 * @param skseed - pointer to input private key seed.
 * @param pk - pointer to output public key.
 * @param ssk - pointer to output private key.
 */
void generate_keys_ssk(const uint8_t* pkseed, const uint8_t* skseed,
                       uint8_t* pk, uint8_t* ssk) {
    snova_key_elems key_elems;
    generate_keys_core(&key_elems, pkseed, skseed);
    pk_pack(pk, &key_elems);
    memcpy(ssk, pkseed, seed_length_public);
    memcpy(ssk + seed_length_public, skseed, seed_length_private);
}

/**
 * Generates public and private key. where private key is the expanded version.
 * @param pkseed - pointer to input public key seed.
 * @param skseed - pointer to input private key seed.
 * @param pk - pointer to output public key.
 * @param esk - pointer to output private key. (expanded)
 */
void generate_keys_esk(const uint8_t* pkseed, const uint8_t* skseed,
                       uint8_t* pk, uint8_t* esk) {
    snova_key_elems key_elems;
    generate_keys_core(&key_elems, pkseed, skseed);
    pk_pack(pk, &key_elems);
    sk_pack(esk, &key_elems, skseed);
}

/**
 * Create salt
 * @param array_salt - pointer to output salt.
 */
void create_salt(uint8_t* array_salt) { randombytes(array_salt, bytes_salt); }

/**
 * Computes signature
 * @param pt_signature - pointer to output signature.
 * @param digest - pointer to input digest.
 * @param array_salt - pointer to input salt.
 * @param Aalpha -
 * @param Balpha -
 * @param Qalpha1 -
 * @param Qalpha2 -
 * @param T12 -
 * @param F11 -
 * @param F12 -
 * @param F21 -
 * @param pt_public_key_seed - pointer to output public key seed.
 * @param pt_private_key_seed - pointer to output private key seed.
 */
void sign_digest_core(uint8_t* pt_signature, uint8_t* digest,
                      uint64_t bytes_digest, uint8_t* array_salt,
                      Aalpha_t Aalpha, Balpha_t Balpha, Qalpha1_t Qalpha1,
                      Qalpha2_t Qalpha2, T12_t T12, F11_t F11, F12_t F12,
                      F21_t F21, uint8_t pt_public_key_seed[seed_length_public],
                      uint8_t pt_private_key_seed[seed_length_private]) {
    gf16_t Gauss[m_SNOVA * lsq_SNOVA][m_SNOVA * lsq_SNOVA + 1];
    gf16_t Temp[lsq_SNOVA][lsq_SNOVA];
    gf16_t t_GF16, solution[m_SNOVA * lsq_SNOVA];

    gf16m_t Left_X_tmp, Right_X_tmp;
    gf16_t *Left_X, *Right_X;

    gf16m_t Left[lsq_SNOVA][v_SNOVA], Right[lsq_SNOVA][v_SNOVA];
    gf16m_t X_in_GF16Matrix[n_SNOVA] = {0};
    gf16m_t Fvv_in_GF16Matrix[m_SNOVA];
    gf16_t hash_in_GF16[m_SNOVA * lsq_SNOVA];
    gf16m_t signature_in_GF16Matrix[n_SNOVA];

    uint8_t hash_input[seed_length_public + bytes_digest + bytes_salt];
    uint8_t signed_hash[bytes_hash];
    uint8_t vinegar_input[seed_length_private + bytes_digest + bytes_salt + 1];
    uint8_t vinegar_in_byte[(v_SNOVA * lsq_SNOVA + 1) >> 1];

    Left_X = Left_X_tmp;
    Right_X = Right_X_tmp;
    int flag_redo = 1;
    int num_sign = 0, counter;

    for (int index = 0; index < seed_length_public; ++index) {
        hash_input[index] = pt_public_key_seed[index];
    }
    for (int index = 0; index < bytes_digest; ++index) {
        hash_input[seed_length_public + index] = digest[index];
    }
    for (int index = 0; index < bytes_salt; ++index) {
        hash_input[(seed_length_public + bytes_digest) + index] =
            array_salt[index];
    }

    shake256(signed_hash, bytes_hash, hash_input, seed_length_public + bytes_digest + bytes_salt);

    // put hash value in GF16 array
    convert_bytes_to_GF16s(signed_hash, hash_in_GF16, GF16s_hash);

    for (int index = 0; index < seed_length_private; ++index) {
        vinegar_input[index] = pt_private_key_seed[index];
    }
    for (int index = 0; index < bytes_digest; ++index) {
        vinegar_input[seed_length_private + index] = digest[index];
    }
    for (int index = 0; index < bytes_salt; ++index) {
        vinegar_input[(seed_length_private + bytes_digest) + index] =
            array_salt[index];
    }

    // printf("hash value in GF16: \n");
    // print_gf16(key->hash_in_GF16, m_SNOVA * lsq_SNOVA);
    // printf("===================\n");
    //--------------------------------------
    do {
        num_sign++;
        flag_redo = 0;
        // put hash value in the last column of Gauss matrix
        for (int index = 0; index < (m_SNOVA * lsq_SNOVA); index++) {
            Gauss[index][m_SNOVA * lsq_SNOVA] = hash_in_GF16[index];
        }
        // generate the vinegar value
        vinegar_input[seed_length_private + bytes_digest + bytes_salt] =
            num_sign;

        shake256(vinegar_in_byte, (v_SNOVA * lsq_SNOVA + 1) >> 1, 
                 vinegar_input, seed_length_private + bytes_digest + bytes_salt + 1);
        counter = 0;

        for (int index = 0; index < v_SNOVA; index++) {
            for (int i = 0; i < rank; ++i) {
                for (int j = 0; j < rank; ++j) {
                    set_gf16m(X_in_GF16Matrix[index], i, j,
                              ((counter & 1)
                                   ? (vinegar_in_byte[counter >> 1] >> 4)
                                   : (vinegar_in_byte[counter >> 1] & 0xF)));
                    counter++;
                }
            }
        }

        // evaluate the vinegar part of central map
        for (int alpha = 0; alpha < lsq_SNOVA; ++alpha) {
            for (int index = 0; index < v_SNOVA; ++index) {
                gf16m_t X_in_GF16Matrix_P, temp;
                gf16m_transpose(X_in_GF16Matrix[index], X_in_GF16Matrix_P);
                gf16m_mul(X_in_GF16Matrix_P, Qalpha1[alpha], temp);
                gf16m_mul(Aalpha[alpha], temp, Left[alpha][index]);
                gf16m_mul(Qalpha2[alpha], X_in_GF16Matrix[index], temp);
                gf16m_mul(temp, Balpha[alpha], Right[alpha][index]);
                /*
                Left[alpha][index] = Aalpha[alpha] *
                (X_in_GF16Matrix[index].transpose()) * Qalpha1[alpha];
                Right[alpha][index] = Qalpha2[alpha] * X_in_GF16Matrix[index] *
                Balpha[alpha];
                */
            }
        }
        for (int i = 0; i < m_SNOVA; ++i) {
            gf16m_set_zero(Fvv_in_GF16Matrix[i]);
            for (int alpha = 0; alpha < lsq_SNOVA; ++alpha) {
                for (int j = 0; j < v_SNOVA; ++j) {
                    for (int k = 0; k < v_SNOVA; ++k) {
                        gf16m_t temp1, temp2;
                        gf16m_mul(Left[alpha][j], F11[i][j][k], temp1);
                        gf16m_mul(temp1, Right[alpha][k], temp2);
                        gf16m_add(Fvv_in_GF16Matrix[i], temp2,
                                  Fvv_in_GF16Matrix[i]);
                        /*
                        Fvv_in_GF16Matrix[i] = Fvv_in_GF16Matrix[i] +
                        Left[alpha][j] * F11[i][j][k] * Right[alpha][k];
                        */
                    }
                }
            }
        }
        // add to the last column of Gauss matrix
        for (int i = 0; i < m_SNOVA; ++i) {
            for (int j = 0; j < rank; ++j) {
                for (int k = 0; k < rank; ++k) {
                    int index1 = i * lsq_SNOVA + j * rank + k;
                    int index2 = m_SNOVA * lsq_SNOVA;
                    Gauss[index1][index2] =
                        gf16_get_add(Gauss[index1][index2],
                                     get_gf16m(Fvv_in_GF16Matrix[i], j, k));
                }
            }
        }

        // compute the coefficients of Xo and put into Gauss matrix and compute
        // the coefficients of Xo^t and add into Gauss matrix
        for (int i = 0; i < m_SNOVA; ++i) {
            for (int index = 0; index < o_SNOVA; ++index) {
                for (int ti = 0; ti < lsq_SNOVA; ++ti) {
                    for (int tj = 0; tj < lsq_SNOVA; ++tj) {
                        Temp[ti][tj] = 0;
                    }
                }

                for (int alpha = 0; alpha < lsq_SNOVA; ++alpha) {
                    for (int j = 0; j < v_SNOVA; ++j) {
                        gf16m_t temp1;
                        gf16m_mul(Left[alpha][j], F12[i][j][index], temp1);
                        gf16m_mul(temp1, Qalpha2[alpha], Left_X_tmp);
                        Left_X = Left_X_tmp;
                        Right_X = Balpha[alpha];
                        /*
                        Left_X = Left[alpha][j] * F12[i][j][k] * Qalpha2[alpha];
                        Right_X = Balpha[alpha];
                        */
                        for (int ti = 0; ti < lsq_SNOVA; ++ti) {
                            for (int tj = 0; tj < lsq_SNOVA; ++tj) {
                                gf16_t temp3 = 0;
                                temp3 = gf16_get_mul(
                                    get_gf16m(Left_X, ti / rank, tj / rank),
                                    get_gf16m(Right_X, tj % rank, ti % rank));
                                Temp[ti][tj] =
                                    gf16_get_add(Temp[ti][tj], temp3);
                            }
                        }
                    }
                }

                for (int alpha = 0; alpha < lsq_SNOVA; ++alpha) {
                    for (int j = 0; j < v_SNOVA; ++j) {
                        gf16m_t temp1;
                        Left_X = Aalpha[alpha];
                        gf16m_mul(Qalpha1[alpha], F21[i][index][j], temp1);
                        gf16m_mul(temp1, Right[alpha][j], Right_X_tmp);
                        Right_X = Right_X_tmp;
                        /*
                        Left_X = Aalpha[alpha];
                        Right_X = Qalpha1[alpha] * F21[i][k][j] *
                        Right[alpha][j];
                        */
                        for (int ti = 0; ti < lsq_SNOVA; ++ti) {
                            for (int tj = 0; tj < lsq_SNOVA; ++tj) {
                                gf16_t temp1 = 0;
                                temp1 = gf16_get_mul(
                                    get_gf16m(Left_X, ti / rank, tj % rank),
                                    get_gf16m(Right_X, tj / rank, ti % rank));
                                Temp[ti][tj] =
                                    gf16_get_add(Temp[ti][tj], temp1);
                            }
                        }
                    }
                }
                for (int ti = 0; ti < lsq_SNOVA; ++ti) {
                    for (int tj = 0; tj < lsq_SNOVA; ++tj) {
                        Gauss[i * lsq_SNOVA + ti][index * lsq_SNOVA + tj] =
                            Temp[ti][tj];
                    }
                }
            }
        }
        //
        // Gauss elimination
        for (int i = 0; i < m_SNOVA * lsq_SNOVA; ++i) {
            if (Gauss[i][i] == 0) {
                for (int j = i + 1; j < m_SNOVA * lsq_SNOVA; ++j) {
                    if (Gauss[j][i] != 0) {
                        for (int k = i; k < m_SNOVA * lsq_SNOVA + 1; ++k) {
                            t_GF16 = Gauss[i][k];
                            Gauss[i][k] = Gauss[j][k];
                            Gauss[j][k] = t_GF16;
                        }
                        break;
                    }
                }
            }
            if (Gauss[i][i] == 0) {
                flag_redo = 1;
                break;
            }

            t_GF16 = inv(Gauss[i][i]);
            for (int k = i; k < m_SNOVA * lsq_SNOVA + 1; ++k) {
                Gauss[i][k] = gf16_get_mul(Gauss[i][k], t_GF16);
            }

            for (int j = i + 1; j < m_SNOVA * lsq_SNOVA; ++j) {
                if (Gauss[j][i] != 0) {
                    t_GF16 = Gauss[j][i];
                    for (int k = i; k < m_SNOVA * lsq_SNOVA + 1; ++k) {
                        Gauss[j][k] = gf16_get_add(
                            Gauss[j][k], gf16_get_mul(Gauss[i][k], t_GF16));
                    }
                }
            }
        }

        if (!flag_redo) {
            for (int i = m_SNOVA * lsq_SNOVA - 1; i >= 0; --i) {
                t_GF16 = 0;
                for (int k = i + 1; k < m_SNOVA * lsq_SNOVA; ++k) {
                    t_GF16 = gf16_get_add(
                        t_GF16, gf16_get_mul(Gauss[i][k], solution[k]));
                }
                solution[i] =
                    gf16_get_add(Gauss[i][m_SNOVA * lsq_SNOVA], t_GF16);
            }
        }

    } while (flag_redo);
    // printf("times of Gauss elimination : %d\n", num_sign);
    for (int index = 0; index < o_SNOVA; ++index) {
        for (int i = 0; i < rank; ++i) {
            for (int j = 0; j < rank; ++j) {
                set_gf16m(X_in_GF16Matrix[index + v_SNOVA], i, j,
                          solution[index * lsq_SNOVA + i * rank + j]);
            }
        }
    }
    for (int index = 0; index < v_SNOVA; ++index) {
        gf16m_clone(signature_in_GF16Matrix[index], X_in_GF16Matrix[index]);
        for (int i = 0; i < o_SNOVA; ++i) {
            gf16m_t temp1;
            gf16m_mul(T12[index][i], X_in_GF16Matrix[v_SNOVA + i], temp1);
            gf16m_add(signature_in_GF16Matrix[index], temp1,
                      signature_in_GF16Matrix[index]);
            /*
            signature_in_GF16Matrix[index] = signature_in_GF16Matrix[index] +
            T12[index][i] * X_in_GF16Matrix[v_SNOVA + i];
            */
        }
    }
    for (int index = 0; index < o_SNOVA; ++index) {
        gf16m_clone(signature_in_GF16Matrix[v_SNOVA + index],
                    X_in_GF16Matrix[v_SNOVA + index]);
    }
    // output signature
    for (int index = 0; index < n_SNOVA * lsq_SNOVA; ++index) {
        ((gf16_t*)signature_in_GF16Matrix)[index] = get_gf16m(
            signature_in_GF16Matrix[index / lsq_SNOVA],
            (index % lsq_SNOVA) / l_SNOVA, (index % lsq_SNOVA) % l_SNOVA);
    }
    convert_GF16s_to_bytes(pt_signature, (gf16_t*)signature_in_GF16Matrix,
                           n_SNOVA * lsq_SNOVA);
    for (int i = 0; i < bytes_salt; ++i) {
        pt_signature[bytes_signature + i] = array_salt[i];
    }
}

/**
 * Compute the signature using ssk (private key seed). some preparatory work
 * before using sign_digest_core()
 * @param pt_signature - pointer to output signature.
 * @param digest - pointer to input digest.
 * @param array_salt - pointer to input salt.
 * @param ssk - pointer to input private key (seed).
 */
void sign_digest_ssk(uint8_t* pt_signature, const uint8_t* digest,
                     uint64_t bytes_digest, uint8_t* array_salt,
                     const uint8_t* ssk) {
    snova_key_elems key_elems;
    uint8_t* pk_seed = ssk;
    uint8_t* sk_seed = ssk + seed_length_public;

    gen_seeds_and_T12(key_elems.T12, sk_seed);
    gen_A_B_Q_P(&(key_elems.map1), pk_seed);
    gen_F(&(key_elems.map2), &(key_elems.map1), key_elems.T12);

    sign_digest_core(pt_signature, digest, bytes_digest, array_salt,
                     key_elems.map1.Aalpha, key_elems.map1.Balpha,
                     key_elems.map1.Qalpha1, key_elems.map1.Qalpha2,
                     key_elems.T12, key_elems.map2.F11, key_elems.map2.F12,
                     key_elems.map2.F21, pk_seed, sk_seed);
}

/**
 * Compute the signature using esk (). some preparatory work before using
 * sign_digest_core()
 * @param pt_signature - pointer to output signature.
 * @param digest - pointer to input digest.
 * @param array_salt - pointer to input salt.
 * @param esk - pointer to input private key (expanded).
 */
void sign_digest_esk(uint8_t* pt_signature, const uint8_t* digest,
                     uint64_t bytes_digest, uint8_t* array_salt,
                     const uint8_t* esk) {
    sk_gf16 sk_upk;
    sk_unpack(&sk_upk, esk);
    sign_digest_core(pt_signature, digest, bytes_digest, array_salt,
                     sk_upk.Aalpha, sk_upk.Balpha, sk_upk.Qalpha1,
                     sk_upk.Qalpha2, sk_upk.T12, sk_upk.F11, sk_upk.F12,
                     sk_upk.F21, sk_upk.pt_public_key_seed,
                     sk_upk.pt_private_key_seed);
}

/**
 * Verifies signature.
 * @param pt_digest - pointer to input digest.
 * @param pt_signature - pointer to output signature.
 * @param pk - pointer to output public key.
 * @returns - 0 if signature could be verified correctly and -1 otherwise
 */
int verify_signture(const uint8_t* pt_digest, uint64_t bytes_digest,
                    const uint8_t* pt_signature, const uint8_t* pk) {
    uint8_t hash_in_bytes[bytes_hash];
    uint8_t hash_input[seed_length_public + bytes_digest + bytes_salt];
    uint8_t signed_hash[bytes_hash];
    uint8_t* pt_salt = pt_signature + bytes_signature;

    gf16m_t Left[lsq_SNOVA][n_SNOVA], Right[lsq_SNOVA][n_SNOVA];
    gf16m_t hash_in_GF16Matrix[m_SNOVA];
    gf16m_t signature_in_GF16Matrix[n_SNOVA];
    gf16m_t P22[m_SNOVA][o_SNOVA][o_SNOVA];

    map_group1 map1;
    gf16m_t temp1, temp2;

    public_key* pk_stru = (public_key*)pk;

    for (int index = 0; index < seed_length_public; ++index) {
        hash_input[index] = pk_stru->pt_public_key_seed[index];
    }
    for (int index = 0; index < bytes_digest; ++index) {
        hash_input[seed_length_public + index] = pt_digest[index];
    }
    for (int index = 0; index < bytes_salt; ++index) {
        hash_input[(seed_length_public + bytes_digest) + index] =
            pt_salt[index];
    }

    shake256(signed_hash, bytes_hash, hash_input, seed_length_public + bytes_digest + bytes_salt);
#if (o_SNOVA * l_SNOVA) & 0x1 == 1
    signed_hash[bytes_hash - 1] &= 0x0f;
#endif

    convert_bytes_to_GF16s(pt_signature, (gf16_t*)signature_in_GF16Matrix,
                           GF16s_signature);
    // generate PRNG part of public key
    gen_A_B_Q_P(&map1, pk_stru->pt_public_key_seed);
    // read  P22
    input_P22((uint8_t*)P22, pk_stru->P22);
    // evaluate signature GF16Matrix array
    for (int alpha = 0; alpha < lsq_SNOVA; ++alpha) {
        for (int index = 0; index < n_SNOVA; ++index) {
            // Left[alpha][index]
            gf16m_transpose(signature_in_GF16Matrix[index], temp1);
            gf16m_mul(temp1, map1.Qalpha1[alpha], temp2);
            gf16m_mul(map1.Aalpha[alpha], temp2, Left[alpha][index]);
            // Right[alpha][index]
            gf16m_mul(map1.Qalpha2[alpha], signature_in_GF16Matrix[index],
                      temp2);
            gf16m_mul(temp2, map1.Balpha[alpha], Right[alpha][index]);

            /*
            Left[alpha][index] = Aalpha[alpha] *
            (signature_in_GF16Matrix[index].transpose()) * Qalpha1[alpha];
            Right[alpha][index] = Qalpha2[alpha] *
            signature_in_GF16Matrix[index] * Balpha[alpha];
            */
        }
    }

    for (int i = 0; i < m_SNOVA; ++i) {
        // gf16m_add(hash_in_GF16Matrix[i], hash_in_GF16Matrix[i],
        // hash_in_GF16Matrix[i]); //歸零
        gf16m_set_zero(hash_in_GF16Matrix[i]);
        for (int alpha = 0; alpha < lsq_SNOVA; ++alpha) {
            for (int dj = 0; dj < v_SNOVA; ++dj) {
                for (int dk = 0; dk < v_SNOVA; ++dk) {
                    gf16m_mul(Left[alpha][dj], map1.P11[i][dj][dk], temp1);
                    gf16m_mul(temp1, Right[alpha][dk], temp2);
                    gf16m_add(hash_in_GF16Matrix[i], temp2,
                              hash_in_GF16Matrix[i]);
                    // hash_in_GF16Matrix[i] = hash_in_GF16Matrix[i] +
                    // Left[alpha][dj] * P11[i][dj][dk] * Right[alpha][dk];
                }
            }
            for (int dj = 0; dj < v_SNOVA; ++dj) {
                for (int dk = 0; dk < o_SNOVA; ++dk) {
                    gf16m_mul(Left[alpha][dj], map1.P12[i][dj][dk], temp1);
                    gf16m_mul(temp1, Right[alpha][dk + v_SNOVA], temp2);
                    gf16m_add(hash_in_GF16Matrix[i], temp2,
                              hash_in_GF16Matrix[i]);
                    // hash_in_GF16Matrix[i] = hash_in_GF16Matrix[i] +
                    // Left[alpha][dj] * P12[i][dj][dk] *
                    // Right[alpha][dk+v_SNOVA];
                }
            }
            for (int dj = 0; dj < o_SNOVA; ++dj) {
                for (int dk = 0; dk < v_SNOVA; ++dk) {
                    gf16m_mul(Left[alpha][dj + v_SNOVA], map1.P21[i][dj][dk],
                              temp1);
                    gf16m_mul(temp1, Right[alpha][dk], temp2);
                    gf16m_add(hash_in_GF16Matrix[i], temp2,
                              hash_in_GF16Matrix[i]);
                    // hash_in_GF16Matrix[i] = hash_in_GF16Matrix[i] +
                    // Left[alpha][dj + v_SNOVA] * P21[i][dj][dk] *
                    // Right[alpha][dk];
                }
            }
            for (int dj = 0; dj < o_SNOVA; ++dj) {
                for (int dk = 0; dk < o_SNOVA; ++dk) {
                    gf16m_mul(Left[alpha][dj + v_SNOVA], P22[i][dj][dk], temp1);
                    gf16m_mul(temp1, Right[alpha][dk + v_SNOVA], temp2);
                    gf16m_add(hash_in_GF16Matrix[i], temp2,
                              hash_in_GF16Matrix[i]);
                    // hash_in_GF16Matrix[i] = hash_in_GF16Matrix[i] +
                    // Left[alpha][dj + v_SNOVA] * P22[i][dj][dk] *
                    // Right[alpha][dk + v_SNOVA];
                }
            }
        }
    }

    // GF16Matrix array ==>  GF16 array ==> byte array
    for (int index = 0; index < m_SNOVA * lsq_SNOVA; ++index) {
        ((gf16_t*)signature_in_GF16Matrix)[index] = get_gf16m(
            hash_in_GF16Matrix[index / lsq_SNOVA],
            (index % lsq_SNOVA) / l_SNOVA, (index % lsq_SNOVA) % l_SNOVA);
    }

    convert_GF16s_to_bytes(hash_in_bytes, (gf16_t*)signature_in_GF16Matrix,
                           m_SNOVA * lsq_SNOVA);
    int result = 0;
    for (int i = 0; i < bytes_hash; ++i) {
        if (hash_in_bytes[i] != signed_hash[i]) {
            result = -1;
            break;
        }
    }

    return result;
}
