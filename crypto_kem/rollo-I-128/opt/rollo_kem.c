#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "typedef.h"
#include "matrix.h"
#include "gf2mn.h"
#include "rng_gf2mn.h"
#include "api.h"
#include "rollo_kem.h"
#include "utils.h"
#include "sha2.h"

/**
 * \fn static u1 GenKey( _SecretKey *SecretKey, _PublicKey *PublicKey)
 * \brief This function generates the private and public keys
 *
 * \param[out]  *SecretKey : address of the private key's structure.
 * \param[out]  *PublicKey : address of the public key's structure.
 */

static u1 GenKey( _SecretKey *SecretKey, _PublicKey *PublicKey)
{

  // Elements place in memory
  pu4 pu4Sk_a_copy = SecretKey->pu4Sk_a + NM_WORDS_PAD;
  pu4 pu4a_inv = pu4Sk_a_copy + NM_WORDS_PAD;
  pu4 pu4modulo = pu4Sk_a_copy + NM_WORDS_PAD;

  // Step 1 : generation of the private key support of rank PARAM_D
  Create_support_using_rng(SecretKey->pu4F, PARAM_D);

  // Step 2 : First key part generation from the support

  Create_gf2mn_using_support(SecretKey->pu4Sk_a,   SecretKey->pu4F,  PARAM_D);

  //Step 3 : Copy of the first part of the key and computation of its inverse

  memcpy(pu4Sk_a_copy, SecretKey->pu4Sk_a, NM_WORDS*4);

  init_modulo(pu4modulo, PARAM_N);
  modular_inverse_GF2mn(pu4a_inv, pu4Sk_a_copy,pu4modulo);
  memset(pu4Sk_a_copy,0x00,NM_WORDS*4);


  // Step 5 : Creation of the second part of the private key from the support pu4F

  Create_gf2mn_using_support(SecretKey->pu4Sk_b, SecretKey->pu4F, PARAM_D);


  // Step 6 : Computation of the public key
  u1 pad_Karatsuba = Karatsuba_Degree(PARAM_N);
  pu4 pu4R_Karatsuba = pu4a_inv + NM_WORDS_PAD;
  pu4 pu4Workspace_Karatsuba =  pu4R_Karatsuba + 2*NM_WORDS_PAD;

  Karatsuba(pu4a_inv,
            SecretKey->pu4Sk_b,
            pu4R_Karatsuba,
            pu4Workspace_Karatsuba,
            pad_Karatsuba,
            PARAM_M_LEN_WORD);

  memset(pu4a_inv,0x00,NM_WORDS*4);

  ModularReduction(pu4R_Karatsuba,
                   2*PARAM_N,
                   PARAM_N,
                   PARAM_M_LEN_WORD);

  //The workspace is cleared
  memcpy(PublicKey->pu4Pk, pu4R_Karatsuba, NM_WORDS*4);
  memset(pu4R_Karatsuba,0x00,4*NM_WORDS*4);

  return 0;
}



/**
 * \fn static pu4 Encaps(_PublicKey PublicKey, pu4 pu4Cipher, pu1 ss, pu4 pu4Workspace)
 * \brief This function encapsulates a secret data and thus generates a cipher = pk*e_1 + e_2  and the shared secret.
 *
 * \param[in]  PublicKey : address of the private key's structure.
 * \param[out]  pu4Cipher : address of the output cipher.
 * \param[out]        ss  : address of the shared secret.
 * \param[in]   pu4Workspace : workspace.
 */

static pu4 Encaps(_PublicKey PublicKey, pu4 pu4Cipher, pu1 ss, pu4 pu4Workspace)
{

  // Step 1 : Generation of the error's support of rank PARAM_R
  pu4 pu4Support = pu4Workspace;
  Create_support_using_rng(pu4Support, PARAM_R);

  // Step 2 : Generation of error e_1 from pu4Support
  pu4 pu4Error      = pu4Support + NM_WORDS;
  Create_gf2mn_using_support(pu4Error, pu4Support, PARAM_R);

  // Step 3 : Multiplication between PublicKey.pu4Pk and pu4Error (pk*e_1)
  u1 pad_Karatsuba = Karatsuba_Degree(PARAM_N);
  pu4 pu4BufferCalculus = pu4Error + NM_WORDS_PAD*PARAM_M_LEN_WORD;
  Karatsuba(pu4Error,
            PublicKey.pu4Pk,
            pu4BufferCalculus,
            pu4BufferCalculus + 2 * NM_WORDS_PAD,
            pad_Karatsuba,
            PARAM_M_LEN_WORD);

  memset(pu4BufferCalculus + 2 * NM_WORDS_PAD,
         0x00,
         NM_WORDS_PAD*4);

  ModularReduction(pu4BufferCalculus,
                   2*PARAM_N,
                   PARAM_N,
                   PARAM_M_LEN_WORD);

 // Step 4 : Generation of the second part of the error (e_2) from pu4Support, it replaces the first part of the error to gain memory usage

  memset(pu4Error, 0x00, NM_WORDS*4);

  Create_gf2mn_using_support(pu4Error, pu4Support, PARAM_R);


  // Step 5 : Addition between the previous multiplication result pk*e_1 and e_2
  Add_GF2mn(PARAM_N, PARAM_M_LEN_WORD, pu4Error, pu4BufferCalculus, pu4Cipher);

  //Step 6 : Computation of the shared secret
  Reduced_Echelon_matrix(pu4Support,PARAM_R,PARAM_M,PARAM_M_LEN_WORD,pu4Support);
  pu1 support = (pu1) (pu4Support + PARAM_R*PARAM_M_LEN_WORD);
  word_to_byte(support, pu4Support, PARAM_M_LEN_WORD,PARAM_R);

  sha512(ss ,support, PARAM_R*PARAM_M_LEN_BYTES);

  // Clear of the workspace
  memset(pu4Workspace, 0x00, (3*NM_WORDS +PARAM_R*PARAM_M_LEN_WORD)*4);


  return 0;
}



/**
 * \fn static u2 Decaps(_SecretKey SecretKey, _PublicKey PublicKey, pu4 pu4Cipher, pu1 ss, pu4 pu4Workspace)
 * \brief This function decapsalutes a secret data from the cipher generated during the encapsulation and computes the shared secret.
 *
 * \param[in]  SecretKey : address of the private key's structure.
 * \param[in]  pu4Cipher : address of the input cipher.
 * \param[out]        ss  : address of the shared secret.
 * \param[in]   pu4Workspace : workspace.
 */


static u2 Decaps(_SecretKey SecretKey, pu4 pu4Cipher, pu1 ss, pu4 pu4Workspace)
{
  pu4 pu4Syndrome            = pu4Workspace;
  u1 pad_Karatsuba = Karatsuba_Degree(PARAM_N);
  pu4 pu4Workspace_tmp = pu4Syndrome + 2* NM_WORDS_PAD;


  // Step 1 : Syndrome computation
  Karatsuba(SecretKey.pu4Sk_a,
            pu4Cipher,
            pu4Syndrome,
            pu4Workspace_tmp,
            pad_Karatsuba,
            PARAM_M_LEN_WORD);
  memset(pu4Workspace_tmp,
         0x00,
         3*NM_WORDS_PAD*4);
  ModularReduction(pu4Syndrome,
                   2*pad_Karatsuba,
                   PARAM_N,
                   PARAM_M_LEN_WORD);


  // Step 2 : Recovery of the support of the error
  pu4Workspace_tmp= pu4Syndrome + NM_WORDS_PAD;
  pu4 pu4Result = pu4Workspace_tmp + 2*PARAM_M_LEN_WORD*PARAM_D*PARAM_R;
  if(Rank_Support_Recover(pu4Result, SecretKey.pu4F, pu4Syndrome, pu4Workspace_tmp))
  {
    return 1;
  }


  // Step 3 : Computation of the shared secret from the recovered support

  pu1 Recover_support = (pu1)(pu4Result + PARAM_R*PARAM_M_LEN_WORD);
  word_to_byte(Recover_support, pu4Result, PARAM_M_LEN_WORD,PARAM_R);
  sha512(ss , Recover_support, PARAM_R*PARAM_M_LEN_BYTES);

  //Clear the workspace
  memset(pu4Workspace, 0x00, (2*PARAM_M_LEN_WORD*PARAM_D+4*PARAM_M_LEN_WORD*PARAM_R+NM_WORDS)*4);

  return 0;
}

u1 CryptoRAM[0x2000];

int crypto_kem_keypair(unsigned char* pk, unsigned char* sk)
{
  _SecretKey SecretKey;
  _PublicKey PublicKey;

  memset(CryptoRAM,0x00,0x1000);
  SecretKey.pu4F = (pu4)CryptoRAM;
  SecretKey.pu4Sk_a = SecretKey.pu4F + PARAM_D  * PARAM_M_LEN_WORD;
  SecretKey.pu4Sk_b = SecretKey.pu4Sk_a + NM_WORDS_PAD;
  PublicKey.pu4Pk = SecretKey.pu4Sk_b + NM_WORDS_PAD;

  GenKey(&SecretKey, &PublicKey);

  memcpy(sk, &SecretKey, sizeof(SecretKey));
  memcpy(pk, &PublicKey, sizeof(PublicKey));

  return 0;
}

int crypto_kem_enc(unsigned char* ct, unsigned char* ss, unsigned char* pk)
{
  u4 workspace[(3*NM_WORDS +PARAM_R*PARAM_M_LEN_WORD)*4];

  Encaps(*(_PublicKey *)pk, (pu4)ct, ss, workspace);

  return 0;
}

int crypto_kem_dec(unsigned char* ss, unsigned char* ct, unsigned char* sk)
{
  u4 workspace[(2*PARAM_M_LEN_WORD*PARAM_D+4*PARAM_M_LEN_WORD*PARAM_R+NM_WORDS)*4];

  return Decaps(*(_SecretKey *)sk, (pu4)ct, ss, workspace);
}
