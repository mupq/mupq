// SPDX-License-Identifier: MIT

#include "api.h"
#include "portable_endian.h"
#include "aimer_instances.h"
#include "tree.h"


int
crypto_sign_keypair(unsigned char *pk, unsigned char *sk)
{
  int ret = 0;
  aimer_publickey_t public_key;
  aimer_privatekey_t secret_key;

  if (!pk || !sk)
  {
    return -1;
  }

  const aimer_instance_t* instance = aimer_instance_get(AIMER_INSTANCE);
  if (!instance)
  {
    return -1;
  }

  ret = aimer_keygen(instance->params, &public_key, &secret_key);
  if (ret == -1)
  {
    return ret;
  }
  
  size_t block_size = instance->aim_params.block_size;

  pk[0] = (uint8_t)instance->params;
  memcpy(pk + 1, public_key.data, 2 * block_size);

  sk[0] = (uint8_t)instance->params;
  memcpy(sk + 1, secret_key.data, block_size);
  memcpy(sk + 1 + block_size, public_key.data, 2 * block_size);

  return 0;
}

int
crypto_sign(unsigned char *sm, size_t *smlen,
            const unsigned char *m, size_t mlen,
            const unsigned char *sk)
{
  int ret = 0;
  size_t signature_len = 0;

  aimer_publickey_t public_key;
  aimer_privatekey_t secret_key;

  public_key.params = sk[0];

  const aimer_instance_t* instance = aimer_instance_get(public_key.params);
  if (!instance)
  {
    return -1;
  }
  size_t block_size = instance->aim_params.block_size;

  memset(&secret_key, 0x00, sizeof(aimer_privatekey_t));

  memcpy(secret_key.data, sk + 1, block_size);
  memcpy(public_key.data, sk + 1 + block_size, 2 * block_size);

  ret = aimer_sign(&public_key, &secret_key, m, mlen, sm + mlen, &signature_len);
  if (ret == -1)
  {
    return ret;
  }

  *smlen = mlen + signature_len;
  memcpy(sm, m, mlen);

  return 0;
}

int
crypto_sign_open(unsigned char *m, size_t *mlen,
                 const unsigned char *sm, size_t smlen,
                 const unsigned char *pk)
{
  int ret = 0;

  aimer_publickey_t public_key;
  memset(&public_key, 0, sizeof(aimer_publickey_t));

  public_key.params = pk[0];

  const aimer_instance_t* instance = aimer_instance_get(public_key.params);
  if (!instance)
  {
    return -1;
  }
  size_t block_size = instance->aim_params.block_size;

  size_t reveal_list_size = ceil_log2(instance->num_MPC_parties) * instance->seed_size;
  size_t signature_len = instance->salt_size + (2 * instance->digest_size) +
    instance->num_repetitions * (reveal_list_size + instance->digest_size +
    block_size + (instance->aim_params.num_input_sboxes + 2) * instance->field_size);

  if (signature_len > CRYPTO_BYTES)
  {
    return -1;
  }

  const size_t message_len = smlen - signature_len;
  const uint8_t* message = sm;
  const uint8_t* sig = sm + message_len;

  memcpy(public_key.data, pk + 1, 2 * block_size);

  ret = aimer_verify(&public_key, sig, signature_len, message, message_len);
  if (ret == -1)
  {
    return ret;
  }

  memcpy(m, message, message_len);
  *mlen = message_len;

  return 0;
}
