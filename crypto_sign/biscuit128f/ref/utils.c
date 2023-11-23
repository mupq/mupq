#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"

// #define shake256_inc_init(h_ctx)                                              \
//   Keccak_HashInitialize_SHAKE256 ((void *) h_ctx)
// #define shake256_inc_absorb(h_ctx, msg, len)                                  \
//   Keccak_HashUpdate ((void *) h_ctx, msg, len << 3)
// #define shake256_inc_finalize(h_ctx)                                          \
//   Keccak_HashFinal ((void *) h_ctx, NULL)
// #define shake256_inc_squeeze(out, len, h_ctx)                                 \
//   Keccak_HashSqueeze ((void *) h_ctx, out, len << 3)

void
commit (uint8_t *output, const uint8_t *salt, uint16_t e, uint16_t i,
        uint8_t *seed, int seclevel)
{
  h_ctx_t h_ctx;

  const uint8_t id[1] = { 0 };
  uint8_t extra[4];
  extra[0] = (uint8_t) e;
  extra[1] = (uint8_t) (e >> 8);
  extra[2] = (uint8_t) i;
  extra[3] = (uint8_t) (i >> 8);

  shake256_inc_init (h_ctx);
  shake256_inc_absorb (h_ctx, id, 1);
  shake256_inc_absorb (h_ctx, salt, seclevel >> 2);
  shake256_inc_absorb (h_ctx, extra, 4);
  shake256_inc_absorb (h_ctx, seed, seclevel >> 3);
  shake256_inc_finalize (h_ctx);
  shake256_inc_squeeze (output, seclevel >> 2, h_ctx);
}

void
expand_init (h_ctx_t h_ctx, const uint8_t *input, size_t inlen)
{
  shake256_inc_init (h_ctx);
  shake256_inc_absorb (h_ctx, input, inlen);
  shake256_inc_finalize (h_ctx);
}

void
expandtape_init (h_ctx_t h_ctx, const uint8_t *salt, uint16_t e, uint16_t i,
                 uint8_t *seed, int seclevel)
{
  const uint8_t id[1] = { 4 };
  uint8_t extra[4];
  extra[0] = (uint8_t) e;
  extra[1] = (uint8_t) (e >> 8);
  extra[2] = (uint8_t) i;
  extra[3] = (uint8_t) (i >> 8);

  shake256_inc_init (h_ctx);
  shake256_inc_absorb (h_ctx, id, 1);
  shake256_inc_absorb (h_ctx, salt, seclevel >> 2);
  shake256_inc_absorb (h_ctx, extra, 4);
  shake256_inc_absorb (h_ctx, seed, seclevel >> 3);
  shake256_inc_finalize (h_ctx);
}

void
expand_extract (uint8_t *output, size_t outlen, h_ctx_t h_ctx)
{
  shake256_inc_squeeze (output, outlen, h_ctx);
}

void
prf_init (h_ctx_t h_ctx, const uint8_t *input, size_t inlen)
{
  shake256_inc_init (h_ctx);
  if (input)
    shake256_inc_absorb (h_ctx, input, inlen);
}

void
prf_update (h_ctx_t h_ctx, const uint8_t *input, size_t inlen)
{
  shake256_inc_absorb (h_ctx, input, inlen);
}

void
prf_ready (h_ctx_t h_ctx)
{
  shake256_inc_finalize (h_ctx);
}

void
prf_generate (uint8_t *output, size_t outlen, h_ctx_t h_ctx)
{
  shake256_inc_squeeze (output, outlen, h_ctx);
}

void
H1_init (h_ctx_t h_ctx, const uint8_t *salt, const uint8_t *msg, size_t msglen,
         int seclevel)
{
  const uint8_t id[1] = { 1 };
  shake256_inc_init (h_ctx);
  shake256_inc_absorb (h_ctx, id, 1);
  shake256_inc_absorb (h_ctx, salt, seclevel >> 2);
  shake256_inc_absorb (h_ctx, msg, msglen);
}

void
H1_update (h_ctx_t h_ctx, const uint8_t *sigma1, size_t sigma1len)
{
  shake256_inc_absorb (h_ctx, sigma1, sigma1len);
}

void
H1_final (uint8_t *output, h_ctx_t h_ctx, int seclevel)
{
  shake256_inc_finalize (h_ctx);
  shake256_inc_squeeze (output, seclevel >> 2, h_ctx);
}

void
H2_init (h_ctx_t h_ctx, const uint8_t *salt, const uint8_t *h1, int seclevel)
{
  const uint8_t id[1] = { 2 };
  shake256_inc_init (h_ctx);
  shake256_inc_absorb (h_ctx, id, 1);
  shake256_inc_absorb (h_ctx, salt, seclevel >> 2);
  shake256_inc_absorb (h_ctx, h1, seclevel >> 2);
}

void
H2_update (h_ctx_t h_ctx, const uint8_t *sigma2, size_t sigma2len)
{
  shake256_inc_absorb (h_ctx, sigma2, sigma2len);
}

void
H2_final (uint8_t *output, h_ctx_t h_ctx, int seclevel)
{
  shake256_inc_finalize (h_ctx);
  shake256_inc_squeeze (output, seclevel >> 2, h_ctx);
}

int
ilog2 (int N)
{
  int lN = 0;
  N = N - 1;
  while (N)
    {
      lN++;
      N >>= 1;
    }
  return lN;
}

static void
child_seeds (uint8_t *output, const uint8_t *parent, const uint8_t *salt,
             uint16_t e, uint16_t layer, uint16_t index, int seclevel)
{
  h_ctx_t h_ctx;

  const uint8_t id[1] = { 3 };
  uint8_t extra[6];
  extra[0] = (uint8_t) e;
  extra[1] = (uint8_t) (e >> 8);
  extra[2] = (uint8_t) layer;
  extra[3] = (uint8_t) (layer >> 8);
  extra[4] = (uint8_t) index;
  extra[5] = (uint8_t) (index >> 8);

  shake256_inc_init (h_ctx);
  shake256_inc_absorb (h_ctx, id, 1);
  shake256_inc_absorb (h_ctx, salt, seclevel >> 2);
  shake256_inc_absorb (h_ctx, extra, 6);
  shake256_inc_absorb (h_ctx, parent, seclevel >> 3);
  shake256_inc_finalize (h_ctx);
  shake256_inc_squeeze (output, seclevel >> 2, h_ctx);
}

void
get_seeds (uint8_t *output, const uint8_t *seed, const uint8_t *salt,
           uint16_t e, uint16_t N, int seclevel)
{
  int l, lN;
  uint8_t buffer[64];

  lN = ilog2 (N);

  memcpy (output, seed, seclevel >> 3);
  for (l = 0; l < lN; l++)
    {
      int idx;
      for (idx = 0; idx < (1 << (l + 1)) && (idx << (lN - l - 1)) < N; idx++)
        {
          if ((idx & 1) == 0)
            {
              child_seeds (buffer, output, salt, e, l, idx >> 1, seclevel);
            }
          memcpy (output, buffer + (idx & 1) * (seclevel >> 3), seclevel >> 3);
          output += (seclevel >> 3) << (lN - l - 1);
        }
      output -= idx * (seclevel >> 3) << (lN - l - 1);
    }
}

void
get_path (uint8_t *output, const uint8_t *seed, const uint8_t *salt,
          uint16_t e, uint16_t ibar, uint16_t N, int seclevel)
{
  int l, lN;
  uint8_t buffer[64];
  uint8_t tmp[32];

  l = lN = ilog2 (N);

  memcpy (tmp, seed, seclevel >> 3);
  while (l--)
    {
      const int idx = ibar >> l;
      child_seeds (buffer, tmp, salt, e, lN - l - 1, idx >> 1, seclevel);
      memcpy (tmp, buffer + (idx & 1) * (seclevel >> 3),
              seclevel >> 3);
      memcpy (output, buffer + ((idx & 1) ^ 1) * (seclevel >> 3),
              seclevel >> 3);
      output += seclevel >> 3;
    }
}

void
get_path_seed (uint8_t *output, const uint8_t *path, const uint8_t *salt,
               uint16_t e, uint16_t i, uint16_t ibar, uint16_t N, int seclevel)
{
  int l, lN;
  uint8_t buffer[64];

  l = lN = ilog2 (N);

  while (l--)
    {
      if ((i >> l) != (ibar >> l))
        break;
      path += seclevel >> 3;
    }

  if (l < 0)
    return;

  memcpy (output, path, seclevel >> 3);
  while (l--)
    {
      const int idx = i >> l;
      child_seeds (buffer, output, salt, e, lN - l - 1, idx >> 1, seclevel);
      memcpy (output, buffer + (idx & 1) * (seclevel >> 3),
              seclevel >> 3);
    }
}

void
get_path_seeds (uint8_t *output, const uint8_t *path, const uint8_t *salt,
                uint16_t e, uint16_t ibar, uint16_t N, int seclevel)
{
  int l, lN;
  uint8_t buffer[64];

  lN = ilog2 (N);

  for (l = 0; l < lN; l++)
    {
      int idx, idx_bar;
      idx_bar = ibar >> (lN - l - 1);
      for (idx = 0; idx < (1 << (l + 1)) && (idx << (lN - l - 1)) < N; idx++)
        {
          if ((idx & 1) == 0)
            {
              if ((idx >> 1) == (idx_bar >> 1))
                {
                  memcpy (buffer + ((idx_bar & 1) ^ 1) * (seclevel >> 3), path,
                          seclevel >> 3);
                  path += seclevel >> 3;
                }
              else
                {
                  child_seeds (buffer, output, salt, e, l, idx >> 1, seclevel);
                }
            }
          memcpy (output, buffer + (idx & 1) * (seclevel >> 3), seclevel >> 3);
          output += (seclevel >> 3) << (lN - l - 1);
        }
      output -= idx * (seclevel >> 3) << (lN - l - 1);
    }
}
