#include <string.h>

#include "fips202_wrapper.h"

#ifndef MEDS_t
#include "params.h"
#endif

#include "seed.h"
#include "log.h"

void t_hash(uint8_t *stree, uint8_t *salt, int h, int i)
{
  keccak_state shake;

  int start = i;
  int end = i+1;

  uint8_t buf[MEDS_st_salt_bytes + MEDS_st_seed_bytes + sizeof(uint32_t)] = {0};

  memcpy(buf, salt, MEDS_st_salt_bytes);

  uint8_t *pos = buf + MEDS_st_salt_bytes + MEDS_st_seed_bytes;

  for (h = h+1; h < MEDS_seed_tree_height+1; h++)
  {
    start = start<<1;
    end = end<<1;

    if ((start << (MEDS_seed_tree_height - h)) >= MEDS_t)
      break;

    for (int i = start; i < end; i+=2)
    {
      if ((i << (MEDS_seed_tree_height - h)) >= MEDS_t)
        break;

      for (int j = 0; j < 4; j++)
        pos[j] = (SEED_TREE_ADDR(h-1, i>>1) >> (j*8)) & 0xff;

      memcpy(buf + MEDS_st_salt_bytes,
          &stree[MEDS_st_seed_bytes * SEED_TREE_ADDR(h-1, i>>1)],
          MEDS_st_seed_bytes);

      shake256_init(&shake);
      shake256_absorb(&shake, buf, MEDS_st_salt_bytes + MEDS_st_seed_bytes + sizeof(uint32_t));
      shake256_finalize(&shake);

      int len = 2*MEDS_st_seed_bytes;

      if (((i+1) << (MEDS_seed_tree_height - h)) >= MEDS_t)
        len = MEDS_st_seed_bytes;

      shake256_squeeze(&stree[MEDS_st_seed_bytes * SEED_TREE_ADDR(h, i)], len, &shake);
    }
  }
}

void stree_to_path_to_stree(uint8_t *stree, uint8_t *h_digest, uint8_t *path, uint8_t *salt, int mode)
{
  int h = 0;
  int i = 0;

  unsigned int id = 0;

  int idx = 0;

  unsigned int indices[MEDS_w] = {0};

  for (int i = 0; i < MEDS_t; i++)
    if (h_digest[i] > 0)
      indices[idx++] = i;


  while (1 == 1)
  {
    int index_leaf = 0;

    while ((indices[id] >> (MEDS_seed_tree_height - h)) == i)
    {
      h += 1;
      i <<= 1;

      if (h > MEDS_seed_tree_height)
      {
        h -= 1;
        i >>= 1;

        if (id+1 < MEDS_w)
          id++;

        index_leaf = 1;

        break;
      }
    }

    if (index_leaf == 0)
    {
      if (mode == STREE_TO_PATH)
      {
        memcpy(path, &stree[MEDS_st_seed_bytes * SEED_TREE_ADDR(h, i)], MEDS_st_seed_bytes);
        path += MEDS_st_seed_bytes;
      }
      else
      {
        memcpy(&stree[MEDS_st_seed_bytes * SEED_TREE_ADDR(h, i)], path, MEDS_st_seed_bytes);
        path += MEDS_st_seed_bytes;

        t_hash(stree, salt, h, i);
      }
    }

    // backtrack
    while ((i & 1) == 1)
    {
      h -= 1;
      i >>= 1;
    }

    // next leaf
    i+=1;

    if ((i << (MEDS_seed_tree_height - h)) >= MEDS_t)
      return;
  }
}

