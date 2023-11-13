#if !defined( ASCON_OFFSETS_H_ )
#define ASCON_OFFSETS_H_

/*
 * Offsets of various fields in the address structure when we use ASCON_HASH as
 * the ASCON_SIGN hash function
 */

#define ASGN_OFFSET_LAYER     3   /* The byte used to specify the Merkle tree layer */
#define ASGN_OFFSET_TREE      8   /* The start of the 8 byte field used to specify the tree */
#define ASGN_OFFSET_TYPE      19  /* The byte used to specify the hash type (reason) */
#define ASGN_OFFSET_KP_ADDR2  22  /* The high byte used to specify the key pair (which one-time signature) */
#define ASGN_OFFSET_KP_ADDR1  23  /* The low byte used to specify the key pair */
#define ASGN_OFFSET_CHAIN_ADDR 27  /* The byte used to specify the chain address (which Winternitz chain) */
#define ASGN_OFFSET_HASH_ADDR 31  /* The byte used to specify the hash address (where in the Winternitz chain) */
#define ASGN_OFFSET_TREE_HGT  27  /* The byte used to specify the height of this node in the FORS or Merkle tree */
#define ASGN_OFFSET_TREE_INDEX 28 /* The start of the 4 byte field used to specify the node in the FORS or Merkle tree */

#endif /* ASCON_OFFSETS_H_ */
