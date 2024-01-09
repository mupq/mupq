#ifndef MQOM_GF251_INTERNAL_H
#define MQOM_GF251_INTERNAL_H

#define PRIME 251
// Work as soon as x <= 128510 = 2^16.97
#define _gf251_reduce16(x) ((x) - 251 * (((uint32_t)(x) * 33421) >> 23))
// Work as soon as x <= 7967475679 = 2^32.89
#define _gf251_reduce32(x) ((x) - 251 * (((uint64_t)(x) * 2190262207) >> 39))

#endif /* MQOM_GF251_INTERNAL_H */
