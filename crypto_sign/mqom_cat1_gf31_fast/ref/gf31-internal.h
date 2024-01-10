#ifndef MQOM_GF31_INTERNAL_H
#define MQOM_GF31_INTERNAL_H

#define PRIME 31
// Work as soon as x <= 63487 = 2^15.95
#define _gf31_reduce16(x) ((x) - 31 * (((uint32_t)(x) * 67651) >> 21))
// Work as soon as x <= 4160749567 = 2^31.95
#define _gf31_reduce32(x) ((x) - 31 * (((uint64_t)(x) * 4433514629) >> 37))

#endif /* MQOM_GF31_INTERNAL_H */
