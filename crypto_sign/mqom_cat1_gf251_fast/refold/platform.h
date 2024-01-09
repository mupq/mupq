#ifndef MQOM_ALL_PARAMETERS_H
#define MQOM_ALL_PARAMETERS_H

/* on Mac OS-X and possibly others, ALIGN(x) is defined in param.h, and -Werror chokes on the redef. */
#ifdef ALIGN
#undef ALIGN
#endif

#if defined(__GNUC__)
#define ALIGN(x) __attribute__ ((aligned(x)))
#elif defined(_MSC_VER)
#define ALIGN(x) __declspec(align(x))
#elif defined(__ARMCC_VERSION)
#define ALIGN(x) __align(x)
#else
#define ALIGN(x)
#endif

#define ALIGNED_ALLOC(alignment,size) aligned_alloc(alignment,size)

#endif /* MQOM_ALL_PARAMETERS_H */
