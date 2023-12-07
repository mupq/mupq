#ifndef ASSERTIONS_H_
#define ASSERTIONS_H_

#include <stdio.h>
#include <stdlib.h>

/**
 * @brief simple requirement: will kill the program if condition is not met
 */
#define REQUIRE_DRAMATICALLY(condition, ...)                                   \
  {                                                                            \
    if (!(condition)) {                                                        \
      fprintf(stderr, "Requirement failed: " __VA_ARGS__);                     \
      fprintf(stderr, "\n");                                                   \
      abort();                                                                 \
    }                                                                          \
  }

#define REQUIRE_ALIGNMENT(variabl, alignmt)                                    \
  {                                                                            \
    if ((uint64_t)(variabl) & ((alignmt)-1)) {                                 \
      fprintf(stderr, "Requirement failed: variable should be aligned\n");     \
      fprintf(stderr, "\n");                                                   \
      abort();                                                                 \
    }                                                                          \
  }

#ifndef NDEBUG
/**
 * @brief assertion executed in debug mode only: will kill the program if
 * condition is not met
 */
#define ASSERT_DRAMATICALLY(condition, ...)                                    \
  {                                                                            \
    if (!(condition)) {                                                        \
      fprintf(stderr, "Assertion failed: " __VA_ARGS__);                       \
      fprintf(stderr, "\n");                                                   \
      abort();                                                                 \
    }                                                                          \
  }
/**
 * @brief alignment assertion in debug mode only: will kill the program if
 * alignment is not met
 */
#define ASSERT_ALIGNMENT REQUIRE_ALIGNMENT
#else
/**
 * @brief assertion (not executed in production mode): would kill the program if
 * condition is not met
 */
#define ASSERT_DRAMATICALLY(condition, ...)                                    \
  {}
#define ASSERT_ALIGNMENT(variabl, alignmt)                                     \
  {}
#endif

#endif // ASSERTIONS_H_
