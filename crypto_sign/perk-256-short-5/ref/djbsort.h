
/**
 * @file djbsort.h
 * @brief Header file for sorting functions
 */

#ifndef DJB_SORT_H
#define DJB_SORT_H

#include <stdint.h>

extern void uint32_sort(uint32_t *, long long) __attribute__((visibility("default")));
extern void int32_sort(int32_t *, long long) __attribute__((visibility("default")));

#endif
