#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <stdint.h>

/*
We define the size of the finite field to use (options are 8, 16, 32, 48, 64 and 80 ),
OIL_VARS, the number oil variables and the number of equations in the UOV system,
VINEGAR_VARS, the number of vinegar variables,
SHAKEVENUM, the version of the shake XOF that is used (either 128 or 256)
and Wherther or not we are using Message Recovery or not
*/

// Enable / Disable message recovery mode
// #define MESSAGE_RECOVERY

// Choose a parameter set

/* SECURITY LEVEL 2 */
	/* SMALL PK
		#define FIELD_SIZE 47
		#define OIL_VARS 42
		#define VINEGAR_VARS 182
		#define SHAKENUM 128
		#define FIRST_PART_TARGET 32 */

	/* SMALL SIG */
		#define FIELD_SIZE 7
		#define OIL_VARS 57
		#define VINEGAR_VARS 197
		#define SHAKENUM 128
		#define FIRST_PART_TARGET 32

/* SECURITY LEVEL 4 */
	/* SMALL PK
		#define FIELD_SIZE 61
		#define OIL_VARS 60
		#define VINEGAR_VARS 261
		#define SHAKENUM 256
		#define FIRST_PART_TARGET 48 */

	/* SMALL SIG
		#define FIELD_SIZE 7
		#define OIL_VARS 83
		#define VINEGAR_VARS 283
		#define SHAKENUM 256
		#define FIRST_PART_TARGET 48 */

/* SECURITY LEVEL 5 */
	/* SMALL PK
		#define FIELD_SIZE 79
		#define OIL_VARS 76
		#define VINEGAR_VARS 341
		#define SHAKENUM 256
		#define FIRST_PART_TARGET 64 */

	/* SMALL SIG
		#define FIELD_SIZE 7
		#define OIL_VARS 110
		#define VINEGAR_VARS 374
		#define SHAKENUM 256
		#define FIRST_PART_TARGET 64 */

/* Custom parameter set */
 /*
	#define FIELD_SIZE
	#define OIL_VARS
	#define VINEGAR_VARS
	#define SHAKENUM
 */

#define SALT_BYTES 16
#define PRNG_KECCAK
//#define PRNG_CHACHA

/* derived parameters */
#define VARS (OIL_VARS+VINEGAR_VARS)
#define FIELD_MASK ((((uint64_t) 1) << FIELD_SIZE) -1)

// The number of coefficients per polynomial of the Public map P
// that is not generated but stored in the public key
#define STORED_COLS_OF_P (OIL_VARS*(OIL_VARS + 1) / 2)

// Defines the field operations
#define __DEFINE_OPERATION(FS,OPERATION) f##FS##OPERATION
#define _DEFINE_OPERATION(FS,OPERATION) __DEFINE_OPERATION(FS,OPERATION)
#define DEFINE_OPERATION(OPERATION) _DEFINE_OPERATION(FIELD_SIZE,OPERATION)

#define FELT DEFINE_OPERATION(FELT)
#define serialize_FELT DEFINE_OPERATION(serialize_FELT)
#define deserialize_FELT DEFINE_OPERATION(deserialize_FELT)
#define printFELT DEFINE_OPERATION(printFELT)
#define isEqual DEFINE_OPERATION(isEqual)
#define multiply DEFINE_OPERATION(multiply)
#define minus(x) (x)
#define add DEFINE_OPERATION(add)
#define subtract(x,y) add(x,y)
#define inverse DEFINE_OPERATION(inverse)
#define ZERO DEFINE_OPERATION(ZERO)
#define ONE DEFINE_OPERATION(ONE)
#define xpown DEFINE_OPERATION(xpown)
#define scalarMultiply DEFINE_OPERATION(scalarMultiply)
#define addInPlace DEFINE_OPERATION(addInPlace)

#endif

