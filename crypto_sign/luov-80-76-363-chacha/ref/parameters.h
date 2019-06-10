#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <stdint.h>

/*
We define the size of the finite field to use (options are 8, 16, 32, 48, 64 and 80 ),
OIL_VARS, the number oil variables and the number of equations in the UOV system,
VINEGAR_VARS, the number of vinegar variables,
SHAKENUM, the version of the shake XOF that is used (either 128 or 256)
and Wherther or not we are using Message Recovery or not
*/

// Enable / Disable message recovery mode
// #define MESSAGE_RECOVERY

// Choose a parameter set

/* SECURITY LEVEL 2 */
	/* SMALL PK
		#define FIELD_SIZE 48
		#define OIL_VARS 43
		#define VINEGAR_VARS 222
		#define SHAKENUM 128 */

	/* SMALL SIG
		#define FIELD_SIZE 8
		#define OIL_VARS 58
		#define VINEGAR_VARS 237
		#define SHAKENUM 128 */

/* SECURITY LEVEL 4 */
	/* SMALL PK
		#define FIELD_SIZE 64
		#define OIL_VARS 61
		#define VINEGAR_VARS 302
		#define SHAKENUM 256 */

	/* SMALL SIG
		#define FIELD_SIZE 8
		#define OIL_VARS 82
		#define VINEGAR_VARS 323
		#define SHAKENUM 256 */

/* SECURITY LEVEL 5 */
	/* SMALL PK */
		#define FIELD_SIZE 80
		#define OIL_VARS 76
		#define VINEGAR_VARS 363
		#define SHAKENUM 256

	/* SMALL SIG
		#define FIELD_SIZE 8
		#define OIL_VARS 107
		#define VINEGAR_VARS 371
		#define SHAKENUM 256  */

/* Custom parameter set */
 /*
	#define FIELD_SIZE
	#define OIL_VARS
	#define VINEGAR_VARS
	#define SHAKENUM
 */

#define SALT_BYTES 16
//#define PRNG_KECCAK
#define PRNG_CHACHA

/* derived parameters */
#define VARS (OIL_VARS+VINEGAR_VARS)

// The number of coefficients per polynomial of the Public map P
// that is not generated but stored in the public key
#define STORED_COLS_OF_P (OIL_VARS*(OIL_VARS + 1) / 2)

#ifdef MESSAGE_RECOVERY
	#define FIRST_PART_TARGET (SHAKENUM/4)
	#define SECOND_PART_TARGET ((OIL_VARS*sizeof(FELT))-FIRST_PART_TARGET)
	#define RECOVERED_PART_MESSAGE (SECOND_PART_TARGET-1)
#else
	#define FIRST_PART_TARGET OIL_VARS
	#define SECOND_PART_TARGET 0
	#define RECOVERED_PART_MESSAGE 0
#endif

#endif
