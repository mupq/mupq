#include "LUOV.h"
#include "randombytes.h"

/*
	Calculates Q_2, the last OIL_VARS*(OIL_VARS+1)/2 columns of the macaulay matrix of the public system and writes it to pk

	T : the secret linear map
	pk : the public key
*/
static void calculateQ2(column *T , unsigned char *pk) {
	int i, j, k;
	column *TempMat = malloc(sizeof(column) * OIL_VARS);
	column r;

	// Absorb the public seed in a sponge object
	ColumnGenerator CG;
	ColumnGenerator_init(&CG, PK_SEED(pk));

	// Allocate memory for temporary matrices that will store the values P_i,1 T + P_i,2 for i from 1 to OIL_VARS.
	// These OIL_VARS matrices are bitsliced into one OIL_VARS by OIL_VARS array of columns.
	// All entries of the matrices are initialized to zero.

	column **TempMat2 = malloc(sizeof(column*) * OIL_VARS);
	for (i = 0; i < OIL_VARS; i++) {
		TempMat2[i] = malloc(sizeof(column)*OIL_VARS);
		for (j = 0; j < OIL_VARS; j++) {
			TempMat2[i][j] = empty;
		}
	}

	// Simultaneously calculate P_i,1*T + P_i,2 for all i from 1 to OIL_VARS
	for (i = 0; i <= VINEGAR_VARS; i++) {
		for (j = 0; j < OIL_VARS; j++) {
			TempMat[j] = empty;
		}

		// Calculates P_i,1*T
		for (j = i; j <= VINEGAR_VARS; j++) {
			r = Next_Column(&CG);
			for (k = 0; k < OIL_VARS; k++) {
				if (getBit(T[j], k)) {
					TempMat[k] = xor(TempMat[k], r);
				}
			}
		}
		// Add P_i,2
		for (j = 0; j < OIL_VARS; j++) {
			r = Next_Column(&CG);
			TempMat[j] = xor(TempMat[j], r);
		}


		// Calculate P_i,3 = Transpose(T)*TempMat_i, and store the result in Q_2
		for (k = 0; k < OIL_VARS; k++) {
			for (j = 0; j < OIL_VARS; j++) {
				if (getBit(T[i] , k ) ) {
					TempMat2[k][j] = xor(TempMat2[k][j], TempMat[j]);
				}
			}
		}
	}

	// Write Q2 to the public key
	writer W = newWriter(PK_Q2(pk));
	for (i = 0; i < OIL_VARS; i++) {
		for (j = i; j < OIL_VARS; j++) {
			column col = TempMat2[i][j];
			if (i != j)
				col = xor(col,TempMat2[j][i]);
			serialize_column(&W, col);
		}
	}

	// Fill the remainder of the last byte of pk with 0's
	while (W.bitsUsed != 0)
		writeBit(&W,0);

	// Free the memory occupied by TempMat & TempMat2
	for (i = 0; i < OIL_VARS; i++) {
		free(TempMat2[i]);
	}
	free(TempMat2);
	free(TempMat);
}

/*
	Generates a key pair

	pk : receives the public key
	sk : receives the secret key
*/
int luov_keygen(unsigned char *pk, unsigned char *sk) {
	//printIntermediateValue("--- Start keygen ---\n");

	// Pick random secret key
	randombytes(sk , 32);

	// Calculate public seed
	Sponge sponge;
	initializeAndAbsorb(&sponge , sk, 32);
	squeezeBytes(&sponge, PK_SEED(pk) , 32);

	// Calculate T
	column T[VINEGAR_VARS+1];
	T[0]=empty;/* makes T linear instead of affine*/
	squeeze_column_array(&sponge , &(T[1]) , VINEGAR_VARS);

	// Calculates Q_2, the part of the public map P that cannot be generated from the public seed
	calculateQ2( T , pk );

	//printIntermediateValue("--- End keygen ---\n");
	return 0;
}


/*
	Builds the augmented matrix for the system F(x) = target , after fixing the vinegar variables

	A                 : Receives the augmented matrix, should be initialized to the zero matrix
	vinegar_variables : An assignment to the vinegar variables
	target            : The target vector to find a solution for
	T                 : The V-by-M matrix that determines the secret linear transformation T
	publicseed            : The public seed that is used to generate the first part of the secret key
*/
static void BuildAugmentedMatrix(Matrix A, const FELT *vinegar_variables , const FELT *target, const column *T, const unsigned char *publicseed){ //Sponge *sponge) {
	int i, j, k, x;
	column **F2;
	uint16_t r;
	FELT prod;

	// Sets the right hand side of the Augmented matrix to the target vector
	for (k = 0; k < OIL_VARS; k++) {
		A.array[k][OIL_VARS] = target[k];
	}

	// Allocate memory for the matrices F_1,2, ... , F_OIL_VARS,2.
	// These matrices are bit sliced together and stored as a OIL_VARS by VINEGAR_VARS array of columns.
	// All entries of these matrices are initialized to zero.
	F2 = malloc(sizeof(column*)*16);
	for (k = 0; k < 16; k++) {
		F2[k] = calloc(VINEGAR_VARS+1,sizeof(column));
	}

	// loop over each set of 16 equations
	for(x=0 ; x< 1+ (OIL_VARS-1)/16; x++){

		// Initialize the PRNG to produce the first part of 16 components of the public map
		PieceGenerator PG;
		PieceGenerator_init(&PG,publicseed,x);

		// Clear F2
		for (i = 0; i < 16; i++) {
			for(j=0; j<=VINEGAR_VARS; j++){
				F2[i][j] = empty;
			}
		}

		// Computes F_i,2 = (P_i,1 + Transpose(P_i,1)*T + P_i,2 ) simultaneously for all i from 1 to OIL_VARS
		// and subtracts the evaluation of P in the vinegar variables from the right hand side
		for (i = 0; i <= VINEGAR_VARS; i++) {
			for (j = i; j <= VINEGAR_VARS; j++) {
				r = Next_Piece(&PG);

				prod = multiply(vinegar_variables[i], vinegar_variables[j]);
				for (k = 0; k < 16; k++) {
					if (x*16 + k >= OIL_VARS)
						break;
					if ( r & ((uint16_t) 1 ) << k ) {
						// subtract the term in v_i*v_j from the right hand side
						A.array[x*16+k][OIL_VARS] = add(A.array[x*16+k][OIL_VARS],prod);

						// add (P_i,1 + Transpose(P_i,1)*T part to F_i,2
						F2[k][j] = xor(F2[k][j],T[i]);
						F2[k][i] = xor(F2[k][i],T[j]);
					}
				}
			}
			// add P_i,2 part to F_i,2
			for (j = 0; j < OIL_VARS; j++) {
				r = Next_Piece(&PG);
				for (k = 0; k < 16; k++) {
					if (x*16 + k >= OIL_VARS)
						break;
					if ( r & ((uint16_t) 1 ) << k ) {
						flipBit(&F2[k][i],j);
					}
				}
			}
		}

		// Calculate v*P_i,2 and assign to the i-th row of the LHS of the augmented matrix
		for (k = 0; k < 16; k++)	{
			if (16*x+k>= OIL_VARS)
				break;
			for (i = 0; i <= VINEGAR_VARS; i++) {
				for (j = 0; j < OIL_VARS; j++)	{
					if (getBit(F2[k][i],j)) {
						A.array[16*x+k][j] = add(A.array[16*x+k][j],vinegar_variables[i]);
					}
				}
			}
		}

	}

	// free the memory occupied by the F_i,2
	for (k = 0; k<16; k++){
		free(F2[k]);
	}
	free(F2);
}

/*
	Solves the system F(x) = target for x

	publicseed : The public seed used to generate the first part of the public map P
	T : The secret linear map
	target : The target vector to find a solution for
	solution : receives a solution
*/
static void solvePrivateUOVSystem(const unsigned char *publicseed, column *T, FELT *target , FELT *solution) {
	Matrix A;
	int solution_found = 0;

    // Repeatedly try an assignment to the vinegar variables until a unique solution is found
	while (solution_found == 0) {
		// Set homogenizing variable to one and pick random vinegar variables
		solution[0] = ONE;
		randombytes((unsigned char *) (solution + 1),VINEGAR_VARS*FIELD_SIZE/8);

		// Print vinegar values if KAT is defined
		//printVinegarValues(&(solution[1]));

		// Build the augmented matrix for the linear system
		A = zeroMatrix(OIL_VARS, OIL_VARS + 1);
		BuildAugmentedMatrix(A, solution , target, T , publicseed);

		// Print augmented matrix if KAT is defined
		//printAugmentedMatrix(A);

		// Try to find a unique solution to the linear system
		solution_found = getUniqueSolution(A,&(solution[1+VINEGAR_VARS]));

		// Report whether a solution is found if KAT is defined
		//reportSolutionFound(solution_found);

		// Free the memory occupied by the augmented matrix
		destroy_matrix(A);
	}
}

#ifndef MESSAGE_RECOVERY

/*
	Computes the target vector by hashing the document, after padding it with a 0x00 byte and a 16-byte salt.
	(Only used in appended signature mode)

	document : The document that is being signed
	len : The number of bytes of the document being signed
	target : receives the target vector
	salt : The 16-byte salt
 */
static void computeTarget(const unsigned char* document , uint64_t len, FELT *target, const unsigned char* salt){
	Sponge sponge;
	unsigned char pad = 0;

    shake_inc_init(&sponge);
    shake_inc_absorb(&sponge, document, len);
    shake_inc_absorb(&sponge, &pad, 1);
    shake_inc_absorb(&sponge, salt, SALT_BYTES);
    shake_inc_finalize(&sponge);
    shake_inc_squeeze((unsigned char *)target, sizeof(FELT)*OIL_VARS, &sponge);
}
#else

/*
	Computes the target vector.
	The document is hashed, after padding it with a 0x01 byteand a 16-byte salt, to get the first part of the target.
	Then, the first part is hashed again and xored with the last part of the padded document to get the second part of the target.
	(Only used in message recovery mode)

	document : The document that is being signed
	len : The number of bytes of the document being signed
	target : receives the target vector
	salt: The 16-byte salt
 */
void computeTarget(const unsigned char* document , uint64_t len, FELT *target, const unsigned char* salt){
	int i,start_recovery;
	Sponge sponge;
	unsigned char *buf = (unsigned char *) target;
	unsigned char pad = 1;

	// Compute first part of the target and put in the first part of the buffer
    shake_inc_init(&sponge);
    shake_inc_absorb(&sponge, document, len);
    shake_inc_absorb(&sponge, &pad, 1);
    shake_inc_absorb(&sponge, salt, SALT_BYTES);
    shake_inc_finalize(&sponge);
	squeezeBytes(&sponge  , buf , FIRST_PART_TARGET );

	// Absorb first part of target into a Sponge object and squeeze into the second part of the buffer
	initializeAndAbsorb(&sponge, buf , FIRST_PART_TARGET);
	squeezeBytes(&sponge  , &(buf[FIRST_PART_TARGET]) , SECOND_PART_TARGET );

	// If not the entire document can be covered from a signature, we xor the last part of the message
	// into the second part of the buffer and we xor the last byte with a 0x01.
	// Otherwise, we xor the entire document into the second part of the buffer, and we xor the next byte with a 0x01
	if(len > RECOVERED_PART_MESSAGE){
		start_recovery = len- RECOVERED_PART_MESSAGE;
		buf[FIRST_PART_TARGET + SECOND_PART_TARGET - 1] ^= 1;
	}
	else{
		start_recovery = 0;
		buf[FIRST_PART_TARGET + len] ^= 1;
	}
	for(i = start_recovery ; i<len ; i++){
		buf[FIRST_PART_TARGET + i-start_recovery] ^= document[i];
	}
}
#endif

#ifdef MESSAGE_RECOVERY
/*
	If message recovery is enabled, this function extracts the last part of the document from the evaluated signature and appends it to the first part of document

	document : Initially this contains the first part of the document, after the call to this function this contains the entire original document
	len      : pointer to the length of document, which is altered appropriately
	evaluation : The evaluation of the public map in the signature
*/
static void extractMessage(unsigned char *document ,unsigned long long *len , FELT *evaluation){
	int i, reading;
	unsigned char buf2[SECOND_PART_TARGET];
	Sponge sponge;

	unsigned char *eval = (unsigned char *) evaluation;

	// Absorb the first part of the buffer into a Sponge object and squeeze into buffer 2
	initializeAndAbsorb(&sponge, eval , FIRST_PART_TARGET );
	squeezeBytes(&sponge  , buf2 , SECOND_PART_TARGET );

	// Xor the secon part of the evaluation into buffer 2
	for(i = 0 ; i<SECOND_PART_TARGET ; i++ ){
		buf2[i] ^= eval[FIRST_PART_TARGET + i];
	}

	// Start searching from the left for the first byte equal to 0x01
	// All bytes before this byte are appended to the document and len is increased
	reading = 0;
	unsigned long long oldlen = *len;
	for (i = SECOND_PART_TARGET-1; i >= 0 ; i--)
	{
		if(reading){
			document[oldlen + i] = buf2[i];
		}
		else{
			if(buf2[i] == 1){
				reading = 1;
				*len += i;
			}
		}
	}
}
#endif

/*
	Generates a signature for a document

	sm : receives the signed message
	smlen : receives the length of the signed message
	m : the message to be signed
	mlen : the length oth the message to be signed
	sk : the secret key
*/
int luov_sign(unsigned char *sm, size_t *smlen, const unsigned char *m , size_t mlen,  const unsigned char *sk) {
	int i, j;
	FELT target[OIL_VARS];

	//printIntermediateValue("--- Start signing ---\n");

	unsigned char *sig = sm;

	// If not the entire mesage can be recovered from a signature, we copy the first part to sm.
	if( mlen > RECOVERED_PART_MESSAGE ){
		memcpy(sm,m,mlen - RECOVERED_PART_MESSAGE);
		sig += mlen - RECOVERED_PART_MESSAGE;
	}

	// pick random salt
	randombytes(SIG_SALT(sig),SALT_BYTES);

	// calculate public seed and T from private seed
	Sponge sponge;
	initializeAndAbsorb(&sponge, sk , 32);
	unsigned char publicseed[32];
	squeezeBytes(&sponge, publicseed , 32);

	column T[VINEGAR_VARS +1];
	T[0] = empty; // this makes the linear map T linear, picking the first row of T random would make T an affine transformation
	squeeze_column_array(&sponge, &(T[1]) , VINEGAR_VARS);

	// compute the target for the public map P
	computeTarget(m, mlen , target, SIG_SALT(sig));

	FELT solution[VARS+1];

	// Generate a solution to F(x) = target
	solvePrivateUOVSystem(publicseed, T, target, solution);

	// Print solution to the equation F(x) = target if KAT is defined
	//printPrivateSolution(solution);

	// Convert into a solution for P(x) = target
	for (i = 0; i <= VINEGAR_VARS; i++) {
		for (j = 0; j < OIL_VARS; j++) {
			if ( getBit(T[i] , j )) {
				solution[i] = subtract(solution[i],solution[VINEGAR_VARS +1+ j]);
			}
		}
	}

	// Copy the solution to the signature
	memcpy(SIG_SOL(sig),(unsigned char *) (solution+1),VARS*sizeof(FELT));

	// Ste the length of the signed message
	*smlen = mlen + CRYPTO_BYTES - RECOVERED_PART_MESSAGE;
	if(*smlen < CRYPTO_BYTES)
		*smlen = CRYPTO_BYTES;

	//printIntermediateValue("--- End signing ---\n");
	return 0;
}

/*
	Evaluated the public map P in a signature

	pk : The public key
	sig : The point that P is evaluated in
	evaluation : Receives the vector P(signature)
*/
static void evaluatePublicMap(const unsigned char *pk, const unsigned char *sig , FELT* evaluation){
	int i,j,k;
	FELT prod;
	column r;

	// Add a homogenizing variable equal to 1
	FELT solution[VARS+1];
	solution[0] = ONE;
	memcpy((unsigned char *)(solution+1), SIG_SOL(sig) , VARS*sizeof(FELT));

	// initialize evaluation to zero
	for(i = 0 ; i<OIL_VARS ; i++){
		evaluation[i] = ZERO;
	}

	// Initialise PRNG to generate the first part of the public map P
	ColumnGenerator CG;
	ColumnGenerator_init(&CG,PK_SEED(pk));

	// Evaluate the part of P that is generated from the public seed
	for (i = 0; i <= VINEGAR_VARS; i++) {
		for (j = i; j <= VARS ; j++) {
			r = Next_Column(&CG);

			prod = multiply(solution[i], solution[j]);
			for (k = 0; k < OIL_VARS; k++) {
				if (getBit(r,k)) {
					evaluation[k] = add(evaluation[k],prod);
				}
			}
		}
	}

	// Evaluate the part of P that is stored in the public key
	reader R = newReader(PK_Q2(pk));
	for (i = VINEGAR_VARS +1 ; i <= VARS; i++) {
		for (j = i; j <= VARS; j++) {
			prod = multiply(solution[i], solution[j]);
			column col = deserialize_column(&R);
			for (k = 0; k < OIL_VARS; k++) {
				if ( getBit(col , k) ) {
					evaluation[k] = add(evaluation[k],prod);
				}
			}
		}
	}

	// prints the evaluation of the public map if KAT is defined
	//printEvaluation(evaluation);
}

/*
	Verifies a signature for a document

	Remark : If we are in message recovery mode, this function does more work than strictly necessary

	m : receives the message
	mlen : receives the length of the message
	sm : the signed message to verify
	smlen : length ot the signed message
	pk : the public key

	returns : 0 if the signature is valid, -1 otherwise
*/
int luov_verify(unsigned char *m , size_t *mlen, const unsigned char *sm, size_t smlen , const unsigned char *pk ) {
	int i;
	FELT evaluation[OIL_VARS];
	FELT target[OIL_VARS];

	//printIntermediateValue("--- Start verifying ---\n");

	// reject if the signed message is too short
	if (smlen < CRYPTO_BYTES)
		return -1;

	// Copy the part of the message that cannot be recovered from the signature into m
	memcpy(m,sm,smlen-CRYPTO_BYTES);
	*mlen = smlen - CRYPTO_BYTES;
	const unsigned char *sig = sm + smlen-CRYPTO_BYTES;

	// Evaluate the public map P at the signature
	evaluatePublicMap(pk, sig , evaluation);

	#ifdef MESSAGE_RECOVERY
	// If we are in message recovery mode, we extracts a part of the document from the signature
	extractMessage(m , mlen , evaluation);
	#endif

	// We compute the target based on the full document
	computeTarget(m, *mlen, target, SIG_SALT(sig));

	// Output 0 if the evaluation of the public map is equal to the target, otherwise output -1
	for(i=0 ; i<OIL_VARS ; i++){
		if (! isEqual(target[i], evaluation[i])){
			return -1;
		}
	}

	//printIntermediateValue("--- End verifying ---\n");

	return 0;
}
