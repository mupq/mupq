/**
 * \file gabidulin.c
 * \brief Implementation of gabidulin.h
 *
 * The decoding algorithm provided is based on q_polynomials reconstruction, see \cite gabidulin:welch and \cite gabidulin:generalized for details.
 *
 */

#include "ffi_elt.h"
#include "ffi_vec.h"
#include "gabidulin.h"
#include "parameters.h"
#include "q_polynomial.h"


/**
 * \fn gabidulin_code gabidulin_code_init(const ffi_poly g, unsigned int k, unsigned int n)
 * \brief This function initializes a gabidulin code
 *
 * \param[in] g Generator polytor defining the code
 * \param[in] k Size of polytors representing messages
 * \param[in] n Size of vetors representing codewords
 * \return Gabidulin code
 */
gabidulin_code gabidulin_code_init(const ffi_poly g, unsigned int k, unsigned int n) {
  gabidulin_code code;

  code.g = g;
  code.k = k;
  code.n = n;

  return code;
}



/**
 * \fn void gabidulin_code_encode(ffi_qre c, gabidulin_code gc, const ffi_vec m)
 * \brief This function encodes a message into a codeword
 *
 * This function assumes <b>FIELD_Q</b> = 2
 *
 * \param[out] c Vector of size <b>n</b> representing a codeword
 * \param[in] gc Gabidulin code
 * \param[in] m Vector of size <b>k</b> representing a message
 */
void gabidulin_code_encode(ffi_qre c, gabidulin_code gc, const ffi_vec m) {

  // Compute generator matrix
  ffi_elt matrix[PARAM_K][PARAM_N];
  for(unsigned int j = 0 ; j < gc.n ; ++j) {
    ffi_elt_set(matrix[0][j], gc.g->v[j]);
    for(unsigned int i = 1 ; i < gc.k ; ++i) {
      ffi_elt_sqr(matrix[i][j], matrix[i-1][j]);
    }
  }

  // Encode message
  ffi_elt tmp;
  ffi_elt_set_zero(tmp);
  ffi_vec_set_zero(c->v, gc.n);
  for(unsigned int i = 0 ; i < gc.k ; ++i) {
    for(unsigned int j = 0 ; j < gc.n ; ++j) {
      ffi_elt_mul(tmp, m[i], matrix[i][j]);
      ffi_elt_add(c->v[j], c->v[j], tmp);
    }
  }
}



/**
 * \fn void gabidulin_code_decode(ffi_vec m, gabidulin_code gc, const ffi_qre y)
 * \brief This function decodes a word
 *
 * As explained in the supporting documentation, the provided decoding algorithm works as follows (see \cite gabidulin:welch and \cite gabidulin:generalized for details):
 *   1. Find a solution (<b>V</b>, <b>N</b>) of the q-polynomial Reconstruction2(<b>y</b>, <b>gc.g</b>, <b>gc.k</b>, (<b>gc.n</b> - <b>gc.k</b>)/2) problem using \cite gabidulin:generalized (section 4, algorithm 5) ;
 *   2. Find <b>f</b> by computing <b>V \ (N.A) + I</b> (see "Polynomials with lower degree" improvement from \cite gabidulin:generalized, section 4.4.2) ;
 *   3. Retrieve the message <b>m</b> as the k first coordinates of <b>f</b>.
 *
 *  This function assumes <b>FIELD_Q</b> = 2
 *
 * \param[out] m Vector of size <b>k</b> representing a message
 * \param[in] gc Gabidulin code
 * \param[in] y Vector of size <b>n</b> representing a word to decode
 */
void gabidulin_code_decode(ffi_vec m, gabidulin_code gc, const ffi_qre y) {

  /*
   *  Step 1: Solving the q-polynomial reconstruction2 problem
   */

  int t = (gc.n - gc.k) / 2;
  int max_degree_N = (gc.n - gc.k) % 2 == 0 ? gc.k + t - 1 : gc.k + t;

  q_polynomial* A = q_polynomial_init(gc.k);
  q_polynomial* I = q_polynomial_init(gc.k - 1);

  q_polynomial* N0 = q_polynomial_init(max_degree_N);
  q_polynomial* N1 = q_polynomial_init(max_degree_N);
  q_polynomial* V0 = q_polynomial_init(t);
  q_polynomial* V1 = q_polynomial_init(t);

  q_polynomial* qtmp1 = q_polynomial_init(max_degree_N);
  q_polynomial* qtmp2 = q_polynomial_init(max_degree_N);
  q_polynomial* qtmp3 = q_polynomial_init(t);
  q_polynomial* qtmp4 = q_polynomial_init(t);

  ffi_vec u0, u1;
  ffi_vec_init(&u0, PARAM_N);
  ffi_vec_init(&u1, PARAM_N);

  ffi_elt e1, e2, tmp1, tmp2;


  // Initialization step

  // A(g[i]) = 0 for 0 <= i <= k - 1
  // I(g[i]) = y[i] for 0 <= i <= k - 1

  q_polynomial_set_interpolate_vect_and_zero(A, I, gc.g, y, gc.k);
  q_polynomial_set_one(N0);
  q_polynomial_set_zero(N1);
  q_polynomial_set_zero(V0);
  q_polynomial_set_one(V1);

  q_polynomial_set_zero(qtmp1);
  q_polynomial_set_zero(qtmp2);
  q_polynomial_set_zero(qtmp3);
  q_polynomial_set_zero(qtmp4);

  // u0[i] = A(g[i]) - V0(y[i])
  // u1[i] = I(g[i]) - V1(y[i])

  for(unsigned int i = 0 ; i < gc.n ; ++i) {
    q_polynomial_evaluate(tmp1, A, gc.g->v[i]);
    q_polynomial_evaluate(tmp2, V0, y->v[i]);
    ffi_elt_add(tmp1, tmp1, tmp2);
    ffi_elt_set(u0[i], tmp1);

    q_polynomial_evaluate(tmp1, I, gc.g->v[i]);
    q_polynomial_evaluate(tmp2, V1, y->v[i]);
    ffi_elt_add(tmp1, tmp1, tmp2);
    ffi_elt_set(u1[i], tmp1);
  }

  // Interpolation step
  int updateType = -1;
  for(unsigned int i = gc.k ; i < gc.n ; ++i) {

    unsigned int j = i;
    while(j < gc.n &&
          ffi_elt_is_zero(u0[j]) == 0 &&
          ffi_elt_is_zero(u1[j]) == 1) {
      j++;
    }

    if(j == gc.n) {
      break;
    } else {
      if(i != j) {
        // Permutation of the coordinates of positions i and j
        ffi_elt_set(tmp1, u0[i]);
        ffi_elt_set(u0[i], u0[j]);
        ffi_elt_set(u0[j], tmp1);

        ffi_elt_set(tmp1, u1[i]);
        ffi_elt_set(u1[i], u1[j]);
        ffi_elt_set(u1[j], tmp1);
      }
    }


    // Update q_polynomials according to discrepancies
    if(ffi_elt_is_zero(u1[i]) != 1) {
      updateType = 1;

      // e1 = - u1[i]^q / u1[i]
      // e2 = - u0[i] / u1[i]
      // N0' = N1^q - e1.N1
      // V0' = V1^q - e1.V1
      // N1' = N0 - e2.N1
      // V1' = V0 - e2.V1

      ffi_elt_inv(tmp1, u1[i]);
      ffi_elt_sqr(e1, u1[i]);
      ffi_elt_mul(e1, e1, tmp1);
      ffi_elt_mul(e2, u0[i], tmp1);

      q_polynomial_scalar_mul(qtmp1, N1, e1);
      q_polynomial_qexp(qtmp2, N1);
      q_polynomial_scalar_mul(qtmp3, V1, e1);
      q_polynomial_qexp(qtmp4, V1);

      q_polynomial_scalar_mul(N1, N1, e2);
      q_polynomial_add(N1, N0, N1);

      q_polynomial_scalar_mul(V1, V1, e2);
      q_polynomial_add(V1, V0, V1);

      q_polynomial_add(N0, qtmp1, qtmp2);
      q_polynomial_add(V0, qtmp3, qtmp4);
    }

    if(ffi_elt_is_zero(u0[i]) == 1 &&
       ffi_elt_is_zero(u1[i]) == 1) {
      updateType = 2;

      // N0' = N1^q
      // V0' = V1^q
      // N1' = N0
      // V1' = V0

      q_polynomial_qexp(qtmp1, N1);
      q_polynomial_qexp(qtmp2, V1);

      q_polynomial_set(N1, N0);
      q_polynomial_set(V1, V0);
      q_polynomial_set(N0, qtmp1);
      q_polynomial_set(V0, qtmp2);
    }


    // Update discrepancies
    for(unsigned int k = i + 1 ; k < gc.n ; ++k) {
      if(updateType == 1) {

        // u0[k]' = u1[k]^q - e1.u1[k]
        // u1[k]' = u0[k] - e2.u1[k]

        ffi_elt_mul(tmp1, e1, u1[k]);
        ffi_elt_sqr(tmp2, u1[k]);
        ffi_elt_add(tmp1, tmp1, tmp2);

        ffi_elt_mul(tmp2, e2, u1[k]);
        ffi_elt_add(tmp2, tmp2, u0[k]);
        ffi_elt_set(u1[k], tmp2);

        ffi_elt_set(u0[k], tmp1);
      }

      if(updateType == 2) {

        // u0[k]' = u0[k]
        // u1[k]' = u1[k]^q

        ffi_elt_sqr(tmp1, u1[k]);
        ffi_elt_set(u1[k], tmp1);
      }
    }
  }

  /*
   *  Step 2: Computing f (qtmp1 variable) using Loidreau's improvement for lower degree polynomials
   */

  q_polynomial_mul(qtmp1, N1, A);
  q_polynomial_left_div(qtmp3, qtmp2, qtmp1, V1);
  q_polynomial_add(qtmp1, qtmp3, I);


  /*
   *  Step 3: Decoding the message as the value of the k first coordinates of f (qtmp1 variable)
   */

  ffi_vec_set(m, qtmp1->values, gc.k);

  q_polynomial_clear(A, gc.k);
  q_polynomial_clear(I, gc.k - 1);

  q_polynomial_clear(N0, max_degree_N);
  q_polynomial_clear(N1, max_degree_N);
  q_polynomial_clear(V0, t);
  q_polynomial_clear(V1, t);

  q_polynomial_clear(qtmp1, max_degree_N);
  q_polynomial_clear(qtmp2, max_degree_N);
  q_polynomial_clear(qtmp3, t);
  q_polynomial_clear(qtmp4, t);

  ffi_vec_clear(u0, PARAM_N);
  ffi_vec_clear(u1, PARAM_N);
}

