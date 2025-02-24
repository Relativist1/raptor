/*
 * Some auxillary functions (some might be repeated?) form ipole
 */

#include "definitions.h"
#include "functions.h"
#include "global_vars.h"
#include "model_definitions.h"
#include "model_functions.h"
#include "model_global_vars.h"

#define delta(i, j) ((i) == (j))

void bl_coord(double *X, double *r, double *th)
{
    // printf("%f/n", R0);
    *r = exp(X[1]);
    *th = X[2];

    return;
}

void gcov_bl(double r, double theta, double g_dd[NDIM][NDIM]){
    double sint = sin(theta);
    double cost = cos(theta);
    double sigma = r * r + a * a * cost * cost;
    double delta = r * r + a * a - 2. * r;
    double A_ = (r * r + a * a) * (r * r + a * a) - delta * a * a * sint * sint;
    LOOP_ij g_dd[i][j] =0.;
    // Covariant metric elements
    g_dd[0][0] = -(1. - 2. * r / sigma);
    g_dd[1][1] = sigma / delta ;
    g_dd[2][2] = sigma;
    g_dd[3][3] = A_ / sigma * sint * sint;
    g_dd[0][3] = -2. * a * r * sint * sint / sigma;
    g_dd[3][0] = g_dd[0][3];
}
void gcon_bl(double r, double theta, double g_uu[NDIM][NDIM]){

    double sint = sin(theta);
    double cost = cos(theta);
    double sigma = r * r + a * a * cost * cost;
    double delta = r * r + a * a - 2. * r;
    double A_ = (r * r + a * a) * (r * r + a * a) - delta * a * a * sint * sint;
    LOOP_ij g_uu[i][j] = 0.;
    // Contravariant metric elements
    g_uu[0][0] = -A_ / (sigma * delta);
    g_uu[1][1] = delta / sigma;
    g_uu[2][2] = 1. / sigma;
    g_uu[3][3] = (delta - a * a * sint * sint) / (sigma * delta * sint * sint);
    g_uu[0][3] = -2. * a * r / (sigma * delta);
    g_uu[3][0] = g_uu[0][3];

}

void bl_raise(double g_uu[NDIM][NDIM], double V_u[4], double V_d[4]) {

    V_u[0] = 0.;
    V_u[1] = 0.;
    V_u[2] = 0.;
    V_u[3] = 0.;

    LOOP_ij V_u[i] += g_uu[i][j] * V_d[j];
}

void bl_lower(double g_dd[NDIM][NDIM], double V_d[4], double V_u[4]) {

    V_d[0] = 0.;
    V_d[1] = 0.;
    V_d[2] = 0.;
    V_d[3] = 0.;

    LOOP_ij V_d[i] += g_dd[i][j] * V_u[j];
}

void normalize(double vcon[NDIM], double Gcov[NDIM][NDIM]){
  int k, l;
  double norm;

  norm = 0.;
  for (k = 0; k < 4; k++)
    for (l = 0; l < 4; l++)
      norm += vcon[k] * vcon[l] * Gcov[k][l];

  norm = sqrt(fabs(norm));
  for (k = 0; k < 4; k++)
    vcon[k] /= norm;

  return;
}

void set_dxdX(double X[NDIM], double dxdX[NDIM][NDIM])
{

  LOOP_ij dxdX[i][j] = delta(i, j);

    dxdX[1][1] = exp(X[1]);
    dxdX[2][2] = M_PI;
        

  }
void set_dXdx(double X[NDIM], double dXdx[NDIM][NDIM]) {
  double dxdX[NDIM][NDIM];
  set_dxdX(X, dxdX);
  invert_matrix(dxdX, dXdx);
}
void vec_from_ks(double X[NDIM], double v_ks[NDIM], double v_nat[NDIM]) {
  double dXdx[NDIM][NDIM];
  set_dXdx(X, dXdx);

  LOOP_i v_nat[i] = 0.;
  LOOP_ij v_nat[i] += dXdx[i][j] * v_ks[j];
}

void BL_to_KS(double X[NDIM], double BLphoton_u[NDIM], double KSphoton_u[NDIM]) {
    // This just changes vec and not coord.
    double trans[4][4];
    double X_u[4], U_u[4];
    double r, th;
    bl_coord(X, &r, &th);

    LOOP_i {
        U_u[i] = BLphoton_u[i];
    }

    LOOP_ij trans[i][j] = 0.;
    LOOP_i trans[i][i] = 1.;

    trans[0][1] = 2. * r / (r * r - 2. * r + a * a);
    trans[3][1] = a / (r * r - 2. * r + a * a);

    double U_u_dummy[4], X_u_dummy[4];
    LOOP_i {
        U_u_dummy[i] = U_u[i];
        U_u[i] = 0.;
    }

    LOOP_ij U_u[i] += trans[i][j] * U_u_dummy[j];

    LOOP_i {
        KSphoton_u[i] = U_u[i];
    }
}

/*
 LU_decompose():
 Performs a LU decomposition of the matrix A using Crout's method
 with partial implicit pivoting.  The exact LU decomposition of the
 matrix can be reconstructed from the resultant row-permuted form via
 the integer array permute[]
 The algorithm closely follows ludcmp.c of "Numerical Recipes
 in C" by Press et al. 1992.
 This will be used to solve the linear system  A.x = B
 Returns (1) if a singular matrix is found,  (0) otherwise.
*/
int LU_decompose(double A[][NDIM], int permute[])
{

  const double absmin = 1.e-30; /* Value used instead of 0 for singular matrices */

  //static double row_norm[NDIM];
  double row_norm[NDIM];
  double absmax, maxtemp;

  int i, j, k, max_row;
  int n = NDIM;

  max_row = 0;

  /* Find the maximum elements per row so that we can pretend later
   we have unit-normalized each equation: */

  for (i = 0; i < n; i++) {
    absmax = 0.;

    for (j = 0; j < n; j++) {

      maxtemp = fabs(A[i][j]);

      if (maxtemp > absmax) {
        absmax = maxtemp;
      }
    }

    /* Make sure that there is at least one non-zero element in this row: */
    if (absmax == 0.) {
      //fprintf(stderr, "LU_decompose(): row-wise singular matrix!\n");
      return (1);
    }

    row_norm[i] = 1. / absmax; /* Set the row's normalization factor. */
  }

  /* The following the calculates the matrix composed of the sum
   of the lower (L) tridagonal matrix and the upper (U) tridagonal
   matrix that, when multiplied, form the original maxtrix.
   This is what we call the LU decomposition of the maxtrix.
   It does this by a recursive procedure, starting from the
   upper-left, proceding down the column, and then to the next
   column to the right.  The decomposition can be done in place
   since element {i,j} require only those elements with {<=i,<=j}
   which have already been computed.
   See pg. 43-46 of "Num. Rec." for a more thorough description.
   */

  /* For each of the columns, starting from the left ... */
  for (j = 0; j < n; j++) {

    /* For each of the rows starting from the top.... */

    /* Calculate the Upper part of the matrix:  i < j :   */
    for (i = 0; i < j; i++) {
      for (k = 0; k < i; k++) {
        A[i][j] -= A[i][k] * A[k][j];
      }
    }

    absmax = 0.0;

    /* Calculate the Lower part of the matrix:  i <= j :   */

    for (i = j; i < n; i++) {

      for (k = 0; k < j; k++) {
        A[i][j] -= A[i][k] * A[k][j];
      }

      /* Find the maximum element in the column given the implicit
       unit-normalization (represented by row_norm[i]) of each row:
       */
      maxtemp = fabs(A[i][j]) * row_norm[i];

      if (maxtemp >= absmax) {
        absmax = maxtemp;
        max_row = i;
      }

    }

    /* Swap the row with the largest element (of column j) with row_j.  absmax
     This is the partial pivoting procedure that ensures we don't divide
     by 0 (or a small number) when we solve the linear system.
     Also, since the procedure starts from left-right/top-bottom,
     the pivot values are chosen from a pool involving all the elements
     of column_j  in rows beneath row_j.  This ensures that
     a row  is not permuted twice, which would mess things up.
     */
    if (max_row != j) {

      /* Don't swap if it will send a 0 to the last diagonal position.
       Note that the last column cannot pivot with any other row,
       so this is the last chance to ensure that the last two
       columns have non-zero diagonal elements.
       */

      if ((j == (n - 2)) && (A[j][j + 1] == 0.)) {
        max_row = j;
      } else {
        for (k = 0; k < n; k++) {

          maxtemp = A[j][k];
          A[j][k] = A[max_row][k];
          A[max_row][k] = maxtemp;

        }

        /* Don't forget to swap the normalization factors, too...
         but we don't need the jth element any longer since we
         only look at rows beneath j from here on out.
         */
        row_norm[max_row] = row_norm[j];
      }
    }

    /* Set the permutation record s.t. the j^th element equals the
     index of the row swapped with the j^th row.  Note that since
     this is being done in successive columns, the permutation
     vector records the successive permutations and therefore
     index of permute[] also indexes the chronology of the
     permutations.  E.g. permute[2] = {2,1} is an identity
     permutation, which cannot happen here though.
     */

    permute[j] = max_row;

    if (A[j][j] == 0.) {
      A[j][j] = absmin;
    }

    /* Normalize the columns of the Lower tridiagonal part by their respective
     diagonal element.  This is not done in the Upper part because the
     Lower part's diagonal elements were set to 1, which can be done w/o
     any loss of generality.
     */
    if (j != (n - 1)) {
      maxtemp = 1. / A[j][j];

      for (i = (j + 1); i < n; i++) {
        A[i][j] *= maxtemp;
      }
    }

  }

  return (0);

  /* End of LU_decompose() */

}

/*
 LU_substitution():
 Performs the forward (w/ the Lower) and backward (w/ the Upper)
 substitutions using the LU-decomposed matrix A[][] of the original
 matrix A' of the linear equation:  A'.x = B.  Upon entry, A[][]
 is the LU matrix, B[] is the source vector, and permute[] is the
 array containing order of permutations taken to the rows of the LU
 matrix.  See LU_decompose() for further details.
 Upon exit, B[] contains the solution x[], A[][] is left unchanged.
*/
void LU_substitution(double A[][NDIM], double B[], int permute[])
{
  int i, j;
  int n = NDIM;
  double tmpvar;

  /* Perform the forward substitution using the LU matrix.
   */
  for (i = 0; i < n; i++) {

    /* Before doing the substitution, we must first permute the
     B vector to match the permutation of the LU matrix.
     Since only the rows above the currrent one matter for
     this row, we can permute one at a time.
     */
    tmpvar = B[permute[i]];
    B[permute[i]] = B[i];
    for (j = (i - 1); j >= 0; j--) {
      tmpvar -= A[i][j] * B[j];
    }
    B[i] = tmpvar;
  }

  /* Perform the backward substitution using the LU matrix.
   */
  for (i = (n - 1); i >= 0; i--) {
    for (j = (i + 1); j < n; j++) {
      B[i] -= A[i][j] * B[j];
    }
    B[i] /= A[i][i];
  }

  /* End of LU_substitution() */
}

int invert_matrix(double Am[][NDIM], double Aminv[][NDIM])
{

  int i, j;
  int n = NDIM;
  int permute[NDIM];
  double dxm[NDIM], Amtmp[NDIM][NDIM];

  for (i = 0; i < NDIM; i++)
    for (j = 0; j < NDIM; j++) {
      Amtmp[i][j] = Am[i][j];
    }

  // Get the LU matrix:
  if (LU_decompose(Amtmp, permute) != 0) {
    //fprintf(stderr, "invert_matrix(): singular matrix encountered! \n");
    return (1);
  }

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      dxm[j] = 0.;
    }
    dxm[i] = 1.;

    /* Solve the linear system for the i^th column of the inverse matrix: :  */
    LU_substitution(Amtmp, dxm, permute);

    for (j = 0; j < n; j++) {
      Aminv[j][i] = dxm[j];
    }

  }

  return (0);
}
