/* dgetf2.f -- translated by f2c (version 19970805).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b6 = -1.;

/* Subroutine */ int dgetf2_(m, n, a, lda, ipiv, info)
integer *m, *n;
doublereal *a;
integer *lda, *ipiv, *info;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    extern /* Subroutine */ int dger_();
    static integer j;
    extern /* Subroutine */ int dscal_(), dswap_();
    static integer jp;
    extern integer idamax_();
    extern /* Subroutine */ int xerbla_();


/*  -- LAPACK routine (version 3.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     June 30, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DGETF2 computes an LU factorization of a general m-by-n matrix A */
/*  using partial pivoting with row interchanges. */

/*  The factorization has the form */
/*     A = P * L * U */
/*  where P is a permutation matrix, L is lower triangular with unit */
/*  diagonal elements (lower trapezoidal if m > n), and U is upper */
/*  triangular (upper trapezoidal if m < n). */

/*  This is the right-looking Level 2 BLAS version of the algorithm. */

/*  Arguments */
/*  ========= */

/*  M       (input) INTEGER */
/*          The number of rows of the matrix A.  M >= 0. */

/*  N       (input) INTEGER */
/*          The number of columns of the matrix A.  N >= 0. */

/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*          On entry, the m by n matrix to be factored. */
/*          On exit, the factors L and U from the factorization */
/*          A = P*L*U; the unit diagonal elements of L are not stored. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array A.  LDA >= max(1,M). */

/*  IPIV    (output) INTEGER array, dimension (min(M,N)) */
/*          The pivot indices; for 1 <= i <= min(M,N), row i of the */
/*          matrix was interchanged with row IPIV(i). */

/*  INFO    (output) INTEGER */
/*          = 0: successful exit */
/*          < 0: if INFO = -k, the k-th argument had an illegal value */
/*          > 0: if INFO = k, U(k,k) is exactly zero. The factorization */
/*               has been completed, but the factor U is exactly */
/*               singular, and division by zero will occur if it is used 
*/
/*               to solve a system of equations. */

/*  ===================================================================== 
*/

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --ipiv;

    /* Function Body */
    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < max(1,*m)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGETF2", &i__1, 6L);
	return 0;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
	return 0;
    }

    i__1 = min(*m,*n);
    for (j = 1; j <= i__1; ++j) {

/*        Find pivot and test for singularity. */

	i__2 = *m - j + 1;
	jp = j - 1 + idamax_(&i__2, &a[j + j * a_dim1], &c__1);
	ipiv[j] = jp;
	if (a[jp + j * a_dim1] != 0.) {

/*           Apply the interchange to columns 1:N. */

	    if (jp != j) {
		dswap_(n, &a[j + a_dim1], lda, &a[jp + a_dim1], lda);
	    }

/*           Compute elements J+1:M of J-th column. */

	    if (j < *m) {
		i__2 = *m - j;
		d__1 = 1. / a[j + j * a_dim1];
		dscal_(&i__2, &d__1, &a[j + 1 + j * a_dim1], &c__1);
	    }

	} else if (*info == 0) {

	    *info = j;
	}

	if (j < min(*m,*n)) {

/*           Update trailing submatrix. */

	    i__2 = *m - j;
	    i__3 = *n - j;
	    dger_(&i__2, &i__3, &c_b6, &a[j + 1 + j * a_dim1], &c__1, &a[j + (
		    j + 1) * a_dim1], lda, &a[j + 1 + (j + 1) * a_dim1], lda);
	}
/* L10: */
    }
    return 0;

/*     End of DGETF2 */

} /* dgetf2_ */

