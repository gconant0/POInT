/* dormqr.f -- translated by f2c (version 20000121).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__65 = 65;

/* Subroutine */ int dormqr_(side, trans, m, n, k, a, lda, tau, c__, ldc, 
	work, lwork, info, side_len, trans_len)
char *side, *trans;
integer *m, *n, *k;
doublereal *a;
integer *lda;
doublereal *tau, *c__;
integer *ldc;
doublereal *work;
integer *lwork, *info;
ftnlen side_len;
ftnlen trans_len;
{
    /* System generated locals */
    address a__1[2];
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__3[2], i__4, 
	    i__5;
    char ch__1[2];

    /* Builtin functions */
    /* Subroutine */ int s_cat();

    /* Local variables */
    static logical left;
    static integer i__;
    static doublereal t[4160]	/* was [65][64] */;
    extern logical lsame_();
    static integer nbmin, iinfo, i1, i2, i3;
    extern /* Subroutine */ int dorm2r_();
    static integer ib, ic, jc, nb, mi, ni;
    extern /* Subroutine */ int dlarfb_();
    static integer nq, nw;
    extern /* Subroutine */ int dlarft_(), xerbla_();
    extern integer ilaenv_();
    static logical notran;
    static integer ldwork, lwkopt;
    static logical lquery;
    static integer iws;


/*  -- LAPACK routine (version 3.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     June 30, 1999 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DORMQR overwrites the general real M-by-N matrix C with */

/*                  SIDE = 'L'     SIDE = 'R' */
/*  TRANS = 'N':      Q * C          C * Q */
/*  TRANS = 'T':      Q**T * C       C * Q**T */

/*  where Q is a real orthogonal matrix defined as the product of k */
/*  elementary reflectors */

/*        Q = H(1) H(2) . . . H(k) */

/*  as returned by DGEQRF. Q is of order M if SIDE = 'L' and of order N */
/*  if SIDE = 'R'. */

/*  Arguments */
/*  ========= */

/*  SIDE    (input) CHARACTER*1 */
/*          = 'L': apply Q or Q**T from the Left; */
/*          = 'R': apply Q or Q**T from the Right. */

/*  TRANS   (input) CHARACTER*1 */
/*          = 'N':  No transpose, apply Q; */
/*          = 'T':  Transpose, apply Q**T. */

/*  M       (input) INTEGER */
/*          The number of rows of the matrix C. M >= 0. */

/*  N       (input) INTEGER */
/*          The number of columns of the matrix C. N >= 0. */

/*  K       (input) INTEGER */
/*          The number of elementary reflectors whose product defines */
/*          the matrix Q. */
/*          If SIDE = 'L', M >= K >= 0; */
/*          if SIDE = 'R', N >= K >= 0. */

/*  A       (input) DOUBLE PRECISION array, dimension (LDA,K) */
/*          The i-th column must contain the vector which defines the */
/*          elementary reflector H(i), for i = 1,2,...,k, as returned by */
/*          DGEQRF in the first k columns of its array argument A. */
/*          A is modified by the routine but restored on exit. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array A. */
/*          If SIDE = 'L', LDA >= max(1,M); */
/*          if SIDE = 'R', LDA >= max(1,N). */

/*  TAU     (input) DOUBLE PRECISION array, dimension (K) */
/*          TAU(i) must contain the scalar factor of the elementary */
/*          reflector H(i), as returned by DGEQRF. */

/*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*          On entry, the M-by-N matrix C. */
/*          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q. */

/*  LDC     (input) INTEGER */
/*          The leading dimension of the array C. LDC >= max(1,M). */

/*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK) */
/*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */

/*  LWORK   (input) INTEGER */
/*          The dimension of the array WORK. */
/*          If SIDE = 'L', LWORK >= max(1,N); */
/*          if SIDE = 'R', LWORK >= max(1,M). */
/*          For optimum performance LWORK >= N*NB if SIDE = 'L', and */
/*          LWORK >= M*NB if SIDE = 'R', where NB is the optimal */
/*          blocksize. */

/*          If LWORK = -1, then a workspace query is assumed; the routine */
/*          only calculates the optimal size of the WORK array, returns */
/*          this value as the first entry of the WORK array, and no error */
/*          message related to LWORK is issued by XERBLA. */

/*  INFO    (output) INTEGER */
/*          = 0:  successful exit */
/*          < 0:  if INFO = -i, the i-th argument had an illegal value */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --tau;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1 * 1;
    c__ -= c_offset;
    --work;

    /* Function Body */
    *info = 0;
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
    lquery = *lwork == -1;

/*     NQ is the order of Q and NW is the minimum dimension of WORK */

    if (left) {
	nq = *m;
	nw = *n;
    } else {
	nq = *n;
	nw = *m;
    }
    if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*k < 0 || *k > nq) {
	*info = -5;
    } else if (*lda < max(1,nq)) {
	*info = -7;
    } else if (*ldc < max(1,*m)) {
	*info = -10;
    } else if (*lwork < max(1,nw) && ! lquery) {
	*info = -12;
    }

    if (*info == 0) {

/*        Determine the block size.  NB may be at most NBMAX, where NBMAX */
/*        is used to define the local array T. */

/* Computing MIN */
/* Writing concatenation */
	i__3[0] = 1, a__1[0] = side;
	i__3[1] = 1, a__1[1] = trans;
	s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
	i__1 = 64, i__2 = ilaenv_(&c__1, "DORMQR", ch__1, m, n, k, &c_n1, (
		ftnlen)6, (ftnlen)2);
	nb = min(i__1,i__2);
	lwkopt = max(1,nw) * nb;
	work[1] = (doublereal) lwkopt;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DORMQR", &i__1, (ftnlen)6);
	return 0;
    } else if (lquery) {
	return 0;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0 || *k == 0) {
	work[1] = 1.;
	return 0;
    }

    nbmin = 2;
    ldwork = nw;
    if (nb > 1 && nb < *k) {
	iws = nw * nb;
	if (*lwork < iws) {
	    nb = *lwork / ldwork;
/* Computing MAX */
/* Writing concatenation */
	    i__3[0] = 1, a__1[0] = side;
	    i__3[1] = 1, a__1[1] = trans;
	    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
	    i__1 = 2, i__2 = ilaenv_(&c__2, "DORMQR", ch__1, m, n, k, &c_n1, (
		    ftnlen)6, (ftnlen)2);
	    nbmin = max(i__1,i__2);
	}
    } else {
	iws = nw;
    }

    if (nb < nbmin || nb >= *k) {

/*        Use unblocked code */

	dorm2r_(side, trans, m, n, k, &a[a_offset], lda, &tau[1], &c__[
		c_offset], ldc, &work[1], &iinfo, (ftnlen)1, (ftnlen)1);
    } else {

/*        Use blocked code */

	if (left && ! notran || ! left && notran) {
	    i1 = 1;
	    i2 = *k;
	    i3 = nb;
	} else {
	    i1 = (*k - 1) / nb * nb + 1;
	    i2 = 1;
	    i3 = -nb;
	}

	if (left) {
	    ni = *n;
	    jc = 1;
	} else {
	    mi = *m;
	    ic = 1;
	}

	i__1 = i2;
	i__2 = i3;
	for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
	    i__4 = nb, i__5 = *k - i__ + 1;
	    ib = min(i__4,i__5);

/*           Form the triangular factor of the block reflector */
/*           H = H(i) H(i+1) . . . H(i+ib-1) */

	    i__4 = nq - i__ + 1;
	    dlarft_("Forward", "Columnwise", &i__4, &ib, &a[i__ + i__ * 
		    a_dim1], lda, &tau[i__], t, &c__65, (ftnlen)7, (ftnlen)10)
		    ;
	    if (left) {

/*              H or H' is applied to C(i:m,1:n) */

		mi = *m - i__ + 1;
		ic = i__;
	    } else {

/*              H or H' is applied to C(1:m,i:n) */

		ni = *n - i__ + 1;
		jc = i__;
	    }

/*           Apply H or H' */

	    dlarfb_(side, trans, "Forward", "Columnwise", &mi, &ni, &ib, &a[
		    i__ + i__ * a_dim1], lda, t, &c__65, &c__[ic + jc * 
		    c_dim1], ldc, &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (
		    ftnlen)7, (ftnlen)10);
/* L10: */
	}
    }
    work[1] = (doublereal) lwkopt;
    return 0;

/*     End of DORMQR */

} /* dormqr_ */

