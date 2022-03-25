/* dhseqr.f -- translated by f2c (version 19970805).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static doublereal c_b9 = 0.;
static doublereal c_b10 = 1.;
static integer c__4 = 4;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__8 = 8;
static integer c__15 = 15;
static logical c_false = FALSE_;
static integer c__1 = 1;

/* Subroutine */ int dhseqr_(job, compz, n, ilo, ihi, h__, ldh, wr, wi, z__, 
	ldz, work, lwork, info, job_len, compz_len)
char *job, *compz;
integer *n, *ilo, *ihi;
doublereal *h__;
integer *ldh;
doublereal *wr, *wi, *z__;
integer *ldz;
doublereal *work;
integer *lwork, *info;
ftnlen job_len;
ftnlen compz_len;
{
    /* System generated locals */
    address a__1[2];
    integer h_dim1, h_offset, z_dim1, z_offset, i__1, i__2, i__3[2], i__4, 
	    i__5;
    doublereal d__1, d__2;
    char ch__1[2];

    /* Builtin functions */
    /* Subroutine */ int s_cat();

    /* Local variables */
    static integer maxb;
    static doublereal absw;
    static integer ierr;
    static doublereal unfl, temp, ovfl;
    static integer i__, j, k, l;
    static doublereal s[225]	/* was [15][15] */, v[16];
    extern /* Subroutine */ int dscal_();
    extern logical lsame_();
    extern /* Subroutine */ int dgemv_();
    static integer itemp;
    extern /* Subroutine */ int dcopy_();
    static integer i1, i2;
    static logical initz, wantt, wantz;
    extern doublereal dlapy2_();
    extern /* Subroutine */ int dlabad_();
    static integer ii, nh;
    extern doublereal dlamch_();
    extern /* Subroutine */ int dlarfg_();
    static integer nr, ns;
    extern integer idamax_();
    static integer nv;
    extern doublereal dlanhs_();
    extern /* Subroutine */ int dlahqr_();
    static doublereal vv[16];
    extern /* Subroutine */ int dlacpy_();
    extern integer ilaenv_();
    extern /* Subroutine */ int dlaset_(), dlarfx_(), xerbla_();
    static doublereal smlnum;
    static logical lquery;
    static integer itn;
    static doublereal tau;
    static integer its;
    static doublereal ulp, tst1;


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

/*  DHSEQR computes the eigenvalues of a real upper Hessenberg matrix H */
/*  and, optionally, the matrices T and Z from the Schur decomposition */
/*  H = Z T Z**T, where T is an upper quasi-triangular matrix (the Schur 
*/
/*  form), and Z is the orthogonal matrix of Schur vectors. */

/*  Optionally Z may be postmultiplied into an input orthogonal matrix Q, 
*/
/*  so that this routine can give the Schur factorization of a matrix A */
/*  which has been reduced to the Hessenberg form H by the orthogonal */
/*  matrix Q:  A = Q*H*Q**T = (QZ)*T*(QZ)**T. */

/*  Arguments */
/*  ========= */

/*  JOB     (input) CHARACTER*1 */
/*          = 'E':  compute eigenvalues only; */
/*          = 'S':  compute eigenvalues and the Schur form T. */

/*  COMPZ   (input) CHARACTER*1 */
/*          = 'N':  no Schur vectors are computed; */
/*          = 'I':  Z is initialized to the unit matrix and the matrix Z 
*/
/*                  of Schur vectors of H is returned; */
/*          = 'V':  Z must contain an orthogonal matrix Q on entry, and */
/*                  the product Q*Z is returned. */

/*  N       (input) INTEGER */
/*          The order of the matrix H.  N >= 0. */

/*  ILO     (input) INTEGER */
/*  IHI     (input) INTEGER */
/*          It is assumed that H is already upper triangular in rows */
/*          and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally */
/*          set by a previous call to DGEBAL, and then passed to SGEHRD */
/*          when the matrix output by DGEBAL is reduced to Hessenberg */
/*          form. Otherwise ILO and IHI should be set to 1 and N */
/*          respectively. */
/*          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0. */

/*  H       (input/output) DOUBLE PRECISION array, dimension (LDH,N) */
/*          On entry, the upper Hessenberg matrix H. */
/*          On exit, if JOB = 'S', H contains the upper quasi-triangular 
*/
/*          matrix T from the Schur decomposition (the Schur form); */
/*          2-by-2 diagonal blocks (corresponding to complex conjugate */
/*          pairs of eigenvalues) are returned in standard form, with */
/*          H(i,i) = H(i+1,i+1) and H(i+1,i)*H(i,i+1) < 0. If JOB = 'E', 
*/
/*          the contents of H are unspecified on exit. */

/*  LDH     (input) INTEGER */
/*          The leading dimension of the array H. LDH >= max(1,N). */

/*  WR      (output) DOUBLE PRECISION array, dimension (N) */
/*  WI      (output) DOUBLE PRECISION array, dimension (N) */
/*          The real and imaginary parts, respectively, of the computed */
/*          eigenvalues. If two eigenvalues are computed as a complex */
/*          conjugate pair, they are stored in consecutive elements of */
/*          WR and WI, say the i-th and (i+1)th, with WI(i) > 0 and */
/*          WI(i+1) < 0. If JOB = 'S', the eigenvalues are stored in the 
*/
/*          same order as on the diagonal of the Schur form returned in */
/*          H, with WR(i) = H(i,i) and, if H(i:i+1,i:i+1) is a 2-by-2 */
/*          diagonal block, WI(i) = sqrt(H(i+1,i)*H(i,i+1)) and */
/*          WI(i+1) = -WI(i). */

/*  Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N) */
/*          If COMPZ = 'N': Z is not referenced. */
/*          If COMPZ = 'I': on entry, Z need not be set, and on exit, Z */
/*          contains the orthogonal matrix Z of the Schur vectors of H. */
/*          If COMPZ = 'V': on entry Z must contain an N-by-N matrix Q, */
/*          which is assumed to be equal to the unit matrix except for */
/*          the submatrix Z(ILO:IHI,ILO:IHI); on exit Z contains Q*Z. */
/*          Normally Q is the orthogonal matrix generated by DORGHR after 
*/
/*          the call to DGEHRD which formed the Hessenberg matrix H. */

/*  LDZ     (input) INTEGER */
/*          The leading dimension of the array Z. */
/*          LDZ >= max(1,N) if COMPZ = 'I' or 'V'; LDZ >= 1 otherwise. */

/*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK) 
*/
/*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */

/*  LWORK   (input) INTEGER */
/*          The dimension of the array WORK.  LWORK >= max(1,N). */

/*          If LWORK = -1, then a workspace query is assumed; the routine 
*/
/*          only calculates the optimal size of the WORK array, returns */
/*          this value as the first entry of the WORK array, and no error 
*/
/*          message related to LWORK is issued by XERBLA. */

/*  INFO    (output) INTEGER */
/*          = 0:  successful exit */
/*          < 0:  if INFO = -i, the i-th argument had an illegal value */
/*          > 0:  if INFO = i, DHSEQR failed to compute all of the */
/*                eigenvalues in a total of 30*(IHI-ILO+1) iterations; */
/*                elements 1:ilo-1 and i+1:n of WR and WI contain those */
/*                eigenvalues which have been successfully computed. */

/*  ===================================================================== 
*/

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

/*     Decode and test the input parameters */

    /* Parameter adjustments */
    h_dim1 = *ldh;
    h_offset = h_dim1 + 1;
    h__ -= h_offset;
    --wr;
    --wi;
    z_dim1 = *ldz;
    z_offset = z_dim1 + 1;
    z__ -= z_offset;
    --work;

    /* Function Body */
    wantt = lsame_(job, "S", 1L, 1L);
    initz = lsame_(compz, "I", 1L, 1L);
    wantz = initz || lsame_(compz, "V", 1L, 1L);

    *info = 0;
    work[1] = (doublereal) max(1,*n);
    lquery = *lwork == -1;
    if (! lsame_(job, "E", 1L, 1L) && ! wantt) {
	*info = -1;
    } else if (! lsame_(compz, "N", 1L, 1L) && ! wantz) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*ilo < 1 || *ilo > max(1,*n)) {
	*info = -4;
    } else if (*ihi < min(*ilo,*n) || *ihi > *n) {
	*info = -5;
    } else if (*ldh < max(1,*n)) {
	*info = -7;
    } else if (*ldz < 1 || wantz && *ldz < max(1,*n)) {
	*info = -11;
    } else if (*lwork < max(1,*n) && ! lquery) {
	*info = -13;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DHSEQR", &i__1, 6L);
	return 0;
    } else if (lquery) {
	return 0;
    }

/*     Initialize Z, if necessary */

    if (initz) {
	dlaset_("Full", n, n, &c_b9, &c_b10, &z__[z_offset], ldz, 4L);
    }

/*     Store the eigenvalues isolated by DGEBAL. */

    i__1 = *ilo - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wr[i__] = h__[i__ + i__ * h_dim1];
	wi[i__] = 0.;
/* L10: */
    }
    i__1 = *n;
    for (i__ = *ihi + 1; i__ <= i__1; ++i__) {
	wr[i__] = h__[i__ + i__ * h_dim1];
	wi[i__] = 0.;
/* L20: */
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }
    if (*ilo == *ihi) {
	wr[*ilo] = h__[*ilo + *ilo * h_dim1];
	wi[*ilo] = 0.;
	return 0;
    }

/*     Set rows and columns ILO to IHI to zero below the first */
/*     subdiagonal. */

    i__1 = *ihi - 2;
    for (j = *ilo; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = j + 2; i__ <= i__2; ++i__) {
	    h__[i__ + j * h_dim1] = 0.;
/* L30: */
	}
/* L40: */
    }
    nh = *ihi - *ilo + 1;

/*     Determine the order of the multi-shift QR algorithm to be used. */

/* Writing concatenation */
    i__3[0] = 1, a__1[0] = job;
    i__3[1] = 1, a__1[1] = compz;
    s_cat(ch__1, a__1, i__3, &c__2, 2L);
    ns = ilaenv_(&c__4, "DHSEQR", ch__1, n, ilo, ihi, &c_n1, 6L, 2L);
/* Writing concatenation */
    i__3[0] = 1, a__1[0] = job;
    i__3[1] = 1, a__1[1] = compz;
    s_cat(ch__1, a__1, i__3, &c__2, 2L);
    maxb = ilaenv_(&c__8, "DHSEQR", ch__1, n, ilo, ihi, &c_n1, 6L, 2L);
    if (ns <= 2 || ns > nh || maxb >= nh) {

/*        Use the standard double-shift algorithm */

	dlahqr_(&wantt, &wantz, n, ilo, ihi, &h__[h_offset], ldh, &wr[1], &wi[
		1], ilo, ihi, &z__[z_offset], ldz, info);
	return 0;
    }
    maxb = max(3,maxb);
/* Computing MIN */
    i__1 = min(ns,maxb);
    ns = min(i__1,15);

/*     Now 2 < NS <= MAXB < NH. */

/*     Set machine-dependent constants for the stopping criterion. */
/*     If norm(H) <= sqrt(OVFL), overflow should not occur. */

    unfl = dlamch_("Safe minimum", 12L);
    ovfl = 1. / unfl;
    dlabad_(&unfl, &ovfl);
    ulp = dlamch_("Precision", 9L);
    smlnum = unfl * (nh / ulp);

/*     I1 and I2 are the indices of the first row and last column of H */
/*     to which transformations must be applied. If eigenvalues only are 
*/
/*     being computed, I1 and I2 are set inside the main loop. */

    if (wantt) {
	i1 = 1;
	i2 = *n;
    }

/*     ITN is the total number of multiple-shift QR iterations allowed. */

    itn = nh * 30;

/*     The main loop begins here. I is the loop index and decreases from 
*/
/*     IHI to ILO in steps of at most MAXB. Each iteration of the loop */
/*     works with the active submatrix in rows and columns L to I. */
/*     Eigenvalues I+1 to IHI have already converged. Either L = ILO or */
/*     H(L,L-1) is negligible so that the matrix splits. */

    i__ = *ihi;
L50:
    l = *ilo;
    if (i__ < *ilo) {
	goto L170;
    }

/*     Perform multiple-shift QR iterations on rows and columns ILO to I 
*/
/*     until a submatrix of order at most MAXB splits off at the bottom */
/*     because a subdiagonal element has become negligible. */

    i__1 = itn;
    for (its = 0; its <= i__1; ++its) {

/*        Look for a single small subdiagonal element. */

	i__2 = l + 1;
	for (k = i__; k >= i__2; --k) {
	    tst1 = (d__1 = h__[k - 1 + (k - 1) * h_dim1], abs(d__1)) + (d__2 =
		     h__[k + k * h_dim1], abs(d__2));
	    if (tst1 == 0.) {
		i__4 = i__ - l + 1;
		tst1 = dlanhs_("1", &i__4, &h__[l + l * h_dim1], ldh, &work[1]
			, 1L);
	    }
/* Computing MAX */
	    d__2 = ulp * tst1;
	    if ((d__1 = h__[k + (k - 1) * h_dim1], abs(d__1)) <= max(d__2,
		    smlnum)) {
		goto L70;
	    }
/* L60: */
	}
L70:
	l = k;
	if (l > *ilo) {

/*           H(L,L-1) is negligible. */

	    h__[l + (l - 1) * h_dim1] = 0.;
	}

/*        Exit from loop if a submatrix of order <= MAXB has split off
. */

	if (l >= i__ - maxb + 1) {
	    goto L160;
	}

/*        Now the active submatrix is in rows and columns L to I. If 
*/
/*        eigenvalues only are being computed, only the active submatr
ix */
/*        need be transformed. */

	if (! wantt) {
	    i1 = l;
	    i2 = i__;
	}

	if (its == 20 || its == 30) {

/*           Exceptional shifts. */

	    i__2 = i__;
	    for (ii = i__ - ns + 1; ii <= i__2; ++ii) {
		wr[ii] = ((d__1 = h__[ii + (ii - 1) * h_dim1], abs(d__1)) + (
			d__2 = h__[ii + ii * h_dim1], abs(d__2))) * 1.5;
		wi[ii] = 0.;
/* L80: */
	    }
	} else {

/*           Use eigenvalues of trailing submatrix of order NS as 
shifts. */

	    dlacpy_("Full", &ns, &ns, &h__[i__ - ns + 1 + (i__ - ns + 1) * 
		    h_dim1], ldh, s, &c__15, 4L);
	    dlahqr_(&c_false, &c_false, &ns, &c__1, &ns, s, &c__15, &wr[i__ - 
		    ns + 1], &wi[i__ - ns + 1], &c__1, &ns, &z__[z_offset], 
		    ldz, &ierr);
	    if (ierr > 0) {

/*              If DLAHQR failed to compute all NS eigenvalues
, use the */
/*              unconverged diagonal elements as the remaining
 shifts. */

		i__2 = ierr;
		for (ii = 1; ii <= i__2; ++ii) {
		    wr[i__ - ns + ii] = s[ii + ii * 15 - 16];
		    wi[i__ - ns + ii] = 0.;
/* L90: */
		}
	    }
	}

/*        Form the first column of (G-w(1)) (G-w(2)) . . . (G-w(ns)) 
*/
/*        where G is the Hessenberg submatrix H(L:I,L:I) and w is */
/*        the vector of shifts (stored in WR and WI). The result is */
/*        stored in the local array V. */

	v[0] = 1.;
	i__2 = ns + 1;
	for (ii = 2; ii <= i__2; ++ii) {
	    v[ii - 1] = 0.;
/* L100: */
	}
	nv = 1;
	i__2 = i__;
	for (j = i__ - ns + 1; j <= i__2; ++j) {
	    if (wi[j] >= 0.) {
		if (wi[j] == 0.) {

/*                 real shift */

		    i__4 = nv + 1;
		    dcopy_(&i__4, v, &c__1, vv, &c__1);
		    i__4 = nv + 1;
		    d__1 = -wr[j];
		    dgemv_("No transpose", &i__4, &nv, &c_b10, &h__[l + l * 
			    h_dim1], ldh, vv, &c__1, &d__1, v, &c__1, 12L);
		    ++nv;
		} else if (wi[j] > 0.) {

/*                 complex conjugate pair of shifts */

		    i__4 = nv + 1;
		    dcopy_(&i__4, v, &c__1, vv, &c__1);
		    i__4 = nv + 1;
		    d__1 = wr[j] * -2.;
		    dgemv_("No transpose", &i__4, &nv, &c_b10, &h__[l + l * 
			    h_dim1], ldh, v, &c__1, &d__1, vv, &c__1, 12L);
		    i__4 = nv + 1;
		    itemp = idamax_(&i__4, vv, &c__1);
/* Computing MAX */
		    d__2 = (d__1 = vv[itemp - 1], abs(d__1));
		    temp = 1. / max(d__2,smlnum);
		    i__4 = nv + 1;
		    dscal_(&i__4, &temp, vv, &c__1);
		    absw = dlapy2_(&wr[j], &wi[j]);
		    temp = temp * absw * absw;
		    i__4 = nv + 2;
		    i__5 = nv + 1;
		    dgemv_("No transpose", &i__4, &i__5, &c_b10, &h__[l + l * 
			    h_dim1], ldh, vv, &c__1, &temp, v, &c__1, 12L);
		    nv += 2;
		}

/*              Scale V(1:NV) so that max(abs(V(i))) = 1. If V
 is zero, */
/*              reset it to the unit vector. */

		itemp = idamax_(&nv, v, &c__1);
		temp = (d__1 = v[itemp - 1], abs(d__1));
		if (temp == 0.) {
		    v[0] = 1.;
		    i__4 = nv;
		    for (ii = 2; ii <= i__4; ++ii) {
			v[ii - 1] = 0.;
/* L110: */
		    }
		} else {
		    temp = max(temp,smlnum);
		    d__1 = 1. / temp;
		    dscal_(&nv, &d__1, v, &c__1);
		}
	    }
/* L120: */
	}

/*        Multiple-shift QR step */

	i__2 = i__ - 1;
	for (k = l; k <= i__2; ++k) {

/*           The first iteration of this loop determines a reflect
ion G */
/*           from the vector V and applies it from left and right 
to H, */
/*           thus creating a nonzero bulge below the subdiagonal. 
*/

/*           Each subsequent iteration determines a reflection G t
o */
/*           restore the Hessenberg form in the (K-1)th column, an
d thus */
/*           chases the bulge one step toward the bottom of the ac
tive */
/*           submatrix. NR is the order of G. */

/* Computing MIN */
	    i__4 = ns + 1, i__5 = i__ - k + 1;
	    nr = min(i__4,i__5);
	    if (k > l) {
		dcopy_(&nr, &h__[k + (k - 1) * h_dim1], &c__1, v, &c__1);
	    }
	    dlarfg_(&nr, v, &v[1], &c__1, &tau);
	    if (k > l) {
		h__[k + (k - 1) * h_dim1] = v[0];
		i__4 = i__;
		for (ii = k + 1; ii <= i__4; ++ii) {
		    h__[ii + (k - 1) * h_dim1] = 0.;
/* L130: */
		}
	    }
	    v[0] = 1.;

/*           Apply G from the left to transform the rows of the ma
trix in */
/*           columns K to I2. */

	    i__4 = i2 - k + 1;
	    dlarfx_("Left", &nr, &i__4, v, &tau, &h__[k + k * h_dim1], ldh, &
		    work[1], 4L);

/*           Apply G from the right to transform the columns of th
e */
/*           matrix in rows I1 to min(K+NR,I). */

/* Computing MIN */
	    i__5 = k + nr;
	    i__4 = min(i__5,i__) - i1 + 1;
	    dlarfx_("Right", &i__4, &nr, v, &tau, &h__[i1 + k * h_dim1], ldh, 
		    &work[1], 5L);

	    if (wantz) {

/*              Accumulate transformations in the matrix Z */

		dlarfx_("Right", &nh, &nr, v, &tau, &z__[*ilo + k * z_dim1], 
			ldz, &work[1], 5L);
	    }
/* L140: */
	}

/* L150: */
    }

/*     Failure to converge in remaining number of iterations */

    *info = i__;
    return 0;

L160:

/*     A submatrix of order <= MAXB in rows and columns L to I has split 
*/
/*     off. Use the double-shift QR algorithm to handle it. */

    dlahqr_(&wantt, &wantz, n, &l, &i__, &h__[h_offset], ldh, &wr[1], &wi[1], 
	    ilo, ihi, &z__[z_offset], ldz, info);
    if (*info > 0) {
	return 0;
    }

/*     Decrement number of remaining iterations, and return to start of */
/*     the main loop with a new value of I. */

    itn -= its;
    i__ = l - 1;
    goto L50;

L170:
    work[1] = (doublereal) max(1,*n);
    return 0;

/*     End of DHSEQR */

} /* dhseqr_ */

