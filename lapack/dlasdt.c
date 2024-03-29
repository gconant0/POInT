/* dlasdt.f -- translated by f2c (version 20000121).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int dlasdt_(n, lvl, nd, inode, ndiml, ndimr, msub)
integer *n, *lvl, *nd, *inode, *ndiml, *ndimr, *msub;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double log();

    /* Local variables */
    static integer maxn;
    static doublereal temp;
    static integer nlvl, llst, i__, ncrnt, il, ir;


/*  -- LAPACK auxiliary routine (version 3.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     June 30, 1999 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLASDT creates a tree of subproblems for bidiagonal divide and */
/*  conquer. */

/*  Arguments */
/*  ========= */

/*   N      (input) INTEGER */
/*          On entry, the number of diagonal elements of the */
/*          bidiagonal matrix. */

/*   LVL    (output) INTEGER */
/*          On exit, the number of levels on the computation tree. */

/*   ND     (output) INTEGER */
/*          On exit, the number of nodes on the tree. */

/*   INODE  (output) INTEGER array, dimension ( N ) */
/*          On exit, centers of subproblems. */

/*   NDIML  (output) INTEGER array, dimension ( N ) */
/*          On exit, row dimensions of left children. */

/*   NDIMR  (output) INTEGER array, dimension ( N ) */
/*          On exit, row dimensions of right children. */

/*   MSUB   (input) INTEGER. */
/*          On entry, the maximum row dimension each subproblem at the */
/*          bottom of the tree can be of. */

/*  Further Details */
/*  =============== */

/*  Based on contributions by */
/*     Ming Gu and Huan Ren, Computer Science Division, University of */
/*     California at Berkeley, USA */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Find the number of levels on the tree. */

    /* Parameter adjustments */
    --ndimr;
    --ndiml;
    --inode;

    /* Function Body */
    maxn = max(1,*n);
    temp = log((doublereal) maxn / (doublereal) (*msub + 1)) / log(2.);
    *lvl = (integer) temp + 1;

    i__ = *n / 2;
    inode[1] = i__ + 1;
    ndiml[1] = i__;
    ndimr[1] = *n - i__ - 1;
    il = 0;
    ir = 1;
    llst = 1;
    i__1 = *lvl - 1;
    for (nlvl = 1; nlvl <= i__1; ++nlvl) {

/*        Constructing the tree at (NLVL+1)-st level. The number of */
/*        nodes created on this level is LLST * 2. */

	i__2 = llst - 1;
	for (i__ = 0; i__ <= i__2; ++i__) {
	    il += 2;
	    ir += 2;
	    ncrnt = llst + i__;
	    ndiml[il] = ndiml[ncrnt] / 2;
	    ndimr[il] = ndiml[ncrnt] - ndiml[il] - 1;
	    inode[il] = inode[ncrnt] - ndimr[il] - 1;
	    ndiml[ir] = ndimr[ncrnt] / 2;
	    ndimr[ir] = ndimr[ncrnt] - ndiml[ir] - 1;
	    inode[ir] = inode[ncrnt] + ndiml[ir] + 1;
/* L10: */
	}
	llst <<= 1;
/* L20: */
    }
    *nd = (llst << 1) - 1;

    return 0;

/*     End of DLASDT */

} /* dlasdt_ */

