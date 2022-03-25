//Gavin Conant 5/26/1999
//Taken from Numerical Recipes in C, 2nd Edition
//Press, W.H., Teukolsky, S.A., Vetterling, W.T., Flannery, B.P.
//Cambridge University Press, 1994.

#include <math.h>
#ifndef ___NRUTILS_H___
#define ___NRUTILS_H___

static double sqrarg;

#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a)) 

static double maxarg1, maxarg2;
#define FMAX(a,b) (maxarg1=(a), maxarg2=(b), maxarg1>maxarg2 ? (maxarg1) : (maxarg2))

#ifdef __cplusplus

extern "C" double *NRvector(long nl, long nh);
extern "C" void free_NRvector(double *v, long nl, long nh);
extern "C" double **matrix (long nrl, long nrh, long ncl, long nch);
extern "C" void free_matrix(double **m, long nrl, long nrh, long ncl, long nch);

#else
double *NRvector(long nl, long nh);
void free_NRvector(double *v, long nl, long nh);
double **matrix (long nrl, long nrh, long ncl, long nch);
void free_matrix(double **m, long nrl, long nrh, long ncl, long nch);
#endif
#endif
