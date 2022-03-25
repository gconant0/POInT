//Numerical recipes library functions
//Gavin Conant 6/2/1999

#include <stdlib.h>
#include "nrutil.h"




#define NR_END 1
#define FREE_ARG char*

double *NRvector(long nl, long nh)
{
  double *v;
  v=(double*)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));

  return(v-nl+NR_END);
}

double **matrix (long nrl, long nrh, long ncl, long nch)
{
  long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
  double **m;

  m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
  
  m[nrl]=(double*)malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i=nrl+1; i<=nrh; i++)
    m[i]=m[i-1]+ncol;
  
  return(m);
}

void free_NRvector(double *v, long nl, long nh)
{
  free((FREE_ARG) (v+nl-NR_END));
}

void free_matrix(double **m, long nrl, long nrh, long ncl, long nch)
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}




