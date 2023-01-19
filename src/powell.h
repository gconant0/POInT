//Gavin Conant 5/26/1999
//Defines maximization routines for the prosrch program.
#include "maxlike.h"
#include "codon_like.h"
#include "nucleotide_like.h"
#include "other_like.h"
#include "structs.h"
#include "nrutil.h"


#ifndef ___POWELL_H___
#define ___POWELL_H___


#define GOLD 1.618034
#define GLIMIT 10.0
#ifndef TINY
#define TINY 1.0e-20
#endif
#define SHFT(a,b,c,d) (a)=(b); (b)=(c);(c)=(d);
#define ITMAX 200
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
//#define TOL tolerance;

class Powell {
 public:
	BOOL swap_brn;	
	Powell() {tolerance=1.0e-5; swap_brn=FALSE;};
	double Init_min(Like_model *cmodel, Exchange *cexchange, BOOL brnopt);
	void set_extern_model(Like_model *cmodel);
	void set_tolerance (double tol) {tolerance=tol;};
	double get_tolerance () {return(tolerance);};
protected:
	double tolerance;
	BOOL opt_branches;
	Exchange *curr_exchange;
	Like_model *curr_model;


  void powell_fn (double p[], double **xi, int n, double ftol, int *iter,
		  double *fret, Like_model *themod);
  
  void linmin( double p[], double xi[], int n, double *fret, 
	       Like_model *themod);
};


void mnbrak(double *ax, double *bx, double *cx, double *fa,
	      double *fb, double *fc, double (*func)(double));

double brent (double ax, double bx, double cx, double (*f)(double), 
		double tol, double *xmin);

double f1dim (double x);

double brn_min (double x); 

double single_min (double trs_trv);
double snp_brn_min (double brnlen);

#endif
