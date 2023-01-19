//Gavin Conant  5/26/1999
//Defines maximization routines
//requires powell.h

#ifndef ___POWELL_H___
#include "powell.h"
#endif


#include <iostream>
#include <iomanip>

//#define SHOW_PROGRESS


///using namespace::std;
#ifdef _OMP_MP_VERSION_
int ncom=0;
double *pcom, *xicom;
Like_model *ext_mod;
//#elif defined(LINUX_BUILD)
#else
int ncom;
double *pcom, *xicom;
Like_model *ext_mod;
//#else 
//extern int ncom;
//extern double *pcom, *xicom;
//extern Like_model *ext_mod;
#endif


double Powell::Init_min(Like_model *cmodel, Exchange *cexchange, BOOL brnopt)
{
  double newlnL=0.0, origlnL, *powell_params, **xdirs, toler=2*tolerance;
  int iters=0,i,j;

  curr_model=cmodel;
  curr_exchange=cexchange;
  opt_branches=brnopt;


  set_extern_model(curr_model);

  if (curr_exchange->get_num_params()!=0)
    {
      powell_params= NRvector(1, curr_exchange->get_num_params());
      xdirs = matrix(1, curr_exchange->get_num_params(), 1, curr_exchange->get_num_params());
      
      for (i=1; i<=curr_exchange->get_num_params(); i++)
	for (j=1; j<=curr_exchange->get_num_params(); j++)
	  {
	    if (i+j==curr_exchange->get_num_params()+1)
	      xdirs[i][j]=1.0;
	    else
	      xdirs[i][j]=0.0;
	    
	  }

      double *swtch = NRvector(1,curr_exchange->get_num_params());
  
      for (i=1; i<=curr_exchange->get_num_params(); i++)
	{
	  swtch[i]=xdirs[i][1];
	  xdirs[i][1]=xdirs[i][2];
	  xdirs[i][2]=swtch[i];
	} 


      if (brnopt==TRUE)
		  for(i=0; i<2; i++)
		curr_model->newton_branch_opt(0.01);
	  //std::cout<<"Setting params\n";
      //Receives the parameters from the like model with the index 
      //offset by 1  (Starts at 1, not 0)
      curr_model->send_params(powell_params, 1); 
      curr_model->recalculate_transprobs();
      origlnL=curr_model->find_appropriate_ln_like();
      
      
#if defined (SHOW_PROGRESS)
      stdstd::cout<<"Starting parameters: "<<powell_params[1];
      for (i=2; i<=curr_exchange->get_num_params(); i++)
	std::cout<<", "<<powell_params[i];
      std::cout<<endl;

      std::cout<<"Tolerance: "<<tolerance<<endl;
      std::cout<<"Initial ln Likelihood:    "<<setw(14)<<setprecision(9)<<origlnL<<endl<<endl;;
#endif    
		if (swap_brn==TRUE) curr_model->swap_brn_opt=TRUE;	
      Powell::powell_fn(powell_params, xdirs, curr_exchange->get_num_params(), 
		toler, &iters, &newlnL, curr_model); 
      if (swap_brn==TRUE) curr_model->swap_brn_opt=FALSE;	
      if (brnopt==TRUE) {
#ifdef SHOW_PROGRESS		  
	std::cout<<"Branch opt\n";
#endif
	for(i=0; i<3; i++)
	  curr_model->newton_branch_opt();
      }
      newlnL=Powell::curr_model->find_appropriate_ln_like();
      //Powell minimizes the -lnL, so we must reset to lnL

#ifdef SHOW_PROGRESS        
      Powell::curr_model->describe_results();
      

      std::cout<<"\nFinal ln likelihood: "<<setw(12)<<setprecision(7)<<newlnL<<endl<<endl;
#endif      
      
      free_NRvector(powell_params, 1, curr_exchange->get_num_params());
      if (curr_exchange->get_num_params()>1)
	free_matrix(xdirs,1, curr_exchange->get_num_params(),1, Powell::curr_exchange->get_num_params());
     
      return(newlnL);
    }
  
}  //End Init_min




void Powell::set_extern_model(Like_model *cmodel)
{
  ext_mod=cmodel;
}




void Powell::powell_fn (double p[], double **xi, int n, double ftol, int *iter,
		    double *fret, Like_model *themod)
{
	int i, ibig, j;
	double del, fp, fptt, t, *pt, *ptt, *xit;
  
	pt=NRvector(1,n);
	ptt=NRvector(1,n);
	xit=NRvector(1,n);
	*fret=themod->min_lnL(p, 1);
 

	for (j=1; j<=n; j++) 
		pt[j]=p[j];
 
	for (*iter=1;;++(*iter))
    {
		fp=(*fret);
		ibig=0;
		del=0.0;
 
		for (i=1; i<=n; i++)
		{
			for (j=1; j<=n; j++)
				xit[j]=xi[j][i];
			fptt=(*fret);
#ifdef CHECK_ARRAY
            std::cout<<"P array pre linmin:\n";
            for (j=1; j<=n; j++) std::cout<<p[j]<<"\t";
            std::cout<<endl;
#endif
			linmin(p, xit, n, fret, themod);
#ifdef CHECK_ARRAY
            std::cout<<"P array post linmin:\n";
            for (j=1; j<=n; j++) std::cout<<p[j]<<"\t";
            std::cout<<endl;
#endif
			if (fabs(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret)))
			{
				del=fabs(fptt-(*fret));
				ibig=i;
			}
		}
     
		if (2.0*fabs(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret)))
		{
			free_NRvector(xit,1,n);
			free_NRvector(ptt,1,n);
			free_NRvector(pt,1,n);
			return;
		}

		if (*iter==ITMAX)
			cerr<<"Powell::powell-fn exceeded maximum iterations\n";

		for (j=1; j<=n; j++)
		{
			ptt[j]=2.0*p[j]-pt[j];
			xit[j]=p[j]-pt[j];
			pt[j]=p[j];
		}
#ifdef CHECK_ARRAY
        std::cout<<"P array rerun prior to last lnL:\n";
        for (j=1; j<=n; j++) std::cout<<p[j]<<"\t";
        
        std::cout<<endl;
        std::cout<<"PrelnL: "<<themod->min_lnL(p, 1)<<endl;
        
#endif
        
#ifdef CHECK_ARRAY
        std::cout<<"Ptt array last lnL:\n";
        for (j=1; j<=n; j++) std::cout<<ptt[j]<<"\t";
        std::cout<<endl;
#endif
		fptt=themod->min_lnL(ptt, 1);
#ifdef CHECK_ARRAY
        std::cout<<"P array rerun after last lnL:\n";
        for (j=1; j<=n; j++) std::cout<<p[j]<<"\t";
        
        std::cout<<endl;
        std::cout<<"PostlnL: "<<themod->min_lnL(p, 1)<<endl;
        
#endif


		if (fptt < fp)
		{
			t=2.0*(fp-2.0*(*fret)+fptt)*SQR(fp-(*fret)-del)-del*SQR(fp-fptt);

			if (t< 0.0)
			{
				linmin(p, xit, n ,fret, themod);
				for (j=1; j<=n; j++)
				{
					xi[j][ibig]=xi[j][n];
					xi[j][n]=xit[j];
				}
			}

		}

    }

}  //End powell_fn




void Powell::linmin( double p[], double xi[], int n, double *fret, 
	       Like_model *themod)
{
  
  int j;
  double xx, xmin, fx, fb, fa, bx, ax;
  double (*thefunc)(double);

  ncom=n;
  pcom=NRvector(1,n);
  xicom=NRvector(1,n);

  thefunc=&f1dim;
  
  for (j=1; j<=n; j++)
    {
      pcom[j]=p[j];
      xicom[j]=xi[j];
    }
 
    ax=0.0;
    xx=1.0;
  
    mnbrak(&ax, &xx, &bx, &fa, &fx, &fb, thefunc);

#if defined (SHOW_PROGRESS) 
    std::cout<<"Min braketed by: "<<ax<<","<<setw(10)<<setprecision(7)<<fa<<" "<<bx<<","
      <<setw(10)<<setprecision(7)<<fb<<"( total params: "<<n<<")"<<endl;
#endif

    *fret=brent(ax,xx,bx, thefunc, tolerance, &xmin);

    for (j=1; j<=n; j++)
      {
	xi[j] *= xmin;
	p[j] += xi[j];
      }

    curr_model->just_found_opt=TRUE;
#if defined (SHOW_PROGRESS)
    std::cout<<"Min found at: "<<xmin<<","<<setw(10)<<setprecision(7)<<*fret<<endl<<"\t";
#endif
    //Powell::curr_model->optimize_branches(); 
    if (opt_branches == TRUE) 
      *fret=curr_model->newton_branch_opt(0.001);
    free_NRvector(xicom, 1, n);
    free_NRvector(pcom, 1, n);
}  //End Powell::linmin




void mnbrak(double *ax, double *bx, double *cx, double *fa, 
		    double *fb, double *fc, double (*func)(double))
  //This function is taken from Numerical Recipes in C--it brackets the minimum 
  //of a function using parabolic interpolation with an increasing step size at 
  //each iteration
{
  double ulim, u,r,q,fu,dum;
    //std::cout<<"Calling min_lnL from mnbrakA\n";
  *fa=(*func)(*ax);
    //std::cout<<"Calling min_lnL from mnbrakB\n";
  *fb=(*func)(*bx);
  if (*fb > *fa)
    {
      SHFT(dum, *ax, *bx, dum);
      SHFT(dum, *fa, *fb, dum);  //Note:book gives dum,fb,fa,dum
    }

 
  *cx=(*bx)+GOLD*(*bx-*ax);
    //std::cout<<"Calling min_lnL from mnbrakC\n";
  *fc=(*func)(*cx);

  while (*fb> *fc)
    {
     r=(*bx-*ax)*(*fb-*fc);
     q=(*bx-*cx)*(*fb-*fa);
     u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
       (2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
     ulim=(*bx)+GLIMIT*(*cx-*bx);
	      
     if((*bx-u)*(u-*cx)>0.0)
      {
		fu=(*func)(u);
		if (fu<*fc)                   //A min was found between b and c
		{
			*ax=(*bx);
			*bx=u;
			*fa=(*fb);
			*fb=fu;
			return;
		}
		else if (fu>*fb)             //A min was found between a and b
		{
			*cx=u;
			*fc=fu;
			return;
		} 
		u=(*cx)+GOLD*(*cx-*bx);      //No min was found, so step-size is increased
		fu=(*func)(u);
       }
    
     else if ((*cx-u)*(u-ulim) >0.0)   
       {
	 fu=(*func)(u);
	 if (fu <*fc)
	   {
	     SHFT(*bx, *cx, u, *cx+GOLD*(*cx-*bx));
	     SHFT(*fb, *fc, fu, (*func)(u));
	   }
       }
     else if ((u-ulim)*(ulim-*cx)>=0.0)
       {
	 u=ulim;
	 fu=(*func)(u);
       }
     else 
       {
	 u=(*cx)+GOLD*(*cx-*bx);
	 fu=(*func)(u);
       }
     SHFT(*ax, *bx, *cx, u);
     SHFT(*fa,*fb, *fc, fu);
    } 
}  //End Powell::mnbrak




double brent (double ax, double bx, double cx, double (*f)(double),  
		      double tol, double *xmin)
{
  int iter;
  double a,b,d, etemp, fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  double e=0.0;

  a=(ax < cx ? ax : cx);
  b=(ax > cx ? ax : cx);
  x=w=v=bx;
  fw=fv=fx=(*f)(x);
  
 
  for (iter=1; iter<=ITMAX; iter++)
    {
      xm=0.5*(a+b);
      tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
      if (fabs(x-xm) <= (tol2-0.5*(b-a)))
	{
	  *xmin=x;
	  return(fx);
	}
      if (fabs(e)> tol1)
	{
	  r=(x-w)*(fx-fv);
	  q=(x-v)*(fx-fw);
	  p=(x-v)*q-(x-w)*r;
	  q=2.0*(q-r);

	  if (q>0.0) 
	    p=-p;

	  q=fabs(q);
	  etemp=e;
	  e=d;
	  if (fabs(p) >= fabs(0.5*q*etemp) || p<=q*(a-x) || p>=q*(b-x))
	    d=CGOLD*(e=(x>=xm ? a-x :b-x));
	  else 
	    {
	      d=p/q;
	      u=x+d;
	      if (u-a < tol2 || b-u < tol2)
		d=SIGN(tol1, xm-x);
	    }
	}
      else      //Else for if (fabs(e)...
	{
	  d=CGOLD*(e=(x>=xm ? a-x: b-x));
	}
      
      u=(fabs(d) >-tol1 ? x+d : x+SIGN(tol1, d));
      fu=(*f)(u);

      if (fu <= fx)
	{
	  if (u >= x)
	    a=x;
	  else
	    b=x;
	  SHFT(v,w,x,u);
	  SHFT(fv,fw,fx,fu);
	}
      else
	{
	  if (u<x) 
	    a=u;
	  else
	    b=u;
	  if (fu<= fw || w==x)
	    {
	      v=w;
	      w=u;
	      fv=fw;
	      fw=fu;
	    }
	  else if (fu<=fv || v==x || v==w)
	    {
	      v=u;
	      fv=fu;
	    }
	}
    }
  cerr<<"Too many iterations in brent"<<endl;
  *xmin=x;
  return (fx);
}  //End brent




double f1dim (double x)
{
    //An external function that allows brent to
    //interact with the like model object when 
    //performing a multi-dimentional optimization
    int j;
    double f, *xt;

    xt=NRvector(1, ncom);

    for (j=1; j<=ncom; j++)
        xt[j]=pcom[j]+(x*xicom[j]);


    f=ext_mod->min_lnL(xt, 1);

    free_NRvector(xt, 1, ncom);
    return(f);
}

 

double brn_min(double x)
{
  //An external function that allows brent
  //to interact with the Like_model object
  //when performing a branch length optimization
  //(You can't pass a pointer to a class function
  //as a pointer to a function
  return(-ext_mod->optimize_a_branch(x));
}




double single_min (double trs_trv)
{
 //An external function that allows brent
  //to interact with the Like_model object
  //when performing only a trs/trv optimization
  //(You can't pass a pointer to a class function
  //as a pointer to a function
  return(ext_mod->optimize_single_param(trs_trv));
}


double snp_brn_min (double brnlen)
{
	//An external function that allows brent
	//to interact with the SNP_model object
	//when performing branch length optimization
	//(You can't pass a pointer to a class function
	//as a pointer to a function
	return(ext_mod-> optimize_a_branch(brnlen));
}
