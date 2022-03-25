//Copyright 1999-2002 Gavin Conant

#include <iostream>
#include <math.h>
#include <iomanip>
#include "maxlike.h"
#include "powell.h"

#ifdef _OPEN_MP_VERSION_
#include "omp.h"
#endif

using namespace::std;

//#define DO_TIME
//#define PRINT_WARNINGS

#ifdef DO_TIME
#include <sys/time.h>
#include <sys/resource.h>

long pre, post;
unsigned long RunTime()
{
    struct rusage r;
    long sec, usec;

    getrusage(RUSAGE_SELF, &r);   // constant defined in include files 
    sec = r.ru_utime.tv_sec;      // whole seconds
    usec = r.ru_utime.tv_usec;    // microseconds

    return((sec * 1000) + (usec / 1000));

} // RunTime
#endif


//#define DEBUG
#define SHOW_PROGRESS
#define MAX_PARAM 100

Like_model::Like_model()
{
}

void Like_model::assemble (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree)
  //Generic function called by all Like_model class inheritors to see up needed arrays etc.
  //Arrays are allocated to do branch lenght optimization, as well as to do numerical optimization
  //(each model has a function that expresses the number of free parameters)
  //The branch lenght optimization expects a tree where the left(right?) child of the
  //root has zero length--this is set here.
{
  int i,j,m;

  curr_exchange=cexchange;
  curr_tree=ctree;
  curr_data=cdata;
  last_codonmat=-1;
	
	//cout<<"Num sites: "<<(*curr_data)[0].Sequence_size()<<endl;
 
  prop_index=0;
  prop_pns_index=0;
  prop_rate_index=0;
  matrix_index=0;
  num_aa_props_per_rate=0;
  swap_brn_opt=FALSE;
  last_brn_remain=0;

  //Allows us not to optimize certain branches if exchange.zero_brn_fixed() ==TRUE
  brn_index=0;

  //Initializes our rates
  for (i=0; i<curr_exchange->get_num_rates(); i++)
    curr_exchange->set_rate(i, 1.0);
 

  just_found_opt=FALSE;
 
  if (curr_exchange->branch_basefreqs()==FALSE)
    {
	  if (curr_exchange->are_optimizing_rate_props() == TRUE)
		lookup_p_tree=new Constrain_Param_Lookup(curr_exchange->get_num_rates());
      
	  num_params_model();
      
      //Holds the parameters the min routine is allowed to change and the identity
      //of those parameters
      params= new double[curr_exchange->get_num_params()];
      param_types= new PARAM_TYPE[curr_exchange->get_num_params()];

      //The next three variables hold ratios of base frequencies that can be altered indepedantly by
      //Powell's routine 
      pur_pyr_split=(get_basefreq(0,(*curr_tree)[0])+get_basefreq(2,(*curr_tree)[0]))/
		(get_basefreq(1,(*curr_tree)[0])+get_basefreq(3,(*curr_tree)[0]));
      
      a_g_split=get_basefreq(0,(*curr_tree)[0])/get_basefreq(2,(*curr_tree)[0]);
      c_t_split=get_basefreq(1,(*curr_tree)[0])/get_basefreq(3,(*curr_tree)[0]);
    }

  if (curr_exchange->get_trs_trv() == 1)
    curr_exchange->set_trs_trv(curr_exchange->get_obs_trs_trv()/
			       ((get_basefreq(0,(*curr_tree)[0])*get_basefreq(2,(*curr_tree)[0])+
				 get_basefreq(1,(*curr_tree)[0])*get_basefreq(3,(*curr_tree)[0]))/
				((get_basefreq(0,(*curr_tree)[0])+get_basefreq(2,(*curr_tree)[0]))*
				 (get_basefreq(1,(*curr_tree)[0])+get_basefreq(3,(*curr_tree)[0])))));
  
  start_swap=0;
  brn_order=new int [curr_exchange->get_num_branches()];
  for (i=0; i<curr_exchange->get_num_branches(); i++) {
    (*curr_tree)[i]->set_brnlen(find_ut((*curr_tree)[i]));
	//cout<<"Brlen: "<<(*curr_tree)[i]->get_brnlen()<<endl;
    brn_order[i] = i;
  }

  if ((curr_tree->rooted_tree() == FALSE) && (curr_tree->find_null_branch_id()->get_brnlen()!=0))
    {
      curr_tree->find_null_branch_id()->get_sibling()->
	set_brnlen(curr_tree->find_null_branch_id()->get_brnlen()+
		  curr_tree->find_null_branch_id()->get_sibling()->get_brnlen());
		    
      curr_tree->find_null_branch_id()->set_brnlen(0.0);
    }

  
  //The array that holds the likelihoods for each site/codon
  site_lnLs= new double[curr_exchange->get_num_localities()];

  brncondprobs1= new double ** [curr_exchange->get_num_localities()];
  brncondprobs2= new double ** [curr_exchange->get_num_localities()];
    
  for (i=0; i<curr_exchange->get_num_localities(); i++)
    {
		brncondprobs1[i]= new double * [curr_exchange->get_num_rates()];
		brncondprobs2[i]= new double * [curr_exchange->get_num_rates()];
		for(j=0; j<curr_exchange->get_num_rates(); j++) {
			brncondprobs1[i][j]= new double [curr_exchange->get_condlike_size()];
			brncondprobs2[i][j]= new double [curr_exchange->get_condlike_size()];
	  }
    }
  
  if (curr_exchange->full_conprobs() == TRUE) {
		//We need tree depth -1 arrays, plus one for the branch and 2 for the scratch arrays
		save_conprobs_size=curr_tree->find_max_tree_depth()+3;
		save_conprobs=new double **[save_conprobs_size];
		for(i=0; i<save_conprobs_size; i++) {
			save_conprobs[i]=new double *[curr_exchange->get_num_localities()];
			for(j=0; j<curr_exchange->get_num_localities(); j++)
				save_conprobs[i][j]=new double [curr_exchange->get_condlike_size()];
		}
  }
  else save_conprobs=0;

  //Sets up the differences between codon model and site model
  initialize_arrays();
 

  //Parameters available to the minimization function are intialized by this function. 
  if (curr_exchange->branch_basefreqs()==FALSE) 
    intialize_parameters(params, param_types);
   
  recalculate_transprobs();
  //Print model name
//  print_model();
  
} //End Like_model::Like_model





void Like_model::initialize_arrays()
  //Virtual function
{
  cerr<<"Wrong initialize_arrays\n";
}




void Like_model::recalculate_transprobs()
  //Using the private function calc_transprobs, calculates the transition probs over the
  //whole tree.  calc_transprobs is a polymorphic (virtual) function determined
  //by the model of evolution chosen--this version is called for nucleotide models--
  //codon models have ther own version (codon_like.cpp)

{
	int i,j,k;

	if (curr_exchange->modeling_codons()==TRUE)
	{
		if( (curr_exchange->branch_basefreqs()==FALSE))
		{
			for(k=0;k<curr_exchange->get_num_rates(); k++)  {
				for(j=0; j<curr_exchange->get_num_nonsyn_patterns(); j++) {
					i=0;
					while((*curr_tree)[i]->get_nonsyn_pattern() != j) i++;
						calc_codon_matrices(k, (*curr_tree)[i]);
					for (i=0; i<curr_exchange->get_num_branches(); i++) {
						if ((*curr_tree)[i]->get_nonsyn_pattern() == j) 
							calc_transprobs((*curr_tree)[i], k);
					}
			
				}
			}
		}
     
		else
		{
			for (i=0; i<curr_exchange->get_num_branches(); i++)
			{ 
				calc_codon_matrices(0, (*curr_tree)[i]);
				for(j=0;j<curr_exchange->get_num_rates(); j++) 
					calc_transprobs((*curr_tree)[i], j);
			}
		}
	}          

	else
		for (i=0; i<curr_exchange->get_num_branches(); i++)
			for(j=0;j<curr_exchange->get_num_rates(); j++)
				calc_transprobs((*curr_tree)[i], j);
 

 
}  //End Like_model::recalculate_transprobs
 

void Like_model::reinit_params()
{
  intialize_parameters(params, param_types);
}

void Like_model::list_transprobs(Branch *taxa)
  //Lists the transition probability matrix for the requested branch
{
  int i,j;
 
  cout<<"ut: "<<taxa->get_brnlen()<<endl;  
  for (i=0; i<curr_exchange->get_condlike_size(); i++)
    {
      for (j=0; j<curr_exchange->get_condlike_size(); j++)
	cout<<setprecision(2)<<setw(4)<<taxa->get_trpb(0, i, j)<<" ";
	
      cout<<endl;
    }
 
} //End list_transprobs




void Like_model::list_rates()
  //Lists the current value of the relative rate in each rate catagory
{
  int i;
  
  for(i=0; i<curr_exchange->get_num_rates(); i++)
    cout<<setw(5)<<curr_exchange->get_rate(i)<<" ";
  cout<<endl;

}  //End Like_model::list_rates




void Like_model::send_params(double remote_params[], int offset)
  //Copys the parameter array.
  //(Note: the powell parameter array is numbered from 1 instead of 0)
  //and hence the params array from Powell and the param_types array are 
  //index 1 differently from each other
{
  int i;

  for (i=0; i< curr_exchange->get_num_params(); i++)
      remote_params[i+offset]=params[i];
    
}  //End Like_model::send_params



PARAM_TYPE Like_model::get_param_type(int param_num)
{
  return(param_types[param_num]);
}


double Like_model::find_appropriate_ln_like() 
//An aliasing function that allows the min_lnL function to handle 
//several different possible setups for site-wise rates
{
	if(curr_exchange->get_num_rates() == 1) {
		if (curr_exchange->get_num_taxa() > 2)
			return(find_ln_like_single_rate());
		else {
			//cout<<"OLD: "<<find_ln_like_single_rate()<<" New: "<<find_two_taxa_ln_like()<<endl;
			return(find_two_taxa_ln_like());
			//return(find_ln_like_single_rate());
		}
	}
	else
		//Note--multiple things may go here eventually: now the choice is a
		//multi-rate model with fixed probs per site for each rate
		return(find_ln_like_fixed_site_rates());
}


double Like_model::find_two_taxa_ln_like() 
{
	int i, j;
	double ln_like=0, condprob;
	Branch *lsib, *rsib;
	
	lsib=(*curr_tree)[0];
	rsib=(*curr_tree)[1];
	
	for(i=0; i<curr_exchange->get_num_localities(); i++) {
		if(( (*curr_data)[lsib->get_taxa_id()][i] != 64) && ( (*curr_data)[rsib->get_taxa_id()][i] != 64))
		   ln_like += log(rsib->get_trpb(0, (*curr_data)[lsib->get_taxa_id()][i] , (*curr_data)[rsib->get_taxa_id()][i])*root_freq((*curr_data)[lsib->get_taxa_id()][i]));
		else {
            if (!(((*curr_data)[lsib->get_taxa_id()][i]==64) && ((*curr_data)[rsib->get_taxa_id()][i]==64))) {
                if ((*curr_data)[lsib->get_taxa_id()][i] == 64) ln_like+=log(root_freq((*curr_data)[rsib->get_taxa_id()][i]));
                else ln_like+=log(root_freq((*curr_data)[lsib->get_taxa_id()][i]));
            }
		}
		//cout<<i<<": "<<ln_like<<" : "<<(*curr_data)[lsib->get_taxa_id()][i]<<" : "<<(*curr_data)[lsib->get_taxa_id()][i]<<endl;
	//		condprob=0;
	//		for(j=0; j<curr_exchange->get_condlike_size(); j++) {
	//			condprob +=root_freq(j);
	//		}
	//		ln_like+=log(condprob);
	//	}
					
		
	}
	 	
	
	return(ln_like);	
}

double Like_model::find_ln_like_single_rate()
  //This function finds the lnL using transprob matrix
  //assuming the whole sequence evolves at a single rate (#0)
{
	int i, j; 
	long double prob, pretime, posttime;
	double final_lnL=0, final_lnL_s=0;
#ifdef DO_TIME
	pretime=RunTime();
#endif
#ifdef MPI_CODE_VERSION
	for(i=mpi_info.start_site; i<mpi_info.end_site; i++)
		if (curr_exchange->likelihood_is_scaled() == FALSE)
			final_lnL_s+=(double)log(prob_w_rate(i,0));
		else {
			//prob=prob_w_rate(i,0);
			final_lnL_s+=(double)(log(prob_w_rate(i,0))-log(LIKE_SCALE));
			//final_lnL+=(double)(log(prob)-log(LIKE_SCALE));
			//cout<<i<<": "<<prob<<"\t"<<final_lnL<<endl;
		}
	MPI_Allreduce(&final_lnL_s, &final_lnL, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#else
	for (i=0; i<curr_exchange->get_num_localities(); i++) {
	   	if (curr_exchange->likelihood_is_scaled() == FALSE)
			final_lnL+=(double)log(prob_w_rate(i,0));
		else {
			//prob=prob_w_rate(i,0);
			final_lnL+=(double)(log(prob_w_rate(i,0))-log(LIKE_SCALE));
			//final_lnL+=(double)(log(prob)-log(LIKE_SCALE));
			//cout<<i<<": "<<prob<<"\t"<<final_lnL<<endl;
		}
	}
   
 
#endif
  last_lnL=final_lnL;   
#ifdef DO_TIME
	posttime=RunTime();
	cout<<"Runtime for likelihood calc: "<<posttime-pretime<<endl;
#endif
  return(final_lnL);

}


double Like_model::find_ln_like_fixed_site_rates() 
{
	int i, j; 
	double site_p, final_lnL=0;


	for (i=0; i<curr_exchange->get_num_localities(); i++) {
		site_p=0;
		for(j=0; j<curr_exchange->get_num_rates(); j++)
			site_p+=curr_exchange->get_site_rate_prob(i, j)*prob_w_rate(i, j);
		if (curr_exchange->likelihood_is_scaled() == FALSE)
			final_lnL+=log(site_p);
		else
			final_lnL+=(log(site_p) - log(LIKE_SCALE));
	}
 
  
	last_lnL=final_lnL;   
	return(final_lnL);

}



double Like_model::find_ln_like_w_rates()
  //This function finds the lnL using transprob matrix
  //and using the site rate assignments stored in the exchange object
  //(There are three possible rate parameterizations:  all sites evolve at 1 rate,
  //each codon position evolves a a different rate (nucleotide evolution models only) and
  //each site evolves at an arbitrary rate set in the exchange)
{
  int i, j, k; 
  double final_lnL=0;
 
  if (curr_exchange->modeling_codons()==TRUE)
    {
      if( (curr_exchange->branch_basefreqs()==FALSE))
		{
		  for(k=0; k<curr_exchange->get_num_rates(); k++)	{
				for(j=0; j<curr_exchange->get_num_nonsyn_patterns(); j++) {
					i=0;
					while((*curr_tree)[i]->get_nonsyn_pattern() != j) i++;	   
					calc_codon_matrices(k, (*curr_tree)[i]);
					for (i=0; i<curr_exchange->get_num_branches(); i++)
						if ((*curr_tree)[i]->get_nonsyn_pattern() == j)
							calc_transprobs((*curr_tree)[i], k);
				}
		  }
		}
      else
		{
			for(k=0; k<curr_exchange->get_num_rates(); k++)	{
		
				for (i=0; i<curr_exchange->get_num_branches(); i++)
					{
						calc_codon_matrices(0, (*curr_tree)[i]);
						calc_transprobs((*curr_tree)[i], k);
					}
			}
		}
	}          
  
  else
	  for(k=0; k<curr_exchange->get_num_rates(); k++)	
		for (i=0; i<curr_exchange->get_num_branches(); i++)
			calc_transprobs((*curr_tree)[i], k);
      
		for (i=0; i<curr_exchange->get_num_localities(); i++) {
			if (curr_exchange->likelihood_is_scaled() == FALSE)
				final_lnL+=log(prob_w_rate(i, curr_exchange->get_site_rate_num(i)));
			else
				final_lnL+=(log(prob_w_rate(i, curr_exchange->get_site_rate_num(i)))- log(LIKE_SCALE));
		}

  last_lnL=final_lnL;   
  return(final_lnL);
     

} //End Like_model::find_ln_like_w_rates




double Like_model::find_lnL_on_tree(Tree *ctree)
  //Calculates the likelihood of the given tree, after first recalculating transprobs
  ///assuming the whole sequence evolves at a single rate (#0)
{
  int i, j; 
  double final_lnL=0, site_p;
  
  curr_tree=ctree;

  if (curr_exchange->modeling_codons()==TRUE)
    {
      if( (curr_exchange->branch_basefreqs()==FALSE))
	{
	  for(j=0; j<curr_exchange->get_num_nonsyn_patterns(); j++) {
	    i=0;
	    while((*curr_tree)[i]->get_nonsyn_pattern() != j) i++;	
	    calc_codon_matrices(0, (*curr_tree)[i]);
	    for (i=0; i<curr_exchange->get_num_branches(); i++)
	      if ((*curr_tree)[i]->get_nonsyn_pattern() == j)
		calc_transprobs((*curr_tree)[i], 0);
	  }
	}
      
      else
	{
	  for (i=0; i<curr_exchange->get_num_branches(); i++)
	    {
	      calc_codon_matrices(0, (*curr_tree)[i]);
	      calc_transprobs((*curr_tree)[i], 0);
	    }
	}
    }          
  
  else
    for (i=0; i<curr_exchange->get_num_branches(); i++)
      calc_transprobs((*curr_tree)[i], 0);
      
    
  for (i=0; i<curr_exchange->get_num_localities(); i++)
    {
      site_p=prob_w_rate(i, 0);
      if (site_p<=0 || site_p>=1.0)
			cout<<"Site: "<<i<<" prob: "<<site_p<<endl;
	  if (curr_exchange->likelihood_is_scaled() == FALSE)
		  final_lnL+=log(site_p);
	  else
		  final_lnL+=(log(site_p)-log(LIKE_SCALE));
  }
  last_lnL=final_lnL;   
  return(final_lnL);
}



void Like_model::return_lnL (int start, int end, int rate_num, double *new_rate, double new_sitelnLs[])
  //Calculates the likelihood of the part of the sequence between start and end using rate # rate_num 
  //and the rate parameters given in the new_rate array.
{
  int i, j;
  double temp;

  
  if (change_rate(rate_num, new_rate)==TRUE)
    {
      if (curr_exchange->modeling_codons()==TRUE)
	{
	  if( (curr_exchange->branch_basefreqs()==FALSE))
	    {
	      for(j=0; j<curr_exchange->get_num_nonsyn_patterns(); j++) {
		i=0;
		while((*curr_tree)[i]->get_nonsyn_pattern() != j) i++;	
		calc_codon_matrices(rate_num, (*curr_tree)[i]);
		for (i=0; i<curr_exchange->get_num_branches(); i++)
		  if ((*curr_tree)[i]->get_nonsyn_pattern() == j)
		      calc_transprobs((*curr_tree)[i], rate_num);
	      }
	    }

	  else
	    {
	      for (i=0; i<curr_exchange->get_num_branches(); i++)
		{
		  calc_codon_matrices(rate_num, (*curr_tree)[i]);
		  calc_transprobs((*curr_tree)[i], rate_num);
		}
	    }
	}          
    
      else
	for (i=0; i<curr_exchange->get_num_branches(); i++)
	  calc_transprobs((*curr_tree)[i], rate_num);
      
    }
  if (curr_exchange->get_num_taxa()>2)
    {
      for (i=start; i<end; i++)
	new_sitelnLs[i]=-log(prob_w_rate(i, rate_num));
      
    }
  else
    for (i=start; i<end; i++)
      {
	temp=0;
	for(j=0; j<curr_exchange->get_condlike_size(); j++)
	  temp+=root_freq(j)*(*curr_tree)[0]->get_trpb(0, j, (*curr_data)[0][i])*
	    (*curr_tree)[0]->get_trpb(0, j, (*curr_data)[1][i]);

	new_sitelnLs[i]=-log(temp);
      }
 
}  //End Like_model::return_lnL




double Like_model::min_lnL (double param[], int offset)
  //This is the function that the Powell class calls to get
  //a function value (i.e. it is the function to be minimized
  //It receives an array containing values for all 
  //of the parameters in the current model.  It resets probablity matrices
  //and rates and finally returns a new lnL
  //The params array and  param_type array are originally from this class
  //but the array of parameters (but not that of types) is sent to the
  //Powell class.  This means that they are 1 off from each other (the type 
  //matrix starts at 0 and the parameters at 1

{

  int i,j, k, kaksnum, pns_num, matrix_num, rate_num;
  double ratio, sum, prev_fix=0.0, prev_loss=0.0, prev_parallel=0.0;
  long double lnL=0.0;
  BOOL one_diff=FALSE, calc_rate_1=FALSE, rerun;

  if(just_found_opt==TRUE)
    {
      //optimize_branches(); 
      just_found_opt=FALSE;
    }
  


  for (i=0; i<curr_exchange->get_num_params(); i++)
    {
      switch (param_types[i])
	{
	case TRS_TRV:
	  if(fabs(fabs(param[i+offset])-get_trs_trv())>=FLOAT_TOL)
	    {
	      curr_exchange->set_trs_trv(fabs(param[i+offset]));
	      one_diff=TRUE;
	      if (curr_exchange->is_mol_clock_3()==TRUE) {
		for(j=0; j<3; j++) {
		  set_ratios((*curr_tree)[j]);
		  if (j==0 || j==1) {
		     if (curr_exchange->get_clock_type() ==KS_CLOCK)
		       (*curr_tree)[j]->set_brnlen(save_k_dupl/(*curr_tree)[j]->get_syn_ratio());
		     else if (curr_exchange->get_clock_type() ==KA_CLOCK)
		       (*curr_tree)[j]->set_brnlen(save_k_dupl/(*curr_tree)[j]->get_nonsyn_ratio());
		  }else{
		    if (curr_exchange->get_clock_type() ==KS_CLOCK)
		      (*curr_tree)[j]->set_brnlen(save_k_common/(*curr_tree)[j]->get_syn_ratio());
		    else if (curr_exchange->get_clock_type() ==KA_CLOCK)
		      (*curr_tree)[j]->set_brnlen(save_k_common/(*curr_tree)[j]->get_nonsyn_ratio());
		  }
		}	      
	      }
	    }
	  break;

	case BRANCH:
		if (curr_exchange->zero_len_brns_fixed() == FALSE) {
            if (fabs(param[i+offset]) < MIN_BRLEN) {
                (*curr_tree)[i-brn_start]->set_brnlen(MIN_BRLEN);
                one_diff=TRUE;
            }
            else {
                if ((fabs((*curr_tree)[i-brn_start]->get_brnlen()-fabs(param[i+offset]))>FLOAT_TOL) || ((fabs(param[i+offset]) > MIN_BRLEN) && (fabs((*curr_tree)[i-brn_start]->get_brnlen()) <=MIN_BRLEN) ))
                {
                    if (curr_exchange->get_model() != DUPL_2_RATE_NOSUBF_NOSTATE)
                        (*curr_tree)[i-brn_start]->set_brnlen(fabs(param[i+offset]));
                    else
                        (*curr_tree)[i-brn_start]->set_brnlen(fabs(param[i+offset])/2.0);
                    // (*curr_tree)[i+brn_start]->get_sibling()->set_brnlen(fabs(param[i+offset]));
                    one_diff=TRUE;
                }
            }
		}
		else 
		{
            if (fabs(param[i+offset]) < MIN_BRLEN) {
                (*curr_tree)[brn_index[i-brn_start]]->set_brnlen(MIN_BRLEN);
                one_diff=TRUE;
            }
            else {
                if ((fabs((*curr_tree)[brn_index[i-brn_start]]->get_brnlen()-fabs(param[i+offset]))>FLOAT_TOL) || ((fabs(param[i+offset]) > MIN_BRLEN) && (fabs((*curr_tree)[i-brn_start]->get_brnlen()) <=MIN_BRLEN) ))
                {
                    (*curr_tree)[brn_index[i-brn_start]]->set_brnlen(fabs(param[i+offset]));
                    // (*curr_tree)[i+brn_start]->get_sibling()->set_brnlen(fabs(param[i+offset]));
                    one_diff=TRUE;
                }	}

		}
	  break;
	  
	case PUR_PYR_SPLIT:
	  if(fabs(fabs(param[i+offset])-pur_pyr_split)>=FLOAT_TOL)
	    {
	      pur_pyr_split=fabs(param[i+offset]);
	      one_diff=TRUE;
	    }
	  break;

	case A_G_SPLIT:
	    if (fabs(fabs(param[i+offset])-a_g_split)>=FLOAT_TOL)
	    {
	      a_g_split=fabs(param[i+offset]);
	      one_diff=TRUE;
	    }
	    break;

	case C_T_SPLIT:
	   if (fabs(fabs(param[i+offset])-c_t_split)>=FLOAT_TOL)
	    {
	      c_t_split=fabs(param[i+offset]);
	      one_diff=TRUE;
	    }
	   break;

	case THIRD_POS_RATE:
	  if (fabs(fabs(param[i+offset])-curr_exchange->get_rate(1))>=FLOAT_TOL)
	    {
	     curr_exchange->set_rate(1, fabs(param[i+offset]));
	      one_diff=TRUE;
	      calc_rate_1=TRUE;
	    }
	   break;

	case RATE_PROB:
		if (fabs(exp(-fabs(param[i+offset]))-lookup_p_tree->get_param_value(i-rate_prob_param_index)) >= FLOAT_TOL) {
			one_diff=TRUE;
			lookup_p_tree->set_param_proportion(i-rate_prob_param_index, exp(-fabs(param[i+offset])));

			for(j=0; j<curr_exchange->get_num_rates(); j++) 
				curr_exchange->set_generic_site_rate_prob(j, lookup_p_tree->get_prob_value(j));
		}
		break;
	case SUB_TYPE:
	  kaksnum = prop_pns_index[i-pns_start];
	  rate_num = prop_rate_index[i-pns_start];
	 
	  if (fabs(curr_exchange->get_p_non_syn(kaksnum, rate_num)-exp(param[i+offset]))>=FLOAT_TOL)
	    {
	      curr_exchange->set_p_non_syn(kaksnum, rate_num, exp(param[i+offset]));
	      one_diff=TRUE;
	    }

	  if (curr_exchange->is_mol_clock_3()==TRUE) {
	  
	    set_ratios((*curr_tree)[kaksnum]);
	    if ((kaksnum==0) || (kaksnum==1)) {
	      if (curr_exchange->get_clock_type() ==KS_CLOCK)
		(*curr_tree)[kaksnum]->set_brnlen(save_k_dupl/(*curr_tree)[kaksnum]->get_syn_ratio());
	      else if (curr_exchange->get_clock_type() == KA_CLOCK)
		(*curr_tree)[kaksnum]->set_brnlen(save_k_dupl/(*curr_tree)[kaksnum]->get_nonsyn_ratio());
	    }
	    else{
	      if (curr_exchange->get_clock_type() == KS_CLOCK) 
		(*curr_tree)[kaksnum]->set_brnlen(save_k_common/(*curr_tree)[kaksnum]->get_syn_ratio());
	      else if (curr_exchange->get_clock_type() == KA_CLOCK) 
		(*curr_tree)[kaksnum]->set_brnlen(save_k_common/(*curr_tree)[kaksnum]->get_nonsyn_ratio());
	    }
	  }

	
	  break;

	case MATRIX_COEFF:
		pns_num=prop_pns_index[i-matcoeff_start];
		matrix_num=matrix_index[i-matcoeff_start];
		
		if (fabs(curr_exchange->get_matrix_coeff(matrix_num, pns_num)-param[i+offset])>=FLOAT_TOL)
	    {
	      curr_exchange->set_matrix_coeff(matrix_num, pns_num, param[i+offset]);
	      one_diff=TRUE;
	    }
	  

	break;
	case CLASS_SUB:
		kaksnum= prop_pns_index[i-pns_start];
		rate_num=prop_rate_index[i-pns_start];
	 
	  if (fabs(curr_exchange->get_p_inter_group(kaksnum, rate_num)-exp(param[i+offset]))>=FLOAT_TOL)
	    {
	      curr_exchange->set_p_inter_group(kaksnum, rate_num, exp(param[i+offset]));
	      one_diff=TRUE;
	    }
	  break;

	case AA_PROP_FAC:
		if (curr_exchange->allow_pos_LCAP_weights()==TRUE) {
			if (fabs(curr_exchange->get_aa_prop_fac(prop_pns_index[i-pns_start], prop_rate_index[i-pns_start], prop_index[i-pns_start])-
				param[i+offset])>=FLOAT_TOL)
			{
				curr_exchange->set_aa_prop_fact(prop_pns_index[i-pns_start], prop_rate_index[i-pns_start], prop_index[i-pns_start], param[i+offset]);
				one_diff=TRUE;
			}
		}
		else {
				if (fabs(curr_exchange->get_aa_prop_fac(prop_pns_index[i-pns_start], prop_rate_index[i-pns_start], prop_index[i-pns_start])+
					fabs(param[i+offset]))>=FLOAT_TOL)
			{
				curr_exchange->set_aa_prop_fact(prop_pns_index[i-pns_start], prop_rate_index[i-pns_start], prop_index[i-pns_start], -fabs(param[i+offset]));
				one_diff=TRUE;
			}

		}
	  break;
	  
	case SCALING_FAC:
		if(curr_exchange->allow_pos_LCAP_weights() ==TRUE) {
			if (fabs(curr_exchange->get_aa_prop_fac(prop_pns_index[i-pns_start], prop_rate_index[i-pns_start], SCALING)-param[i+offset])>=FLOAT_TOL)
			{
				curr_exchange->set_aa_prop_fact(prop_pns_index[i-pns_start], prop_rate_index[i-pns_start], SCALING, param[i+offset]);
				one_diff=TRUE;
			}
		}
		else {
			if (fabs(curr_exchange->get_aa_prop_fac(prop_pns_index[i-pns_start], prop_rate_index[i-pns_start], SCALING)+fabs(param[i+offset]))>=FLOAT_TOL)
			{
				curr_exchange->set_aa_prop_fact(prop_pns_index[i-pns_start], prop_rate_index[i-pns_start], SCALING, -fabs(param[i+offset]));
				one_diff=TRUE;
			}

		}
	  break;

	case THREE_TREE_COMMON:  
	  if (fabs(save_k_common-fabs(param[i+offset]))>=FLOAT_TOL)
	    {
	      //cout<<"Changing common Ks: "<<save_k_common<<" to "<<fabs(param[i+offset])<<" KaKs: "<<(*curr_tree)[2]->get_p_nonsyn()<<endl;
	      save_k_common=fabs(param[i+offset]);
	      (*curr_tree)[2]->set_brnlen(save_k_common/(*curr_tree)[2]->get_syn_ratio());
	      
	      one_diff=TRUE;
	   }
	  break;

	case THREE_TREE_SPLIT: 
	 	  if (fabs(save_k_dupl-fabs(param[i+offset]))>=FLOAT_TOL)
	   {
	     save_k_dupl=fabs(param[i+offset]);
	     (*curr_tree)[0]->set_brnlen(save_k_dupl/(*curr_tree)[0]->get_syn_ratio());
	     (*curr_tree)[1]->set_brnlen(save_k_dupl/(*curr_tree)[1]->get_syn_ratio());
	    
	     one_diff=TRUE;
	   }
	  break;

	case THREE_TREE_COMMON_KA:  
	  if (fabs(save_k_common-fabs(param[i+offset]))>=FLOAT_TOL)
	    { 
	      //cout<<"Changing common Ka: "<<save_k_common<<" to "<<fabs(param[i+offset])<<" KaKs: "<<(*curr_tree)[2]->get_p_nonsyn()<<endl;
	      save_k_common=fabs(param[i+offset]);
	      (*curr_tree)[2]->set_brnlen(save_k_common/(*curr_tree)[2]->get_nonsyn_ratio());
	      
	      one_diff=TRUE;
	   }
	  break;

	case THREE_TREE_SPLIT_KA: 	 
	  if (fabs(save_k_dupl-fabs(param[i+offset]))>=FLOAT_TOL)
	   {
	     save_k_dupl=fabs(param[i+offset]); 
	     (*curr_tree)[0]->set_brnlen(save_k_dupl/(*curr_tree)[0]->get_nonsyn_ratio());
	     (*curr_tree)[1]->set_brnlen(save_k_dupl/(*curr_tree)[1]->get_nonsyn_ratio());

	     one_diff=TRUE;
	   }
	  break;
	case DUPL_PARALLEL_RATE:
		if (fabs(curr_exchange->get_dupl_parallel_rate() - fabs(param[i+offset]))>=FLOAT_TOL) {
			prev_parallel=curr_exchange->get_dupl_parallel_rate();
			if (fabs(param[i+offset]) < MAX_PARAM)
				curr_exchange->set_dupl_parallel_rate(fabs(param[i+offset]));
			else
				curr_exchange->set_dupl_parallel_rate(MAX_PARAM);
			one_diff=TRUE;
		}
		break;
	case DUPL_FIX_RATE:
		if (fabs(curr_exchange->get_dupl_fix_rate() - fabs(param[i+offset]))>=FLOAT_TOL) {
			prev_fix=curr_exchange->get_dupl_fix_rate();
			if (fabs(param[i+offset]) < MAX_PARAM)
				curr_exchange->set_dupl_fix_rate(fabs(param[i+offset]));
			else
				curr_exchange->set_dupl_fix_rate(MAX_PARAM);
			one_diff=TRUE;
		}
		break;
	case STRAND_SWITCH_PROB:
		if (fabs(curr_exchange->get_strand_switch_prob() - exp(-fabs(param[i+offset]))) >= FLOAT_TOL) {
			curr_exchange->set_strand_switch_prob(exp(-fabs(param[i+offset])));
			one_diff=TRUE;
		}
		break;
	case CON_REL_LOSS_RATE:
		if (fabs(curr_exchange->get_loss_rate_scale() - fabs(param[i+offset]))>=FLOAT_TOL) {
			prev_loss=curr_exchange->get_loss_rate_scale();
			if (fabs(param[i+offset]) < MAX_PARAM)
				curr_exchange->set_loss_rate_scale(fabs(param[i+offset]));
			else
				curr_exchange->set_loss_rate_scale(MAX_PARAM);
			one_diff=TRUE;
		}
		break;
	case CON_REL_FIX_RATE:
		if (fabs(curr_exchange->get_fix_rate_scale() - fabs(param[i+offset]))>=FLOAT_TOL) {
			if (fabs(param[i+offset]) < MAX_PARAM)
				curr_exchange->set_fix_rate_scale(fabs(param[i+offset]));
			else
				curr_exchange->set_fix_rate_scale(MAX_PARAM);
			one_diff=TRUE;
		}
		break;
	case CON_FIX_LOSS_RATE:
		if (fabs(curr_exchange->get_fix_loss_rate() - fabs(param[i+offset])) >= FLOAT_TOL) {
			curr_exchange->set_fix_loss_rate(fabs(param[i+offset]));
			one_diff=TRUE;
		}
		break;
	case SNP_GAIN_LOSS_RATIO:
		if (fabs(curr_exchange->get_snp_gains_to_losses() - exp(param[i+offset])) >= FLOAT_TOL) {
			curr_exchange->set_snp_gains_to_losses(exp(param[i+offset]));
			one_diff=TRUE;
		}
		break;
	case SNP_LOSS_SWITCH_RATIO:
		if (fabs(curr_exchange->get_snp_loss_to_switch() - exp(param[i+offset])) >= FLOAT_TOL) {
			curr_exchange->set_snp_loss_to_switch(exp(param[i+offset]));
			one_diff=TRUE;
		}
		break;
    case ARBITRARY_PARAM:
        set_arb_param(i, param[i+offset]);
        one_diff=TRUE;
        break;
	}
#if defined (SHOW_PROGRESS)
 
     cout<<"("<<setprecision(8)<<param[i+offset]<<", "<<param_types[i]<<") ";
#endif
    } //End for loop through parameters


  if(one_diff==TRUE)
    {
        
#if defined (SHOW_PROGRESS)
        cout<<"Parameter change--recalculating transprobs\n";
#endif
        
      curr_exchange->set_basefreqs(pur_pyr_split/((pur_pyr_split+1)*
						  (1+(1.0/a_g_split))), 0);
      curr_exchange->set_basefreqs(pur_pyr_split/
				   ((pur_pyr_split+1)*(1+a_g_split)),2);
      curr_exchange->set_basefreqs(1.0/((pur_pyr_split+1)*(1+(1.0/c_t_split))),1);
      curr_exchange->set_basefreqs(1.0/((pur_pyr_split+1)*(1+c_t_split)),3);

#if defined (SHOW_PROGRESS_BASEFREQ)
      cout<<"   Pur/Pyr Split: "<<pur_pyr_split<<" AG Split: "<<a_g_split<<" CT split: "
	  <<c_t_split<<" Freqs: "<<get_basefreq(0,0)<<" "<<get_basefreq(1,0)<<" "
	  <<get_basefreq(2,0)<<" "<<get_basefreq(3,0)<<" "<<" kappa: "<<get_trs_trv()<<" ";
#endif
      
      if (curr_exchange->modeling_codons()==TRUE)
		{
			if( (curr_exchange->branch_basefreqs()==FALSE))
			{
				for(k=0; k<curr_exchange->get_num_rates(); k++) {
					for(j=0; j<curr_exchange->get_num_nonsyn_patterns(); j++) {
						i=0;
						while((*curr_tree)[i]->get_nonsyn_pattern() != j) i++;

						calc_codon_matrices(k, (*curr_tree)[i]);
					
						for (i=0; i<curr_exchange->get_num_branches(); i++)
							if ((*curr_tree)[i]->get_nonsyn_pattern()==j)
								calc_transprobs((*curr_tree)[i], k);
					}
				}
			}

			else
			{
				for (i=0; i<curr_exchange->get_num_branches(); i++)
				{
					calc_codon_matrices(0, (*curr_tree)[i]);
					calc_transprobs((*curr_tree)[i], 0);
				}
			}
	     
        }
    
      else if (curr_exchange->get_num_taxa()>2)
		{
			for (i=0; i<curr_exchange->get_num_branches(); i++)
				calc_transprobs((*curr_tree)[i], 0);
			if (calc_rate_1==TRUE)
				for (i=0; i<curr_exchange->get_num_branches(); i++)
					calc_transprobs((*curr_tree)[i], 1);
        }
      else if ((curr_exchange->get_num_taxa() ==2) && (curr_exchange->is_rooted_tree() == TRUE)) {
          for (i=0; i<curr_exchange->get_num_branches(); i++)
              calc_transprobs((*curr_tree)[i], 0);
      }

     lnL=find_appropriate_ln_like();
	  //Kludge for when Powell's routinue takes too big a step
	  if(!((lnL<0) && (lnL >-1e32))) {
#ifdef PRINT_WARNINGS
			cerr<<"WARNING: Likelihood function returned invalid number: "<<lnL<<endl;
#endif
		  lnL=-1e32;
			
		  rerun=FALSE;
		  if (prev_fix != curr_exchange->get_dupl_fix_rate()) {
			curr_exchange->set_dupl_fix_rate(prev_fix);
			rerun=TRUE;
		  }
		  if (prev_parallel != curr_exchange->get_dupl_parallel_rate()) {
			  rerun=TRUE;
			  curr_exchange->set_dupl_parallel_rate(prev_parallel);
		  }
		  if (prev_loss != curr_exchange->get_loss_rate_scale()) {
			  rerun=TRUE;
			  curr_exchange->set_loss_rate_scale(prev_loss);
		  }
		  if (rerun== TRUE) {
			lnL=find_appropriate_ln_like();
			if(!((lnL<0) && (lnL >-1e32))) 
			 lnL=-1e32;
		  }
	  }
      last_lnL=lnL;
    }

  else
    lnL=last_lnL;

#if defined (SHOW_PROGRESS)

  cout<<"\tlnL: "<<setprecision(5)<<setw(7)<<lnL<<" (offset: "<<offset<<")"<<endl;
#endif

  return(-lnL);

}  //End Like_model::min_lnL




void Like_model::min_single_param (PARAM_TYPE parameter, double tol)
  //Called by Powell for models with only a single global parameter 
  //(These models are the K2P model and the MG94 codon model with
  //a Juke-Cantor model of nucleotide substitution)
{
  int i;
  double low_t, guess_t, actual_t, low_lnL, guess_lnL, actual_lnL, final_min, last_lnL;
  double (*thefunc)(double);

  thefunc=&single_min;

  newton_branch_opt();
  newton_branch_opt();
  
  single_param_id=parameter;

  switch(single_param_id)
    {
    case TRS_TRV:
      low_t=0;
      guess_t=10;
      break;
      
    case SUB_TYPE:
      low_t=1e-10;
      guess_t=0;
      break;
    
    default:
        low_t=1e-10;
        guess_t=0;
        break;
    }

 
  mnbrak(&low_t, &guess_t, &actual_t, &low_lnL, &guess_lnL, &actual_lnL, thefunc);
  last_lnL =actual_lnL+10;
  while (last_lnL-actual_lnL >FLOAT_TOL) {
    last_lnL=actual_lnL;
    actual_lnL=brent(low_t, actual_t, guess_t, thefunc, tol, &final_min); 
    actual_lnL=newton_branch_opt(0.001);
  }
 
  cout<<"Best trs/trv: "<<fabs(final_min)<<" final like: "<<setw(12)<<setprecision(10)<<actual_lnL<<endl;
 
  for (i=0; i<5; i++)
    newton_branch_opt();
  describe_results();
}  //End Like_model::min_trs_trv




double Like_model::optimize_single_param (double param_val)
  //The objective function called by Powell's routine when using the above
  //function to optimize a singel paramete
{
  int i, j;  

  
  switch(single_param_id)
    {
    case TRS_TRV:
      if (fabs(get_trs_trv()-fabs(param_val))>FLOAT_TOL) 
	{
	  curr_exchange->set_trs_trv(fabs(param_val));
      
	  for (i=0; i<curr_exchange->get_num_branches(); i++)
	    calc_transprobs((*curr_tree)[i], 0);
	}
      break;
      
    case SUB_TYPE:
      
      if (fabs(curr_exchange->get_p_non_syn(0)-exp(-fabs(param_val)))>=FLOAT_TOL)
	    {
	      curr_exchange->set_p_non_syn(0, exp(-fabs(param_val)));
	          
	      if( (curr_exchange->branch_basefreqs()==FALSE))
			{ 
				for(j=0; j<curr_exchange->get_num_nonsyn_patterns(); j++) {
					i=0;
					while((*curr_tree)[i]->get_nonsyn_pattern() != j) i++;	
					calc_codon_matrices(0, (*curr_tree)[i]);
					for (i=0; i<curr_exchange->get_num_branches(); i++)
						if ((*curr_tree)[i]->get_nonsyn_pattern() == j)
							calc_transprobs((*curr_tree)[i], 0);
				}
			}
	      
	      else
		{
		  for (i=0; i<curr_exchange->get_num_branches(); i++)
		    {
		      calc_codon_matrices(0, (*curr_tree)[i]);
		      calc_transprobs((*curr_tree)[i], 0);
		    }
		}
	      
	      
	     
	      
	    }
      break;
        
    default:
        break;
        
    }
   
 
    
  //cout<<"New trs/trv: "<<fabs(param_val)<<" lnL: "<<find_ln_like_single_rate()<<endl;
  return(-find_ln_like_single_rate());
}  //End Like_model::optimize_trs_trv


double Like_model::newton_branch_opt()
{
  return(newton_branch_opt(FLOAT_TOL));
}

double Like_model::newton_branch_opt(double tol)
{
  int i, temp;
  double retval;

  for(i=0; i<curr_exchange->get_num_branches(); i++)
    {
        if ((*curr_tree)[brn_order[i]]->is_fixed_zero()==FALSE ) {
            if (((*curr_tree)[brn_order[i]]!=curr_tree->find_null_branch_id()) &&
                ((*curr_tree)[brn_order[i]]!=curr_tree->find_root()))
                    //cout<<"pre lnL: "<<find_appropriate_ln_like()<<": "<<(*curr_tree)[brn_order[i]]->get_name()<<endl;
                retval=newton_opt_a_branch((*curr_tree)[brn_order[i]], tol);
                //cout<<"Post lnL: "<<find_appropriate_ln_like()<<endl;
        }
    
	
    }
  
  if (start_swap == 0) {
    i=0;
    while (i<curr_exchange->get_num_branches()-1) {
      temp=brn_order[i];
      brn_order[i]=brn_order[i+1];
      brn_order[i+1]=temp;
      i+=2;
    }
    start_swap=1;
  }
  else {
    i=curr_exchange->get_num_branches()-1;
    while(i>0) {
      temp=brn_order[i];
      brn_order[i]=brn_order[i-1];
      brn_order[i-1]=temp;
      i-=2;
    }
    start_swap=0;
  }
  
  return(retval);
} //End Like_model::newton_branch_opt



void Like_model::setup_branch_opt(Branch *branch)
{
  int i,j,k, l;
  Branch *left_root;


  if (curr_exchange->get_num_branches()>3) {
    curr_tree->re_root_tree(branch);
   
   
    //brn_op_cbrn is what we would call the right child of the root
    brn_op_cbrn=branch;
    left_root=curr_tree->find_null_branch_id();

	
    left_root->set_brnlen(0.0);
    
	

    if (brn_op_cbrn->is_tip()==FALSE)  {
		for(j=0; j<curr_exchange->get_num_localities(); j++)
		{
		  for (l=0; l<curr_exchange->get_num_rates(); l++) {
			partial_prob_w_rate(j, curr_tree->get_leftmost_tip(), l, curr_tree->find_root());
		
			for (k=0; k<curr_exchange->get_condlike_size(); k++)
			{
				brncondprobs1[j][l][k]=left_root->get_cond_prob(k)*root_freq(k);
				brncondprobs2[j][l][k]=brn_op_cbrn->get_cond_prob(k);
			} 
		  }
		}
		
	}
    
    else
      for(j=0; j<curr_exchange->get_num_localities(); j++) {
		  for(l=0; l<curr_exchange->get_num_rates(); l++) {
				partial_prob_w_rate( j, curr_tree->get_leftmost_tip(), 
					   l, curr_tree->find_root());
	  
				for (k=0; k<curr_exchange->get_condlike_size(); k++)
					brncondprobs1[j][l][k]=left_root->get_cond_prob(k)*root_freq(k);
			}
	  }
  }
  else 
    brn_op_cbrn=branch;
  
}


double Like_model::newton_opt_a_branch(Branch *branch, double tol) 
  //This is the new branch lenght optimizer that uses the first and second deriviatives of
  //the likelihood to optimize branch lengths using Newton's method.  Details are
  //described in Yang 2000, J. Mol. Evol. 51:423-432
{
  int i,zero_cnt=0;
  double  next_brn, last_brn, new_brn, alpha=1.0;
  double like, like_prime, like_double_prime, last_like, last_like_p, last_like_double_p;
  BOOL newt_continue;

 
  setup_branch_opt(branch);

  //Newton optimization of branch lengths
  last_brn=next_brn=brn_op_cbrn->get_brnlen();
  newt_continue=TRUE;
  branch_opt_like_and_deriv(next_brn, like, like_prime, like_double_prime, TRUE); 
#ifdef SHOW_PROGRESS
  cout<<setw(12)<<setprecision(8)<<"Initial brn: "<<brn_op_cbrn->get_brnlen()<<" Like: "<<like<<" Likep: "<<like_prime
      <<" Like double prime: "<<like_double_prime<<" tol: "<<tol<<endl<<flush;
#endif
  if (((like_prime < 1e200) && (like_prime > -1e200)) && (like_double_prime != 0) &&
      ((like_double_prime < 1e200) && (like_double_prime > -1e200)))
    {
      new_brn=next_brn+alpha*(-like_prime/like_double_prime);
      if ((new_brn < MIN_BRLEN) && (new_brn > 0))
	 new_brn=MIN_BRLEN;
	
      if (new_brn < 0) { 
	while(new_brn < 0) {
	  alpha *= 0.5;
	  new_brn=next_brn+alpha*(-like_prime/like_double_prime); 
	}		
	if (new_brn < MIN_BRLEN) {
	    new_brn=MIN_BRLEN;
	    zero_cnt++;
	}
	  
      }
      next_brn=new_brn;
      last_like=like;
      last_like_p=like_prime;
      last_like_double_p=like_double_prime;
     
      do {
	
			branch_opt_like_and_deriv(next_brn, like, like_prime, like_double_prime, FALSE);
			if (((like >-1e200) && (like <=0)) &&  ((like_prime < 1e200) && 
				(like_prime > -1e200)) && (like_double_prime != 0) &&
				((like_double_prime < 1e200) && (like_double_prime > -1e200)))
			{
				new_brn=next_brn+alpha*(-like_prime/like_double_prime);
#ifdef SHOW_PROGRESS
				cout<<"New like:" <<like<<" like p: "<<like_prime<<" like dp: "<<like_double_prime<<" Next: "<<next_brn
					<<" New brn: "<<new_brn<<endl;
#endif
				if ((new_brn < MIN_BRLEN) && (new_brn > 0))
				{ 
					zero_cnt++;
					new_brn=MIN_BRLEN;  
			
					if (zero_cnt == 3)
						newt_continue=FALSE;
				}
	    
				if (new_brn < 0) { 
					while(new_brn < 0) {
						alpha *= 0.5;
						new_brn=next_brn+alpha*(-like_prime/like_double_prime); 
					}		
					if (new_brn < MIN_BRLEN)
					{
						zero_cnt++;
						new_brn=MIN_BRLEN;
						if (zero_cnt == 3)
							newt_continue=FALSE;
					}
				}
				else if ((last_like>like) && (fabs(last_like-like)>tol)) {
	    
					alpha *= 0.5;
					new_brn=last_brn+alpha*(-last_like_p/last_like_double_p);   
	    
				}
				else if (fabs(last_like-like)<tol) {
 
					newt_continue=FALSE;
					if (last_like > like)
						new_brn=last_brn;
					else
						new_brn=next_brn;
				}		
				else {
					last_brn=next_brn;
					next_brn=new_brn;
					last_like=like;
					last_like_p=like_prime;
					last_like_double_p=like_double_prime;
				}
				next_brn=new_brn;
			}
			else {
				if ((last_like-like)>tol)
				{  
					next_brn=last_brn;
					like=last_like;
				}
				newt_continue=FALSE;
			}
      } while (newt_continue);
	  }
  brn_op_cbrn->set_brnlen(next_brn);

  for(i=0; i<curr_exchange->get_num_rates(); i++) {
	if (curr_exchange->modeling_codons() ==TRUE)
		calc_codon_matrices(i, brn_op_cbrn);    
	calc_transprobs(brn_op_cbrn, i);
  }
#ifdef SHOW_PROGRESS	
    cout<<"Brn: "<<brn_op_cbrn->get_brn_num()<<" New: "<<brn_op_cbrn->get_brnlen()<<" Like: "<<like<<endl;
#endif
	return(like);
}



double Like_model::newton_single_iter(Branch *branch, double step_size, double tol)
{
  int i;
  double  next_brn, last_brn, new_brn;	cout<<"Is not tip\n"<<flush;
  double last_like, like, last_like_p, like_prime, last_like_dp, like_double_prime;
 

  setup_branch_opt(branch);
 
  
  //Newton optimization of branch lengths
  last_brn=next_brn=brn_op_cbrn->get_brnlen();
 
  branch_opt_like_and_deriv(next_brn, like, like_prime, like_double_prime, TRUE); 
 
 last_like=like;
 //cout<<"Initial brn: "<<brn_op_cbrn->get_brnlen()<<" Like: "<<like<<" Likep: "<<like_prime
 // <<" Like dp: "<<like_double_prime<<endl;
  if (((like_prime < 1e200) && (like_prime > -1e200)) && (like_double_prime != 0) &&
	    ((like_double_prime < 1e200) && (like_double_prime > -1e200)))
    {
      new_brn=next_brn+step_size*(-like_prime/like_double_prime);
      //cout<<"New "<<new_brn;
      if ((new_brn < MIN_BRLEN) && (new_brn > 0))
	 new_brn=MIN_BRLEN;
	
      if (new_brn < 0) { 
	while(new_brn < 0) {
	  step_size *= 0.5;
	  new_brn=next_brn+step_size*(-like_prime/like_double_prime); 
	}		
	if (new_brn < MIN_BRLEN) 
	    new_brn=MIN_BRLEN;
      }
	  
      next_brn=new_brn; 
      branch_opt_like_and_deriv(next_brn, like, like_prime, like_double_prime, TRUE);
      //cout<< "Calc like: "<<like<<endl;
    }
  
  if ((last_like-like)>tol) {
    next_brn=last_brn;
    like=last_like;
  }
  branch->set_brnlen(next_brn);
  calc_transprobs(branch, 0);
  //cout<<" Brn: "<<next_brn<<"Final like: "<<like<<endl;
  return(like);
}




void Like_model::optimize_branches()
  //This is an old branch length optimizer that has be superseeded by the above function
{
  int i, j, k;
  double low_brn, guess_brn, actual_brn, low_lnL, guess_lnL, actual_lnL, final_min;
  double (*thefunc)(double);
  Branch *left_root;


  thefunc=&brn_min;
 
  for(i=0; i<curr_exchange->get_num_branches()-1; i++)
    {
      setup_branch_opt((*curr_tree)[i]);
	  
      low_brn=0;
      guess_brn=1.25*brn_op_cbrn->get_brnlen()+2*BRN_TOL;
      
      actual_brn=brn_op_cbrn->get_brnlen();
      
      actual_lnL=brent(low_brn, actual_brn, guess_brn, thefunc, 20*BRN_TOL, &final_min);
      actual_lnL=optimize_a_branch(final_min);
    } 

}  //End Like_model::optimize_branches




double Like_model::optimize_a_branch (double new_ut)
  //This is the objective function for the optimize_branches() function--it is also
  //obsolete
{
  int i,j,k;
  double retval=0, site_p=0;
  if (fabs(brn_op_cbrn->get_brnlen()-fabs(new_ut))>FLOAT_TOL) 
    {
      brn_op_cbrn->set_brnlen(fabs(new_ut));
      calc_transprobs(brn_op_cbrn, 0);
    
      if (brn_op_cbrn->is_tip()==FALSE) 
	{
	  for(i=0; i<curr_exchange->get_num_localities(); i++)
	    {
	      site_p=0;
	      for(j=0; j<curr_exchange->get_condlike_size(); j++)
		for(k=0; k<curr_exchange->get_condlike_size(); k++)
		  site_p+=brn_op_cbrn->get_trpb(0, j, k)*brncondprobs2[i][0][k]*brncondprobs1[i][0][j];
	      retval+=log(site_p);
	    }
	  last_lnL=retval;
	}
      
      else
	for (i=0; i<curr_exchange->get_num_localities(); i++)
	  {
	    site_p=0;
	    for(j=0; j<curr_exchange->get_condlike_size(); j++)
	      site_p+=brn_op_cbrn->get_trpb(0, j, (*curr_data)[brn_op_cbrn->get_taxa_id()][i])*
		brncondprobs1[i][0][j];
	    retval+=log(site_p);
	  }
      last_lnL=retval;
    }
  else
    retval=last_lnL;
 
  return(retval);
}  //End Like_model::optimize_a_branch




void Like_model::branch_opt_like_and_deriv(double ut, double &likelihood, double &likelihood_prime, double &likelihood_double_prime, BOOL first)
  //This is the objective function for the newton_branch_opt() branch length optimizer
{ //Note that I have not implemented likelihood scaling for this function -- GCC 10/1/08
  int i, j, k, l;
  double site_rate_p, site_rate_p_prime, site_rate_p_double_prime,
	  site_p, site_p_prime, site_p_double_prime;

  likelihood=0;
  likelihood_prime=0;
  likelihood_double_prime=0;
  
	
  
  brn_op_cbrn->set_brnlen(fabs(ut));
  

  if (curr_exchange->modeling_codons()==TRUE) {
	if ((curr_exchange->get_num_rates() >1) ) {
	  for(i=0; i<curr_exchange->get_num_rates(); i++) {
		  //Calc transprobs for each rate
			calc_codon_matrices(i, brn_op_cbrn);   

			calc_transprobs(brn_op_cbrn, i);
			calc_first_derivatives(brn_op_cbrn, fabs(ut), i);
			calc_second_derivatives(brn_op_cbrn, fabs(ut), i);
	  }
		
	}
	else {
		if (last_codonmat != brn_op_cbrn->get_nonsyn_pattern())
			calc_codon_matrices(0, brn_op_cbrn);   

		calc_transprobs(brn_op_cbrn, 0);
		calc_first_derivatives(brn_op_cbrn, fabs(ut));
		calc_second_derivatives(brn_op_cbrn, fabs(ut));
	}
  }
  else {
	    //Nucleotide model
		for(i=0; i<curr_exchange->get_num_rates(); i++) {
			calc_transprobs(brn_op_cbrn, i);
			calc_first_derivatives(brn_op_cbrn, fabs(ut), i);
			calc_second_derivatives(brn_op_cbrn, fabs(ut), i);
		}

  }
   if (curr_exchange->get_num_branches() > 3) {
    	   
   if (brn_op_cbrn->is_tip()==FALSE) 
      {
		for(i=0; i<curr_exchange->get_num_localities(); i++)
		{
			site_p=site_p_prime=site_p_double_prime=0;

			for(l=0; l<curr_exchange->get_num_rates(); l++) {
				site_rate_p=site_rate_p_prime=site_rate_p_double_prime=0;
				for(j=0; j<curr_exchange->get_condlike_size(); j++){
					for(k=0; k<curr_exchange->get_condlike_size(); k++) {
						site_rate_p+=brn_op_cbrn->get_trpb(l, j, k)*brncondprobs2[i][l][k]*
							brncondprobs1[i][l][j];
						
						site_rate_p_prime+=brn_op_cbrn->get_trpb_prime(j, k, l)*
							brncondprobs2[i][l][k]*brncondprobs1[i][l][j];		
						site_rate_p_double_prime+=brn_op_cbrn->get_trpb_double_prime(j, k, l)*
							brncondprobs2[i][l][k]*brncondprobs1[i][l][j];
					}
				}

				site_p+=curr_exchange->get_site_rate_prob(i, l)*site_rate_p;
				site_p_prime+=curr_exchange->get_site_rate_prob(i, l)*site_rate_p_prime;
				site_p_double_prime+=curr_exchange->get_site_rate_prob(i, l)*site_rate_p_double_prime;
			}
			
		
			likelihood+=log(site_p);
		    likelihood_prime+=(site_p_prime/site_p);
			likelihood_double_prime+=(site_p*site_p_double_prime-site_p_prime*site_p_prime)/(site_p*site_p);
		}
	
     }
      
   else {
	for (i=0; i<curr_exchange->get_num_localities(); i++)
	{
	  site_p=site_p_prime=site_p_double_prime=0.0;
	  for (l=0; l<curr_exchange->get_num_rates(); l++) {
		  site_rate_p=site_rate_p_prime=site_rate_p_double_prime=0;
		for(j=0; j<curr_exchange->get_condlike_size(); j++) {
			site_rate_p+=brn_op_cbrn->get_trpb(l, j, (*curr_data)[brn_op_cbrn->get_taxa_id()][i])
				*brncondprobs1[i][l][j];
	      
			if ( (*curr_data)[brn_op_cbrn->get_taxa_id()][i]<curr_exchange->get_condlike_size()) {
	     	    
				site_rate_p_prime+=brn_op_cbrn->get_trpb_prime(j, 
					(*curr_data)[brn_op_cbrn->get_taxa_id()][i], l)*brncondprobs1[i][l][j];
		
				site_rate_p_double_prime+=brn_op_cbrn->get_trpb_double_prime(j, 
					(*curr_data)[brn_op_cbrn->get_taxa_id()][i], l)*brncondprobs1[i][l][j];
			}
		}
		site_p+=curr_exchange->get_site_rate_prob(i, l)*site_rate_p;
		site_p_prime+=curr_exchange->get_site_rate_prob(i, l)*site_rate_p_prime;
		site_p_double_prime+=curr_exchange->get_site_rate_prob(i, l)*site_rate_p_double_prime;
	  }
	  
	  likelihood+=log(site_p);
	  if (site_p > 0) {
	    likelihood_prime+=(site_p_prime/site_p);
	    //cout<<"P: "<<site_p<<" site_p_prime: "<<site_p_prime<<" like p: "<<likelihood_prime<<endl;
	    likelihood_double_prime+=(site_p*site_p_double_prime-site_p_prime*site_p_prime)/(site_p*site_p);
	  }
	}
   }
  }
  else {
   
    for (i=0; i<curr_exchange->get_num_localities(); i++)
      {
		site_p=site_p_prime=site_p_double_prime=0;
		
		if ((*curr_data)[brn_op_cbrn->get_sibling()->get_taxa_id()][i] < curr_exchange->get_condlike_size()) 
			site_p=brn_op_cbrn->get_trpb(curr_exchange->get_site_rate_num(i), (*curr_data)[brn_op_cbrn->get_sibling()->get_taxa_id()][i],
				(*curr_data)[brn_op_cbrn->get_taxa_id()][i])*
				root_freq((*curr_data)[brn_op_cbrn->get_sibling()->get_taxa_id()][i]);
		
		else 
			//If there is a gap in the "root" site, just switch our perspective to the other taxa
			site_p=brn_op_cbrn->get_trpb(curr_exchange->get_site_rate_num(i), (*curr_data)[brn_op_cbrn->get_taxa_id()][i],
				(*curr_data)[brn_op_cbrn->get_sibling()->get_taxa_id()][i])*
				root_freq((*curr_data)[brn_op_cbrn->get_taxa_id()][i]);
		
		if (((*curr_data)[brn_op_cbrn->get_sibling()->get_taxa_id()][i]<curr_exchange->get_condlike_size()) &&
			((*curr_data)[brn_op_cbrn->get_taxa_id()][i]<curr_exchange->get_condlike_size())) {
			
			site_p_prime=brn_op_cbrn->get_trpb_prime((*curr_data)[brn_op_cbrn->get_sibling()->get_taxa_id()][i], 
				(*curr_data)[brn_op_cbrn->get_taxa_id()][i], curr_exchange->get_site_rate_num(i))*
				root_freq((*curr_data)[brn_op_cbrn->get_sibling()->get_taxa_id()][i]);
	  
			site_p_double_prime=brn_op_cbrn->get_trpb_double_prime((*curr_data)[brn_op_cbrn->get_sibling()->get_taxa_id()][i], 
				(*curr_data)[brn_op_cbrn->get_taxa_id()][i], curr_exchange->get_site_rate_num(i))*
				root_freq((*curr_data)[brn_op_cbrn->get_sibling()->get_taxa_id()][i]);
	}	
	//cout<<"P: "<<site_p<<" P`: "<<site_p_prime<<" P`` "<<site_p_double_prime<<endl;
	likelihood+=log(site_p);	
	likelihood_prime+=(site_p_prime/site_p);
	likelihood_double_prime+=(site_p*site_p_double_prime-(site_p_prime*site_p_prime))/(site_p*site_p);
      }
    
  }  

  

}



void Like_model::set_rate(double new_rate)
{
  curr_exchange->set_rate(1, new_rate);
}  //End Like_model::set_rate



void Like_model::set_expect_subs()
  //Using the current model of evolution, calculates the expected number of
  //subsitutions per site
{
  int i;
  
  for(i=0; i<curr_exchange->get_num_branches(); i++)
    {
      cout<<"Branch: "<<i<<" Len: "<<(*curr_tree)[i]->get_brnlen()<<" Expect: "<<get_expect_sub_site((*curr_tree)[i]->get_brnlen())<<endl;
    (*curr_tree)[i]->set_expect_subs_site(get_expect_sub_site((*curr_tree)[i]->get_brnlen()));
    }
    
    cout<<"BranchNames and sites\n";
    for(i=0; i<curr_exchange->get_num_branches(); i++)
    {
        (*curr_tree)[i]->name_this_branch();
        cout<<(*curr_tree)[i]->get_name()<<"\t"<<(*curr_tree)[i]->expect_subs_site()<<endl;
    }
}


void Like_model::get_rate_posts(int site, double *posts)
{
	int i;
	double total=0;
		
	for(i=0; i<curr_exchange->get_num_rates(); i++) {
		posts[i]=prob_w_rate(site, i);
		total+=posts[i];
	}

	for(i=0; i<curr_exchange->get_num_rates(); i++)
		posts[i]/=total;
	
}


double Like_model::get_obs_trs_trv_ratio()
  //Calcutates the observed trs/trv ratio for the kappa and base freqs
  //currently in the Exchange object
{
  cerr<<"Wrong get_obs_trs_trv_ration\n";
  return(0);
}


double Like_model::get_basefreq(int base, Branch *tree_branch)         //Virtual function
{
 cerr<<"Wrong get_basefreq\n";
 return(0);
}

double Like_model::get_basefreq(int base, int codon_pos, Branch *tree_branch)         //Virtual function
{
 cerr<<"Wrong get_basefreq\n";
 return(0);
}


double Like_model::root_freq(int site)                                //Virtual function
{
  cerr<<"Wrong root_freq\n";
  return(0);
}


double Like_model::get_expect_sub_site(double brlen)                  //Virtual function   
{
  cerr<<"get_expect_sub_site: Wrong one\n";
  return(0);
}

void Like_model::calc_first_derivatives(Branch *taxa, double ut)     //Virtual function
{
  cerr<<"Wrong calc_first_derivatives\n";
}

void Like_model::calc_second_derivatives(Branch *taxa, double ut)    //Virtual function 
{
  cerr<<"Wrong calc_second_derivatives\n";
}

void Like_model::calc_first_derivatives(Branch *taxa, double ut, int rate)     //Virtual function
{
  cerr<<"Wrong calc_first_derivatives\n";
}

void Like_model::calc_second_derivatives(Branch *taxa, double ut, int rate)    //Virtual function 
{
  cerr<<"Wrong calc_second_derivatives\n";
}
void Like_model::print_model()                                       //Virtual function
{
  //A virtual function implimented in derived classes
  cerr<<"print_model: Wrong function\n";
}


void Like_model::describe_results()                                 //Virtual function
{
  //This is a virtual function that is implemented in the specfic model classes
  cerr<<"describe_results: Wrong one!\n"; 
}

void Like_model::describe_nuc_results()                            //Virtual function
{
  //This is a virtual function that is implemented in the specfic model classes
  cerr<<"describe_nuc_results: Wrong one!\n"; 
}


int Like_model::rate_param_size()                                  //Virtual function
{
  cerr<<"Wrong rate_param_size\n";
  return(0);
}




Like_model::~Like_model()
  //Destructor for the Like_model class
{
  int i, j;
	

  delete[] site_lnLs;
  delete[] params;
  delete[] param_types;
  
  for (i=0; i<curr_exchange->get_num_localities(); i++)
    {
	  for (j=0; j<curr_exchange->get_num_rates(); j++) {
		delete[] brncondprobs1[i][j];
		delete[] brncondprobs2[i][j];
	  }
	  delete[] brncondprobs1[i];
	  delete[] brncondprobs2[i];
	
    }

  delete[] brncondprobs1;
  delete[] brncondprobs2;

  if (save_conprobs != 0) {
	  for (i=0; i<save_conprobs_size; i++) {
		for(j=0; j<curr_exchange->get_num_localities(); j++)
			delete[] save_conprobs[i][j];
		delete[] save_conprobs[i];
	  }
	  delete[] save_conprobs;
  }

  if (prop_index != 0)
	  delete[] prop_index;

  if(prop_pns_index != 0)
	  delete[] prop_pns_index;

  if(prop_rate_index != 0)
	  delete[] prop_rate_index;
	

  if(matrix_index != 0)
	  delete[] matrix_index;

  if(brn_index != 0)
	  delete[] brn_index;

  if(num_aa_props_per_rate != 0)
	delete[] num_aa_props_per_rate;

}  //End ~Like_model
//End of Like_model Public Functions





// Like_model Private Functions:
void Like_model::calc_codon_matrices (int rate_num, Branch *taxa)
  //This is virtual function implement only in codon models--
  //it allows them to use the lapack routines to calculate
  //the instantaneous transition probablity matrix for a set
  //of model parameters
{
  cerr<<"Wrong calc_codon_matrices\n";
}  //End Like_model::calc_codon_matrices



//Protected virtual functions implimented in derived classes
void Like_model::calc_transprobs(Branch *taxa, int rate_num)
{
  cerr<<"Wrong calc_transprobs\n";
}


void Like_model::num_params_model()                                        //Virtual function
{
  cerr<<"Wrong num_params_model\n";
}


void Like_model::intialize_parameters (double par[], PARAM_TYPE types[])   //Virtual function
{
  cerr<<"Wrong intialize_parameters\n";
}

double Like_model::q_matrix_entry (int start_base, int end_base, int codon_pos, Branch *tree_branch)        //Virtual function
  //Another function used only in codon models--this function calculates the instanteous
  //probablity of a given nucleotide substitution
{
  cerr<<"Wrong q_matrix_entry\n";
  return(0);
}


double Like_model::find_ut(Branch *taxa)                                    //Virtual function
{
  cerr<<"Wrong find_ut\n";
  return(0);
}


double Like_model::get_trs_trv()                                            //Virtual function
{
 cerr<<"Wrong get_trs_trv\n";
 return(0);
}


BOOL Like_model::change_rate(int rate_num, double *rate_info)             //Virtual function
{
  cerr<<"Wrong change_rate\n";
  return(FALSE);
}


double Like_model::prob_nonsyn(int start[3], int end[3], Branch *taxa, int rate_num)    //Virtual function
{
  cerr<<"Wrong prob_nonsyn\n";
  return(0);
}
//End of protected virtual functions




void Like_model::copy_transprobs(Branch *from, Branch *to, int rate_num)
{
  int i,j;

  for (i=0; i<curr_exchange->get_condlike_size(); i++)
	for (j=0; j<curr_exchange->get_condlike_size(); j++)	
	      to->set_trpb(rate_num, i,j, from->get_trpb(rate_num, i,j));
} //End Like_model::copy_transprobs




long double Like_model::prob_w_rate (int locale, int rate_num)
  //This is the key likelihood calculating function.  It calls partial_prob_w_rate
  //to do the post-order tree transversal and then uses the root frequencies to
  //calculate the overall likelihood
{
  int i ,j;
  double retval=0;

  partial_prob_w_rate(locale, curr_tree->get_leftmost_tip(), rate_num, 
				  curr_tree->find_root()); 

  if (curr_tree->rooted_tree() ==FALSE) {
	  for (i=0; i<curr_exchange->get_condlike_size(); i++) {
#ifdef _OPEN_MP_VERSION_
		retval+=root_freq(i)*curr_tree->find_root()->get_cond_prob_locale(omp_get_thread_num(), i);
#else
		retval+=root_freq(i)*curr_tree->find_root()->get_cond_prob(i);
#endif
		//  cout<<locale<<": "<<i<<": "<<curr_tree->find_root()->get_cond_prob(i)<<" : "<<root_freq(i)<<endl;
	  }
  }
  else {
	for (i=0; i<curr_exchange->get_condlike_size(); i++)
		for (j=0; j<curr_exchange->get_condlike_size();j++)
#ifdef _OPEN_MP_VERSION_
			retval+=root_freq(j)*curr_tree->find_root()->get_trpb(rate_num, j, i)*
			curr_tree->find_root()->get_cond_prob_locale(omp_get_thread_num(), i);
#else
			retval+=root_freq(j)*curr_tree->find_root()->get_trpb(rate_num, j, i)*
			curr_tree->find_root()->get_cond_prob(i);
#endif
  }
 
 return(retval);
  
}  //End Like_model::prob_site_w_rate




void Like_model::partial_prob_w_rate(int locale, Branch *lsib, int rate_num, Branch *stop_id)
  //The heart of the likelihood calculation: uses the tranisition probablity matrices stored in
  //the tree to calculate the likelihood of part of that tree.  Uses post-order tree transversal
{
  int i, j, k;
  double temp_p1, temp_p2;
  Branch *rsib, *parent;
  
  rsib=lsib->get_sibling();
  parent=lsib->get_parent();
  
  if (lsib->is_tip()==FALSE)
    {
      if (rsib->is_tip()==FALSE) 
	{
	  partial_prob_w_rate(locale, curr_tree->find_left_tip(rsib), rate_num, rsib);
	  for(i=0; i<curr_exchange->get_condlike_size(); i++)
	    {
	      temp_p1=temp_p2=0;
	      for(j=0; j<curr_exchange->get_condlike_size(); j++)
		{
		  temp_p1+=lsib->get_cond_prob(j)*lsib->get_trpb(rate_num, i, j);
		  temp_p2+=rsib->get_cond_prob(j)*rsib->get_trpb(rate_num, i, j);
		} 
	      parent->set_cond_prob(i, temp_p1*temp_p2);
	    }   
	}
    
      else
	{
	  for(i=0; i<curr_exchange->get_condlike_size(); i++)
	    {
	      temp_p1=0;
	      for(j=0; j<curr_exchange->get_condlike_size(); j++)
		temp_p1+=lsib->get_cond_prob(j)*lsib->get_trpb(rate_num, i, j);
	      parent->set_cond_prob(i, temp_p1*
		rsib->get_trpb(rate_num, i, (*curr_data)[rsib->get_taxa_id()][locale])); 
	    }
	}
    }

  else
    {
      if (rsib->is_tip()==FALSE) 
	{
	  partial_prob_w_rate(locale, curr_tree->find_left_tip(rsib), rate_num, rsib);
	
	  for(i=0; i<curr_exchange->get_condlike_size(); i++)
	    {
	      temp_p2=0;
	      for(j=0; j<curr_exchange->get_condlike_size(); j++)
		temp_p2+=rsib->get_cond_prob(j)*rsib->get_trpb(rate_num, i, j);  
	      parent->set_cond_prob(i, temp_p2*
				    lsib->get_trpb(rate_num, i, (*curr_data)[lsib->get_taxa_id()][locale]));    
	    }   
	}

      else
	  for(i=0; i<curr_exchange->get_condlike_size(); i++)
	      parent->set_cond_prob(i,
		lsib->get_trpb(rate_num, i, (*curr_data)[lsib->get_taxa_id()][locale])*
		rsib->get_trpb(rate_num, i, (*curr_data)[rsib->get_taxa_id()][locale])); 
	  
    }

  if (parent!=stop_id)
      partial_prob_w_rate(locale, parent, rate_num, stop_id);
 
}  //End Like_model::partial_prob_w_rate






void Like_model::array_copy(int array1[], int array2[], int size)
{
  int i;

  for(i=0; i<size; i++)
    array2[i]=array1[i];
}  //End Like_model::array_copy



int Like_model::num_rate_prob_params()
//Returns the number of parameters to optimize if we are allowing
//different rates of evolution and different probabilities of those rates
{
	if (curr_exchange->are_optimizing_rate_props() == TRUE)
		return(lookup_p_tree->get_num_cat_params());
	else
		return(0);
}

int Like_model::initialize_rate_prob_params(int start_param_index, double *par, PARAM_TYPE *types)
{
	int i;

	rate_prob_param_index=start_param_index;

	for(i=0; i<lookup_p_tree->get_num_cat_params(); i++) {
		par[start_param_index]=log(lookup_p_tree->get_param_value(i));
		types[start_param_index++]=RATE_PROB;
	}
	return(start_param_index);
}



//End of Like_model functions





//JC_model class functions
//JC_model::Public functions


double JC_model::get_expect_sub_site(double brlen)
{
 if (fabs(brlen)>FLOAT_TOL)
    return((3.0/4.0)*brlen);
  else
    return(0.0); 
}




void JC_model::print_model()
{
  cout<<"\n\nThe current model is the Jukes-Cantor (1969) model. (Equal base frequencies"<<endl
      <<"\tno transition/transversion rate difference.)"<<endl
      <<"\nNumber of Parameters: "<<curr_exchange->get_num_params()<<endl<<endl<<endl;
}  //End JC_model::print_model



 
void JC_model::describe_nuc_results()
{
  int i;

  cout<<"\nResults of Analysis:\n";
  
  for (i=0; i<curr_exchange->get_num_branches(); i++) 
	  cout<<"Branch "<<i<<": Length = "<<(*curr_tree)[i]->get_brnlen()<<endl;
  
}  //End JC_model::describe_results




//JC_model Protected functions
double JC_model::get_basefreq(int base, Branch *tree_branch)
{
  return(0.25);
}

double JC_model::get_basefreq(int base, int codon_pos, Branch *tree_branch)
{
  return(0.25);
}

double JC_model::get_trs_trv()
{
  return(1.0);
}


double JC_model::q_matrix_entry (int start_base, int end_base, int codon_pos, Branch *tree_branch)
{
    return(0.25);
} //End JC_model::q_matrix_entry




double JC_model::find_ut(Branch *taxa)
{
  if (fabs(taxa->expect_subs_site())>FLOAT_TOL)
    return((4.0/3.0)*taxa->expect_subs_site());
  else
    return(0.0);
}



//End JC_model class private functions





//K2P_model class functions
//K2P_model::Public functions

double K2P_model::get_expect_sub_site(double brlen)
{
 if (fabs(brlen)>FLOAT_TOL)
    return(fabs(brlen)/(0.5+get_trs_trv()*0.25));
  else
    return(0.0); 
}





void K2P_model::print_model()
{
  cout<<"\n\nThe current model is the Kimura 2 Parameter model. (Equal base frequencies"<<endl
      <<"\nwith a transition/transversion rate difference.)"<<endl
      <<"\nNumber of Parameters: "<<curr_exchange->get_num_params()<<endl;
  if(get_trs_trv()==1)
    cout<<"\nNote: trs/trv ratio set to 1--this model reduces to the Juke-Cantor model with these settings."<<endl;
  else
    cout<<"\nTrs/trv ratio (kappa)= "<<get_trs_trv()<<endl;
  cout<<endl<<endl;
}  //End K2P_model::print_model

 


void K2P_model::describe_nuc_results()
{
  int i;

  cout<<"\nResults of Analysis:\n";
  
  for (i=0; i<curr_exchange->get_num_branches(); i++)
	  cout<<"Branch "<<i<<": Length = "<<(*curr_tree)[i]->get_brnlen()<<endl;
  
  
  cout<<"Kappa parameter: "<<get_trs_trv()<<endl;

}  //End K2P_model::describe_results




//K2P_model Protected functions
double K2P_model::get_basefreq(int base, Branch *tree_branch)
{
  return(0.25);
}

double K2P_model::get_basefreq(int base, int codon_pos, Branch *tree_branch)
{
  return(0.25);
}

double K2P_model::get_trs_trv()
{
  return(curr_exchange->get_trs_trv());
}



double K2P_model::q_matrix_entry (int start_base, int end_base, int codon_pos, Branch *tree_branch)
{
  if (start_base==end_base)
    return(0.25);
  else
    if (is_transition(start_base, end_base)==TRUE)
      return(0.25*get_trs_trv());
    else
      return(0.25);
} //End K2P_model::q_matrix_entry




double K2P_model::find_ut(Branch *taxa)
{
  if (fabs(taxa->expect_subs_site())>FLOAT_TOL)
    return(fabs(taxa->expect_subs_site())*(0.5+get_trs_trv()*0.25));
  else
    return(0.0);
}
//End K2P Functions





//HKY_model class functions
//HKY_model::Public functions

double HKY_model::get_expect_sub_site(double brlen)
{

  double PurPyrfreqs[4];
 

  PurPyrfreqs[0]=PurPyrfreqs[2]=get_basefreq(0,(*curr_tree)[0])+get_basefreq(2,(*curr_tree)[0]);
  PurPyrfreqs[1]=PurPyrfreqs[3]=get_basefreq(1,(*curr_tree)[0])+get_basefreq(3,(*curr_tree)[0]);
 
  if (fabs(brlen)>FLOAT_TOL)
    return(fabs(brlen)*
	   (2*(PurPyrfreqs[0]*PurPyrfreqs[1]+get_trs_trv()*
	       (get_basefreq(0,(*curr_tree)[0])*get_basefreq(2,(*curr_tree)[0])+
		get_basefreq(1,(*curr_tree)[0])*get_basefreq(3,(*curr_tree)[0])))));
  else
    return(0.0);

}



void HKY_model::print_model()
{
  int i;
  BOOL all_same_freqs=TRUE;

  for(i=0; i<4; i++)
    {
      if (get_basefreq(i,0)!=0.25)
	all_same_freqs=FALSE;
    }

  cout<<"\n\nThe current model is the Hasegawa-Kishino-Yano (1985) model. (Unequal base frequencies"<<endl
      <<"\twith a transition/transversion rate difference.)"<<endl
      <<"\nNumber of Parameters: "<<curr_exchange->get_num_params()<<endl;
  if (all_same_freqs==TRUE)
    {
      if (get_trs_trv()!=1)
	cout<<"\nNote: Base frequencies are all equal--this model reduces to the"<<endl
	    <<"\tKimura 2 parameter model with these settings.  Trs/trv ratio (kappa)= "<<get_trs_trv()<<endl;

      else
	cout<<"\nNote: Trs/trv=1 and base frequencies are all equal--this model reduces to the"<<endl
	    <<"\tJukes-Cantor model with these settings."<<endl;
    }
  else
    {
      cout<<"\nTrs/trv ratio (kappa) = "<<get_trs_trv();
      cout<<".  Base frequencies:  A: "<<get_basefreq(0,0)<<" C: "
	  <<get_basefreq(1,0)
	  <<" G: "<<get_basefreq(2,0)<<" T: "<<get_basefreq(3,0)<<endl;
    }
  cout<<endl<<endl;
}  //End HKY_model::print_model



 
void HKY_model::describe_nuc_results()
{
  int i;

  cout<<"\nResults of Analysis:\n";
  
  for (i=0; i<curr_exchange->get_num_branches(); i++)
    cout<<"Branch "<<i<<": Length = "<<(*curr_tree)[i]->get_brnlen()<<endl;
	
  cout<<"Base frequencies: "<<endl;
  for (i=0; i<4; i++)
    cout<<get_basefreq(i,0)<<" ";

  cout<<"Kappa parameter: "<<get_trs_trv()<<endl;

}  //End HKY_model::describe_results




//HKY_model Protected functions
double HKY_model::get_basefreq(int base, Branch *tree_branch)
{
  return(curr_exchange->return_basefreq(base));
}


double HKY_model::get_basefreq(int base, int codon_pos, Branch *tree_branch)
{
  if (curr_exchange->using_codon_basefreqs()==TRUE)
    return(curr_exchange->return_codon_basefreq(codon_pos, base));
  else
    return(curr_exchange->return_basefreq(base));
}


double HKY_model::get_trs_trv()
{
  return(curr_exchange->get_trs_trv());
}


double HKY_model::q_matrix_entry (int start_base, int end_base, int codon_pos, Branch *tree_branch)
{ 
  if (start_base==end_base)
	return(get_basefreq(end_base, codon_pos, tree_branch));
  else {
    switch (start_base)
      {
      case(0):
        if (is_transition(start_base, end_base)==TRUE)
          return(get_basefreq(2,codon_pos, tree_branch)*get_trs_trv());
        else
          return(get_basefreq(end_base, codon_pos, tree_branch));
        break;
      case(1):
        if (is_transition(start_base, end_base)==TRUE)
          return(get_basefreq(3, codon_pos, tree_branch)*get_trs_trv());
        else
          return(get_basefreq(end_base, codon_pos, tree_branch));
        break;
      case(2):
        if (is_transition(start_base, end_base)==TRUE)
          return(get_basefreq(0, codon_pos, tree_branch)*get_trs_trv());
        else
          return(get_basefreq(end_base,codon_pos, tree_branch));
        break;
      case(3):
        if (is_transition(start_base, end_base)==TRUE)
          return(get_basefreq(1,codon_pos, tree_branch)*get_trs_trv());
        else
          return(get_basefreq(end_base,codon_pos,  tree_branch));
        break;
          default:
          return(0.0);
      }
  }
} //End HKY_model::q_matrix_entry




double HKY_model::find_ut(Branch *taxa)
{
 double PurPyrfreqs[4];
 

  PurPyrfreqs[0]=PurPyrfreqs[2]=get_basefreq(0,taxa)+get_basefreq(2,taxa);
  PurPyrfreqs[1]=PurPyrfreqs[3]=get_basefreq(1,taxa)+get_basefreq(3,taxa);
 
  if (fabs(taxa->expect_subs_site())>FLOAT_TOL)
    return(fabs(taxa->expect_subs_site())/
	   (2*(PurPyrfreqs[0]*PurPyrfreqs[1]+get_trs_trv()*
	       (get_basefreq(0,taxa)*get_basefreq(2,taxa)+get_basefreq(1,taxa)*get_basefreq(3,taxa)))));
  else
    return(0.0);
}


//GG_98__model class functions
//GG_98_model::Public functions


void GG_98_model::print_model()
{
  int i, j;
  BOOL all_same_freqs=TRUE;

  for(j=0; j<curr_exchange->get_num_branches(); j++)
    for(i=0; i<4; i++)
      {
	if (get_basefreq(i,(*curr_tree)[j])!=0.25)
	  all_same_freqs=FALSE;
    }

  cout<<"\n\nThe current model is the Galtier and Gaut (1998) model. (Base frequencies vary among branches"<<endl
      <<"\twith a transition/transversion rate difference.)"<<endl;
  if (all_same_freqs==TRUE)
    {
      if (get_trs_trv()!=1)
	cout<<"\nNote: Base frequencies are all equal--this model reduces to the"<<endl
	    <<"\tKimura 2 parameter model with these settings.  Trs/trv ratio (kappa)= "<<get_trs_trv()<<endl;

      else
	cout<<"\nNote: Trs/trv=1 and base frequencies are all equal--this model reduces to the"<<endl
	    <<"\tJukes-Cantor model with these settings."<<endl;
    }
  else
    {
      cout<<"\nTrs/trv ratio (kappa) = "<<get_trs_trv();
    }
  cout<<endl<<endl;
}  //End HKY_model::print_model



 
void GG_98_model::describe_nuc_results()
{
  int i;

  cout<<"\nResults of Analysis:\n";
  
  for (i=0; i<curr_exchange->get_num_branches(); i++)
    {
      cout<<"Branch "<<i<<": Length = "<<(*curr_tree)[i]->get_brnlen()<<endl;
      cout<<"Base frequencies A: "<<get_basefreq(0, (*curr_tree)[i])<<" C: "<<get_basefreq(1, (*curr_tree)[i])
	  <<" G: "<<get_basefreq(2, (*curr_tree)[i])<<" T: "<<get_basefreq(3, (*curr_tree)[i])<<endl;
    }

  cout<<"Kappa parameter: "<<get_trs_trv()<<endl;

}  //End GG_98_model::describe_results




//GG_98_model Protected functions
double GG_98_model::get_basefreq(int base, Branch *tree_branch)
{
  return(tree_branch->get_branch_basefreq(base));
}


double GG_98_model::get_basefreq(int base, int codon_pos, Branch *tree_branch)
{
  if (curr_exchange->using_codon_basefreqs()==TRUE)
    return(curr_exchange->return_codon_basefreq(codon_pos, base));
  else
    return(tree_branch->get_branch_basefreq(base));
}


double GG_98_model::get_trs_trv()
{
  return(curr_exchange->get_trs_trv());
}


double GG_98_model::q_matrix_entry (int start_base, int end_base, int codon_pos, Branch *tree_branch)
{ 
  switch (start_base)
    {
    case(0):
        if (start_base==end_base)
            return(1-(get_basefreq(2,tree_branch)*get_trs_trv()+
                    (get_basefreq(1,tree_branch)+get_basefreq(3,tree_branch))));
        else {
            if (is_transition(start_base, end_base)==TRUE)
                return(get_basefreq(2,tree_branch)*get_trs_trv());
            else {
                if (end_base==1)
                    return(get_basefreq(1, tree_branch));
                else
                    return(get_basefreq(3, tree_branch));
            }
        }
    case(1):
        if (start_base==end_base)
            return(1-(get_basefreq(3,tree_branch)*get_trs_trv()+
                      (get_basefreq(0,tree_branch)+get_basefreq(2,tree_branch))));
        else {
            if (is_transition(start_base, end_base)==TRUE)
                return(get_basefreq(3,tree_branch)*get_trs_trv());
            else {
                if (end_base==0)
                    return(get_basefreq(0,tree_branch));
                else
                    return(get_basefreq(2,tree_branch));
            }
        }
    case(2):
        if (start_base==end_base)
            return(1-(get_basefreq(0,tree_branch)*get_trs_trv()+
                      (get_basefreq(1,tree_branch)+get_basefreq(3,tree_branch))));
        else {
            if (is_transition(start_base, end_base)==TRUE)
                return(get_basefreq(0,tree_branch)*get_trs_trv());
            else {
                if (end_base==1)
                    return(get_basefreq(1,tree_branch));
                else
                    return(get_basefreq(3,tree_branch));
            }
        }
    case(3):
        if (start_base==end_base)
            return(1-(get_basefreq(1,tree_branch)*get_trs_trv()+
                  (get_basefreq(0,tree_branch)+get_basefreq(2,tree_branch))));
        else {
            if (is_transition(start_base, end_base)==TRUE)
                return(get_basefreq(1,tree_branch)*get_trs_trv());
            else {
                if (end_base==0)
                    return(get_basefreq(0,tree_branch));
                else
                    return(get_basefreq(2,tree_branch));
            }
        }
    default:
        return(0.0);
        break;
    }
 
} //End GG_98_model::q_matrix_entry




double GG_98_model::find_ut(Branch *taxa)
{
 double PurPyrfreqs[4];
 

  PurPyrfreqs[0]=PurPyrfreqs[2]=get_basefreq(0, taxa)+get_basefreq(2, taxa);
  PurPyrfreqs[1]=PurPyrfreqs[3]=get_basefreq(1, taxa)+get_basefreq(3, taxa);
 
  if (fabs(taxa->expect_subs_site())>FLOAT_TOL)
    return(fabs(taxa->expect_subs_site())/
	   (2*(PurPyrfreqs[0]*PurPyrfreqs[1]+get_trs_trv()*
	       (get_basefreq(0,taxa)*get_basefreq(2,taxa)+get_basefreq(1,taxa)*get_basefreq(3,taxa)))));
  else
    return(0.0);
}



