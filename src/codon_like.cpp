//Copyright 1999-2002 Gavin Conant

#include <iostream>
#include <math.h>
#include <iomanip>
#include "codon_like.h"
#include "powell.h"
#include "nrutil.h"

using namespace::std;

//Protected Codon_model

void Codon_model::initialize_arrays()
  //This function is used by all codon models to initialize the q_matrix array
  //which is used for interaction with the Linear_algebra class and ultimately,
  //the lapack routines
{ 
  int i;
 
  q_matrix=new double *[curr_exchange->get_condlike_size()];
  t_matrix=new double *[curr_exchange->get_condlike_size()];
  
  for(i=0; i<curr_exchange->get_condlike_size(); i++) {
    q_matrix[i]=new double [curr_exchange->get_condlike_size()];
    t_matrix[i]=new double [curr_exchange->get_condlike_size()];
  }

   if(curr_exchange->branch_basefreqs()==FALSE)
       calc_codon_matrices(0,(*curr_tree)[0]);
   
  else
     for (i=0; i<curr_exchange->get_num_branches(); i++)
       {
	 calc_codon_matrices(0,(*curr_tree)[i]);
	 calc_transprobs((*curr_tree)[i], 0);
       }
}

void Codon_model::set_expect_subs()
{
	int i, j, k, l, start_codon, end_codon, start[3], end[3];
	double val, val2, freq, diag_sum, diag_sum2, so_far, so_far2, stop_freq;

	stop_freq=0;

	for (start_codon=0; start_codon<curr_exchange->get_condlike_size(); start_codon++)
	{ 
		if (curr_code->is_stop(start_codon)==TRUE) {
			for (j=0; j<3; j++)
				start[j]=curr_code->get_pos_n(start_codon, j);

			stop_freq+=get_basefreq(start[0], 0, (*curr_tree)[0])*
			get_basefreq(start[1], 1, (*curr_tree)[0])*get_basefreq(start[2], 2, (*curr_tree)[0]);
		}   
	}
	stop_excess=1/(1-stop_freq);


	for(i=0; i<curr_exchange->get_num_branches(); i++) {      
		prop_syn=prop_nsyn=neu_prop_syn=neu_prop_nsyn=0;
		diag_sum=diag_sum2=0;
		for (start_codon=0; start_codon<curr_exchange->get_condlike_size(); start_codon++) {
			if (curr_code->is_stop(start_codon)!=TRUE) { 
				so_far=so_far2=0;
				
				for (j=0; j<3; j++)
					start[j]=curr_code->get_pos_n(start_codon, j);

				freq= get_basefreq(start[0], 0, (*curr_tree)[0])*
					get_basefreq(start[1], 1, (*curr_tree)[0])*
					get_basefreq(start[2], 2, (*curr_tree)[0])*stop_excess;

				for(j=0; j<3; j++) {
					for (k=0; k<3; k++) {
						end_codon=curr_code->one_diff(start_codon, j, k);

						for (l=0; l<3; l++)
							end[l]=curr_code->get_pos_n(end_codon, l);

						if (curr_code->is_stop(end_codon)!=TRUE) {
							val=q_matrix_entry(start[0], end[0], 0, (*curr_tree)[0])* 
								q_matrix_entry(start[1], end[1], 1, (*curr_tree)[0])*
								q_matrix_entry(start[2], end[2], 2, (*curr_tree)[0])*
								prob_nonsyn(start, end,(*curr_tree)[i], 0)*stop_excess;
							
							val2=	q_matrix_entry(start[0], end[0], 0, (*curr_tree)[0])* 
								q_matrix_entry(start[1], end[1], 1, (*curr_tree)[0])*
								q_matrix_entry(start[2], end[2], 2, (*curr_tree)[0])*stop_excess;

							if (curr_code->is_non_synon(start, end)==FALSE) { 
								prop_syn+=freq*val;
								neu_prop_syn+=freq*val2;
							}

							so_far +=val;
							so_far2+=val2;
						}
					}
				}
				diag_sum+=freq*so_far;
				diag_sum2+=freq*so_far2;
			}    
		}


		prop_syn/=diag_sum;
		neu_prop_syn/=diag_sum2;
		prop_nsyn=1.0-prop_syn;
		neu_prop_nsyn=1.0-neu_prop_syn;

		// cout<<"Prop syn: "<<prop_syn<<" neu : "<<neu_prop_syn<<endl;
		if ((*curr_tree)[i]->get_brnlen()> MIN_BRLEN) {
			(*curr_tree)[i]->set_expect_subs_site((*curr_tree)[i]->get_brnlen()*prop_syn/(3.0*neu_prop_syn));
			(*curr_tree)[i]->set_expect_nsyn_subs_site((*curr_tree)[i]->get_brnlen()*prop_nsyn/(3.0*neu_prop_nsyn));
		}
		else {
			(*curr_tree)[i]->set_expect_subs_site(0);
			(*curr_tree)[i]->set_expect_nsyn_subs_site(0);
		}
        
        (*curr_tree)[i]->name_this_branch();
        cout<<(*curr_tree)[i]->get_name()<<"\t"<<(*curr_tree)[i]->expect_subs_site()<<endl;
	}
  

  
}


double Codon_model::get_obs_trs_trv_ratio()
{
  int i, j, k, l, start_codon, start[3], end_codon, end[3], diff_a, diff_b;
  double freq, val, transition_prob=0, transversion_prob=0, stop_freq;


  for (start_codon=0; start_codon<curr_exchange->get_condlike_size(); start_codon++) {   
    if (curr_code->is_stop(start_codon)==TRUE) {
      for (j=0; j<3; j++)
	start[j]=curr_code->get_pos_n(start_codon, j);
      stop_freq+=get_basefreq(start[0], 0, (*curr_tree)[0])*get_basefreq(start[1], 1, (*curr_tree)[0])*get_basefreq(start[2], 2, (*curr_tree)[0]);
    }  
  }

  stop_excess=1/(1-stop_freq);

  for (start_codon=0; start_codon<curr_exchange->get_condlike_size(); start_codon++)
    {
      if (curr_code->is_stop(start_codon)!=TRUE)
	{
	  for (j=0; j<3; j++)
	    start[j]=curr_code->get_pos_n(start_codon, j);
	  
	  freq= get_basefreq(start[0], 0, (*curr_tree)[0])*
	    get_basefreq(start[1], 1, (*curr_tree)[0])*
	    get_basefreq(start[2], 2, (*curr_tree)[0])*stop_excess;
	  
	  for(j=0; j<3; j++)
	    {
	      for (k=0; k<3; k++)
		{
		  end_codon=curr_code->one_diff(start_codon, j, k);
		  
		  for (i=0; i<3; i++)
		    end[i]=curr_code->get_pos_n(end_codon, i);
		  
		  if (curr_code->is_stop(end_codon)!=TRUE) {
		    val=q_matrix_entry(start[0], end[0], 0, (*curr_tree)[0])* 
		      q_matrix_entry(start[1], end[1], 1, (*curr_tree)[0])*
		      q_matrix_entry(start[2], end[2], 2, (*curr_tree)[0])*
		      prob_nonsyn(start, end,(*curr_tree)[0], 0)*stop_excess;
		    
		    for(l=0; l<3; l++)
		      {
			if (end[l] != start[l])
			  {
			    diff_a=start[l];
			    diff_b=end[l];
			  }
		      }

		    if (is_transition(diff_a, diff_b) ==TRUE)
		      transition_prob += val*freq;
		    else
		      transversion_prob += val*freq;
		    
		  }
		}
	    }
	
	}    
    }
  return(transition_prob/transversion_prob);
	
}



void Codon_model::set_ratios(Branch *taxa)
{
  int i, j, k, l, start_codon, end_codon, start[3], end[3];
  double val, val2, freq, diag_sum, diag_sum2, so_far, so_far2, stop_freq;
 
  stop_freq=0;
  
  for (start_codon=0; start_codon<curr_exchange->get_condlike_size(); start_codon++)
    { 
      if (curr_code->is_stop(start_codon)==TRUE) {
	for (j=0; j<3; j++)
	  start[j]=curr_code->get_pos_n(start_codon, j);
	    
	stop_freq+=get_basefreq(start[0], 0, (*curr_tree)[0])*
	  get_basefreq(start[1], 1, (*curr_tree)[0])*get_basefreq(start[2], 2, (*curr_tree)[0]);
      }   
    }
  stop_excess=1/(1-stop_freq);
  

      
  prop_syn=prop_nsyn=neu_prop_syn=neu_prop_nsyn=0;
  diag_sum=diag_sum2=0;
  for (start_codon=0; start_codon<curr_exchange->get_condlike_size(); start_codon++)
    {
      if (curr_code->is_stop(start_codon)!=TRUE)
	{ 
	  so_far=so_far2=0;
	  for (j=0; j<3; j++)
	    start[j]=curr_code->get_pos_n(start_codon, j);
	  
	      freq= get_basefreq(start[0], 0, (*curr_tree)[0])*
		get_basefreq(start[1], 1, (*curr_tree)[0])*
		get_basefreq(start[2], 2, (*curr_tree)[0])*stop_excess;
	      
	      for(j=0; j<3; j++)
		{
		  for (k=0; k<3; k++)
		    {
		      end_codon=curr_code->one_diff(start_codon, j, k);
		      
		      for (l=0; l<3; l++)
			end[l]=curr_code->get_pos_n(end_codon, l);
		      
		      if (curr_code->is_stop(end_codon)!=TRUE) {
			val=q_matrix_entry(start[0], end[0], 0, (*curr_tree)[0])* 
			  q_matrix_entry(start[1], end[1], 1, (*curr_tree)[0])*
			  q_matrix_entry(start[2], end[2], 2, (*curr_tree)[0])*
			  prob_nonsyn(start, end, taxa, 0)*stop_excess;
			val2=	q_matrix_entry(start[0], end[0], 0, (*curr_tree)[0])* 
			  q_matrix_entry(start[1], end[1], 1, (*curr_tree)[0])*
			  q_matrix_entry(start[2], end[2], 2, (*curr_tree)[0])*stop_excess;
			
			if (curr_code->is_non_synon(start, end)==FALSE)
			  { 
			    prop_syn+=freq*val;
			    neu_prop_syn+=freq*val2;
			  }
			
			so_far +=val;
			so_far2+=val2;
		      }
		    }
		}
	      diag_sum+=freq*so_far;
	      diag_sum2+=freq*so_far2;
	}    
    }
  
  
  prop_syn/=diag_sum;
  neu_prop_syn/=diag_sum2;
  prop_nsyn=1.0-prop_syn;
  neu_prop_nsyn=1.0-neu_prop_syn;
  //Is this wrong??? GCC 6/11/2004
  //prop_syn/(3.0*neu_prop_syn);


  taxa->set_syn_ratio(prop_syn/(3.0*neu_prop_syn));
  taxa->set_nonsyn_ratio(prop_nsyn/(3.0*neu_prop_nsyn));

}



double Codon_model::root_freq(int site)
  //Used to calculate likelihood under codon-models--this function
  //specifys the probablity of each codon in the model
{
	return(stop_excess*get_basefreq(curr_code->get_pos_n(site,0), 0, curr_tree->find_root())*
			   get_basefreq(curr_code->get_pos_n(site,1), 1, curr_tree->find_root())*
			   get_basefreq(curr_code->get_pos_n(site,2), 2, curr_tree->find_root())); 
} //End Codon_model::root_freq



Codon_model::~Codon_model()
{
  int i;

  for (i=0; i<curr_exchange->get_condlike_size(); i++){
    delete[] q_matrix[i];
    delete[] t_matrix[i];
  }
  delete[] q_matrix;
  delete[] t_matrix;
}



void Codon_model::calc_codon_matrices (int rate_num, Branch *taxa)
  //Creates the instanteous transition probablity matrix to be evaluated by linear
  //algebra routines
{
  int i, j, k, start_codon, end_codon, start[3], end[3];
  double so_far, diag_sum=0, stop_freq=0;

  for (start_codon=0; start_codon<curr_exchange->get_condlike_size(); start_codon++) {   
    if (curr_code->is_stop(start_codon)==TRUE) {
      for (j=0; j<3; j++)
	start[j]=curr_code->get_pos_n(start_codon, j);
      stop_freq+=get_basefreq(start[0], 0, taxa)*get_basefreq(start[1], 1, taxa)*get_basefreq(start[2], 2, taxa);
    }  
  }

  stop_excess=1/(1-stop_freq);
  
  for (start_codon=0; start_codon<curr_exchange->get_condlike_size(); start_codon++) {
      if (curr_code->is_stop(start_codon)!=TRUE) {
          so_far=0;
        
          for (j=0; j<curr_exchange->get_condlike_size(); j++)
            q_matrix[start_codon][j]=0.0;
          
          for (j=0; j<3; j++)
            start[j]=curr_code->get_pos_n(start_codon, j);
          
          for(j=0; j<3; j++) {
              for (k=0; k<3; k++) {
                  end_codon=curr_code->one_diff(start_codon, j, k);
              
                  for (i=0; i<3; i++)
                      end[i]=curr_code->get_pos_n(end_codon, i);
              
                  if (curr_code->is_stop(end_codon)!=TRUE) {
                      q_matrix[start_codon][end_codon]=
                      q_matrix_entry(start[0], end[0], 0, taxa)*
                      q_matrix_entry(start[1], end[1], 1, taxa)*
                      q_matrix_entry(start[2], end[2], 2, taxa)*stop_excess*
                      prob_nonsyn(start, end, taxa, rate_num);
              
             
                      so_far+=q_matrix[start_codon][end_codon];
                  }
              }
          }
          q_matrix[start_codon][start_codon]=-so_far;
        
          diag_sum+=(so_far*stop_excess*get_basefreq(start[0], 0, taxa)*get_basefreq(start[1], 1, taxa)*get_basefreq(start[2], 2, taxa));
      }
      else {
          for (j=0; j<curr_exchange->get_condlike_size(); j++) q_matrix[start_codon][j]=0;
      }
  }

  
  
  for (i=0; i<curr_exchange->get_condlike_size(); i++) {
    for (j=0; j<curr_exchange->get_condlike_size(); j++)
      if (q_matrix[i][j] != 0)
	  q_matrix[i][j]/=diag_sum;
      
    
  }
  
 curr_lin_alg->find_eigen_matrix(q_matrix);
 last_codonmat=taxa->get_nonsyn_pattern();
}  //End Codon_model::calc_codon_matrices




void Codon_model::calc_transprobs (Branch *taxa, int rate_num)
  //Once we have used calc_codon_matrices matrices, the transition probablity matrix
  //for every branch can be calculated with 2 matrix-matrix multiplications
{
  int i,j;
  

  //Sets Transprobs for gaps to 1
  for (i=0; i<curr_exchange->get_condlike_size()+1; i++)
    {
      taxa->set_trpb(rate_num, i, curr_exchange->get_condlike_size(), 1.0);
      taxa->set_trpb(rate_num, curr_exchange->get_condlike_size(), i, 1.0);
    }

 
  curr_lin_alg->find_scaled_matrix(t_matrix, taxa->get_brnlen(), 0);

 
  if (fabs(taxa->get_brnlen())>MIN_BRLEN)
    {
      for (i=0; i<curr_exchange->get_condlike_size(); i++)
	for (j=0; j<curr_exchange->get_condlike_size(); j++)
	  {	
	    if (curr_code->is_stop(i)==FALSE && curr_code->is_stop(j)==FALSE)
		taxa->set_trpb(rate_num, i, j, t_matrix[i][j]);
				  
	    else
	      taxa->set_trpb(rate_num, i, j, 0.0);
	  }
    }
  else
    {
      for (i=0; i<curr_exchange->get_condlike_size(); i++)
	for (j=0; j<curr_exchange->get_condlike_size(); j++)
	  {	
	    if (i==j) 
	      taxa->set_trpb(rate_num, i, j, 1.0);
	    else
	      taxa->set_trpb(rate_num, i, j, 0.0);
	  }
    }
  
  
}  //End Codon_model::calc_transprobs


void Codon_model::calc_first_derivatives(Branch *taxa, double ut)
{
	calc_first_derivatives(taxa, ut, 0);
}


void Codon_model::calc_first_derivatives(Branch *taxa, double ut, int rate)
  //Following Yang 2000, this function calculates the first deriviative
  //of a transition probablity matrix
{
 int i,j;
  

 //t_matrix reused to pass ut-dependant derivative of transprobs matrix
 curr_lin_alg->find_scaled_matrix(t_matrix, taxa->get_brnlen(), 1);
 
 for (i=0; i<curr_exchange->get_condlike_size(); i++)
   for (j=0; j<curr_exchange->get_condlike_size(); j++)
     {	
       if (curr_code->is_stop(i)==FALSE && curr_code->is_stop(j)==FALSE)
	 taxa->set_trpb_prime(i, j, rate, t_matrix[i][j]);
       
       else
	 taxa->set_trpb_prime(i, j, rate, 0.0);
     }
}

void Codon_model::calc_second_derivatives(Branch *taxa, double ut) 
{
	calc_second_derivatives(taxa, ut, 0);
}
void Codon_model::calc_second_derivatives(Branch *taxa, double ut, int rate) 
  //Following Yang 2000, this function calculates the second deriviative
  //of a transition probablity matrix
{
int i,j;
  

//t_matrix reused to pass ut-dependant derivative of transprobs matrix
 curr_lin_alg->find_scaled_matrix(t_matrix, taxa->get_brnlen(), 2);
 
 for (i=0; i<curr_exchange->get_condlike_size(); i++)
   for (j=0; j<curr_exchange->get_condlike_size(); j++)
     {	
       if (curr_code->is_stop(i)==FALSE && curr_code->is_stop(j)==FALSE)
	 taxa->set_trpb_double_prime(i, j, rate, t_matrix[i][j]);
       
       else
	 taxa->set_trpb_double_prime(i, j, rate, 0.0);
     }
}




//MG_94_model functions

void MG_94_model::describe_results()
{
  int i;
  describe_nuc_results();
  cout<<"Instantaneous probablities of a non-synonymous substitutions:\n";
    for(i=0; i<curr_exchange->get_num_p_non_syn(); i++)
      cout<<i<<": "<<curr_exchange->get_p_non_syn(i)<<endl;
}


int MG_94_model::rate_param_size()
{
  return(1);
}


//Protected MG_94_model
BOOL MG_94_model::change_rate(int rate_num, double *rate_info)
{
  BOOL retval=FALSE;
  if (fabs(curr_exchange->get_p_non_syn(0, rate_num)-rate_info[0])>FLOAT_TOL)
    {
      curr_exchange->set_p_non_syn(0, rate_num, rate_info[0]);
      retval=TRUE;
    } 
  return(retval);
}



double MG_94_model::prob_nonsyn(int start[3], int end[3], Branch *taxa,  int rate_num)
{
  return(integer_power(curr_exchange->get_p_non_syn(taxa->get_p_nonsyn_num(), rate_num), curr_code->is_non_synon(start, end)));
}  //End MG_94_model::prob_nonsyn

void MG_94_model::set_num_pns_used()
{
	int i,j;

	num_pns_coeff=0;

	for(i=0; i<curr_exchange->get_num_p_non_syn(); i++) {
		for(j=0; j<curr_exchange->get_num_rates(); j++) {
			if(curr_exchange->p_non_syn_fixed(i,j) ==FALSE)
				num_pns_coeff++;
		}
	}
}
  

int MG_94_model::set_pns_id_indices(double par[], PARAM_TYPE types[])
{
  int i, j, param_num=0;
  pns_start=0;

  prop_pns_index = new int [num_pns_coeff];
  prop_rate_index = new int [num_pns_coeff];
  
  for(i=0; i<curr_exchange->get_num_rates(); i++) {
	 for(j=0; j<curr_exchange->get_num_p_non_syn(); j++) {
		if ((curr_exchange->is_p_non_syn_used(j) == TRUE) && (curr_exchange->p_non_syn_fixed(j,i) ==FALSE)) {
			par[param_num]=log(curr_exchange->get_p_non_syn(j,i));
			types[param_num]=SUB_TYPE;
			prop_rate_index[param_num]=i;
			prop_pns_index[param_num++]=j;
		}
	}
  }
  
  return(param_num);	
}


//C_00_model functions


void C_00_model::describe_results()
{
  describe_nuc_results();
  cout<<"Instantaneous probablity of a non-synonymous extra group substitution: "
      <<curr_exchange->get_p_non_syn(0, 0)<<endl;
  cout<<"Instantaneous probablity of a non-synonymous within group substitution: "
      <<curr_exchange->get_p_inter_group(0)<<endl;
}


int C_00_model::rate_param_size()
{
  return(2);
}


//Protected C_00_model
BOOL C_00_model::change_rate(int rate_num, double *rate_info)
{
  BOOL retval=FALSE;
  if (fabs(curr_exchange->get_p_non_syn(0, rate_num)-rate_info[0])>FLOAT_TOL)
    {
      curr_exchange->set_p_non_syn(rate_num, rate_info[0]);
      retval=TRUE;
    }
  if (fabs(curr_exchange->get_p_inter_group(0, rate_num)-rate_info[1])>FLOAT_TOL)
    {
      curr_exchange->set_p_inter_group(rate_num, rate_info[1]);
      retval=TRUE;
    }
  return(retval);
}



double C_00_model::prob_nonsyn(int start[3], int end[3], Branch *taxa,  int rate_num)
{
  if (curr_exchange->curr_groups->get_group(readchar_to_aa(curr_code->return_AA(start)))==
      curr_exchange->curr_groups->get_group(readchar_to_aa(curr_code->return_AA(end))))
    return(integer_power(curr_exchange->get_p_inter_group(taxa->get_p_intergroup_num(), rate_num), curr_code->is_non_synon(start, end)));
  else
    return(integer_power(curr_exchange->get_p_non_syn(taxa->get_p_nonsyn_num(), rate_num), curr_code->is_non_synon(start, end)));
}



void C_00_model::set_num_pns_used()
{
	int i, j;

	num_pns_coeff=0;
	num_pi_coeff=0;
	
	for(i=0; i<curr_exchange->get_num_p_non_syn(); i++)
	{
		for(j=0; j<curr_exchange->get_num_rates(); j++) {
			if ((curr_exchange->is_p_non_syn_used(i) == TRUE) && (curr_exchange->p_non_syn_fixed(i, j) ==FALSE))
				num_pns_coeff++;
			if ((curr_exchange->is_p_inter_group_used(i) == TRUE)  && (curr_exchange->p_inter_group_fixed(i, j) ==FALSE))
				num_pi_coeff++;
		}
	}

}


int C_00_model::set_pns_id_indices(double par[], PARAM_TYPE types[])
{
	int i, j, param_num=0;
	pns_start=0;

	prop_pns_index = new int [num_pns_coeff+num_pi_coeff];
	prop_rate_index = new int [num_pi_coeff+num_pns_coeff];

	for(i=0; i<curr_exchange->get_num_rates(); i++) {
		for(j=0; j<curr_exchange->get_num_p_non_syn(); j++) {
			if ((curr_exchange->is_p_non_syn_used(j) == TRUE) && (curr_exchange->p_non_syn_fixed(j, i) ==FALSE)) {
				par[param_num]=log(curr_exchange->get_p_non_syn(j, i));
				types[param_num]=SUB_TYPE;
				prop_rate_index[param_num]=i;
				prop_pns_index[param_num++]=j;
			}
		
			if ((curr_exchange->is_p_inter_group_used(j) == TRUE)  && (curr_exchange->p_inter_group_fixed(j, i) ==FALSE)) {
				par[param_num]=log(curr_exchange->get_p_inter_group(j, i));
				types[param_num]=CLASS_SUB;
				prop_rate_index[param_num]=i;
				prop_pns_index[param_num++]=j;
			}
		}
	
	}
	return(param_num);	

}

//Public LCAP_model

void LCAP_model::describe_results()
{
  describe_nuc_results();
  cout<<"Chem. Sim. coefficent  : "<<curr_exchange->get_aa_prop_fac(0, CHEM_COMP)<<endl;
  cout<<"Polarity coefficent    : "<<curr_exchange->get_aa_prop_fac(0, POLARITY)<<endl;
  cout<<"Volume coefficent      : "<<curr_exchange->get_aa_prop_fac(0, VOLUME)<<endl;
  cout<<"Iso. Elec. coefficent  : "<<curr_exchange->get_aa_prop_fac(0, ISO_ELEC)<<endl;
  cout<<"Hydorpahty coefficent  : "<<curr_exchange->get_aa_prop_fac(0, HYDROPATHY)<<endl;
  cout<<"Residual coefficent    : "<<curr_exchange->get_aa_prop_fac(0, SCALING)<<endl;
}



int LCAP_model::rate_param_size()
{
  return(curr_exchange->get_num_aa_props());
}



//Protected LCAP_model
BOOL LCAP_model::change_rate(int rate_num, double *rate_info)
{
  int i;
  BOOL retval=FALSE;
  
  
  for(i=0; i<curr_exchange->get_num_aa_props(); i++) {
	if (fabs(curr_exchange->get_aa_prop_fac(rate_num, i)-
		rate_info[i])>FLOAT_TOL)
		{
			curr_exchange->set_aa_prop_fact(rate_num, i, 
				      rate_info[i]);
			retval=TRUE;
		} 
  }
  
  if (fabs(curr_exchange->get_aa_prop_fac(rate_num, SCALING)-
	   rate_info[curr_exchange->get_prop_index_num(SCALING)])>FLOAT_TOL)
    {
      curr_exchange->set_aa_prop_fact(rate_num, SCALING, 
				      rate_info[curr_exchange->get_prop_index_num(SCALING)]);
      retval=TRUE;
    } 
 
  return(retval);
}




double LCAP_model::prob_nonsyn(int start[3], int end[3], Branch *taxa, int rate_num)
{
	int i;
	double net_p=1;

	if(curr_code->is_non_synon(start, end)==TRUE) {
		for(i=0; i<curr_exchange->get_num_aa_props(); i++)
			net_p*=exp(curr_exchange->get_aa_prop_fac(taxa->get_aa_prop_num(i), rate_num, i)
				*curr_code->aa_pair_diff(start, end, i));
		net_p*=exp(curr_exchange->get_aa_prop_fac(taxa->get_aa_prop_num(SCALING), rate_num, SCALING));
	}





    return(net_p);
}

void LCAP_model::set_num_aa_props_used()
{
	int i, j, k;

	num_aa_props_per_rate=new int [curr_exchange->get_num_rates()];

	num_aa_props_per_rate[0]=0;
	for(k=0; k<curr_exchange->get_num_p_non_syn(); k++) {
		for(j=0; j<curr_exchange->get_num_aa_props(); j++) {
			if ((curr_exchange->aa_prop_allowed(j) == TRUE) && 
					(curr_exchange->is_aa_prop_used(k, j) == TRUE)) {
				num_aa_props_per_rate[0]++;
			}
		}
			if ((curr_exchange->aa_prop_allowed(SCALING) == TRUE) &&
					(curr_exchange->is_aa_prop_used(k, SCALING) == TRUE))
					num_aa_props_per_rate[0]++;
		
	}

	for(i=1; i<curr_exchange->get_num_rates(); i++)
	{
		num_aa_props_per_rate[i]=0;
		for(j=0; j<curr_exchange->get_num_aa_props()+1; j++)
		{
			for(k=0; k<curr_exchange->get_num_p_non_syn(); k++) {
				if ((curr_exchange->aa_prop_allowed(j) == TRUE) && 
					(curr_exchange->is_aa_prop_used(k, j) == TRUE)  &&
					(curr_exchange->get_rate_aa_prop_rate(i,j) != 0)) {
					num_aa_props_per_rate[i]++;
				}
			}
				
		}
	}

	num_aa_prop_coeff=0;
	
	for(i=0; i<curr_exchange->get_num_rates(); i++) 
		num_aa_prop_coeff+=num_aa_props_per_rate[i];
}


//CHECK THIS!!!!!!!
int LCAP_model::set_aa_prop_id_indices(double par[], PARAM_TYPE types[])
{
	int i, j, k, param_num=0;	

	pns_start=0;
	prop_index=new AA_PROPERTIES[num_aa_prop_coeff];
	prop_pns_index = new int [num_aa_prop_coeff];
	prop_rate_index = new int [num_aa_prop_coeff];
	

	//Rate 0
	for(j=0; j<curr_exchange->get_num_p_non_syn(); j++) {
			for(i=0; i<curr_exchange->get_num_aa_props(); i++)
			{
				if((curr_exchange->aa_prop_allowed(i) == TRUE) && 
					(curr_exchange->is_aa_prop_used(j, i) == TRUE)) {
					par[param_num]=curr_exchange->get_aa_prop_fac(j, 0, i);
					types[param_num]=AA_PROP_FAC;
					prop_index[param_num]=curr_exchange->get_prop_num_n(i);
					prop_rate_index[param_num]=0;
					prop_pns_index[param_num++]=j;
				}
			}

			if (curr_exchange->is_aa_prop_used(j, SCALING) == TRUE) {
				par[param_num]=curr_exchange->get_aa_prop_fac(j, 0, SCALING);
				types[param_num]=SCALING_FAC;
				prop_index[param_num]=SCALING;
				prop_rate_index[param_num]=0;
				prop_pns_index[param_num++]=j;
			}
		}

		for(k=1; k<curr_exchange->get_num_rates(); k++) {
			for(j=0; j<curr_exchange->get_num_p_non_syn(); j++) {
				for(i=0; i<curr_exchange->get_num_aa_props(); i++)
				{
					if((curr_exchange->aa_prop_allowed(i) == TRUE) && 
						(curr_exchange->is_aa_prop_used(j, i) == TRUE) &&
						(curr_exchange->get_rate_aa_prop_rate(k,i)!= 0)) {
						par[param_num]=curr_exchange->get_aa_prop_fac(j, k, i);
						types[param_num]=AA_PROP_FAC;
						prop_index[param_num]=curr_exchange->get_prop_num_n(i);
						prop_rate_index[param_num]=k;
						prop_pns_index[param_num++]=j;
					}
				}

				if (curr_exchange->is_aa_prop_used(j, SCALING) == TRUE) {
					par[param_num]=curr_exchange->get_aa_prop_fac(j, k, SCALING);
					types[param_num]=SCALING_FAC;
					prop_index[param_num]=SCALING;
					prop_rate_index[param_num]=k;
					prop_pns_index[param_num++]=j;
				}
			}
		}
	
	return(param_num);
}



//MultiMatrix model functions

void MultiMatrix_model::describe_results()
{
	int i;

	describe_nuc_results();
	cout<<"Coefficients for AA substitution rate matrices: ";
	for(i=0; i<curr_exchange->get_num_matrices(); i++)
	cout<<"Matrix "<<i<<": "<<curr_exchange->get_matrix_coeff(i,0)<<endl;
}


int MultiMatrix_model::rate_param_size()
{
  return(curr_exchange->get_num_matrices());
}


//Protected MultiMatrix_model
BOOL MultiMatrix_model::change_rate(int rate_num, double *rate_info)
{
	int i;
	BOOL retval=FALSE;
  
	if (fabs(curr_exchange->get_p_non_syn(0, rate_num)-rate_info[0])>FLOAT_TOL)
    {
		for(i=0; i<curr_exchange->get_num_matrices(); i++)
			curr_exchange->set_matrix_coeff(i, 0, rate_info[i]);
      retval=TRUE;
    }
  
  return(retval);
}



double MultiMatrix_model::prob_nonsyn(int start[3], int end[3], Branch *taxa, int rate_num)
{
	int i;
	double net_score=1, mat_score;
	//Uses the stored matrix "the_matrix", which is assumed to be an log-odds type matrix to determine substitution 
	//rates.  ***NOTE: GGCC SHOULD CHECK--6/7/04***

	for(i=0; i<curr_exchange->get_num_matrices(); i++) {
		mat_score=exp(curr_exchange->get_matrix_coeff(i, taxa->get_p_nonsyn_num())*
			(*the_matrices[i])[readchar_to_aa(curr_code->return_AA(start))][readchar_to_aa(curr_code->return_AA(end))]); 
			mat_score=integer_power(mat_score,	curr_code->is_non_synon(start, end)); 
			net_score *=mat_score;
	}
	
	return(net_score); 
}


void MultiMatrix_model::set_num_pns_used()
{
	int i, j;

	num_pns_coeff=0;
	
	for(i=0; i<curr_exchange->get_num_p_non_syn(); i++)
	{
		for(j=0; j<curr_exchange->get_num_matrices(); j++)
			if (curr_exchange->is_matrix_coeff_used(i,j) == TRUE)
				num_pns_coeff++;
	}

}
 

int MultiMatrix_model::set_pns_id_indices(double par[], PARAM_TYPE types[])
{

	int i, j, param_num=0;
	pns_start=0;

	prop_pns_index = new int [num_pns_coeff];
	matrix_index=new int [num_pns_coeff];
	prop_rate_index = new int [num_pns_coeff];

	for(i=0; i<curr_exchange->get_num_p_non_syn(); i++) {
		for(j=0; j<curr_exchange->get_num_matrices(); j++)
			if (curr_exchange->is_matrix_coeff_used(i,j)) {
				par[param_num]=curr_exchange->get_matrix_coeff(j,i);
				types[param_num]=MATRIX_COEFF;
				matrix_index[param_num]=j;
				prop_pns_index[param_num++]=i;
			}
	}
	
	
	return(param_num);	

}




//Public MG_94_JC_model
MG_94_JC_model::MG_94_JC_model ()
{
}


MG_94_JC_model::MG_94_JC_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree, Genetic_code *ccode) 
{
  curr_code=ccode;
  curr_lin_alg=new  Linear_Algebra(ccode);
  assemble (cexchange, cdata, ctree);
}



void MG_94_JC_model::num_params_model()
  //Another Polymorphic function that determines the number of parameters to be
  //minimized based on the model of evolution selected
{
  int mol_clock_params=0;
  
  if (curr_exchange->is_mol_clock_3()==TRUE)
    mol_clock_params=2;

  set_num_pns_used();
  curr_exchange->set_num_params(num_pns_coeff+mol_clock_params+
	  num_rate_prob_params());
  
} //End HKY_model::num_params_model




void MG_94_JC_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
  int i, mc_start, param_num;

  param_num=set_pns_id_indices(par, types);

  if (curr_exchange->are_optimizing_rate_props() == TRUE)
	  param_num=initialize_rate_prob_params(param_num, par, types);
   
  if (curr_exchange->is_mol_clock_3()==TRUE)
    {
      mc_start=param_num;
      switch (curr_exchange->get_clock_type() ) {
      case KS_CLOCK: 
		par[mc_start]=save_k_common;
		types[mc_start]=THREE_TREE_COMMON;
	
		par[mc_start+1]=save_k_dupl;
		types[mc_start+1]=THREE_TREE_SPLIT;
		break;
      case KA_CLOCK:
		par[mc_start]=save_k_common;
		types[mc_start]=THREE_TREE_COMMON_KA;
	
		par[mc_start+1]=save_k_dupl;
		types[mc_start+1]=THREE_TREE_SPLIT_KA;
      }
    }
} //End initialize_parameters



//Public MG_94_K2P_model
MG_94_K2P_model::MG_94_K2P_model () 
{
}


MG_94_K2P_model::MG_94_K2P_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree, Genetic_code *ccode) 
{
  curr_code=ccode;
  curr_lin_alg=new  Linear_Algebra(ccode);
  assemble (cexchange, cdata, ctree);
}



void MG_94_K2P_model::num_params_model()
  //Another Polymorphic function that determines the number of parameters to be
  //minimized based on the model of evolution selected
{ 
  int mol_clock_params=0;
  
  if (curr_exchange->is_mol_clock_3()==TRUE)
    mol_clock_params=2;
  
  set_num_pns_used();
  
  curr_exchange->set_num_params(1+num_pns_coeff
	  +mol_clock_params+num_rate_prob_params());
  
} //End HKY_model::num_params_model





void MG_94_K2P_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
  int i, mc_start, param_num;

  param_num=set_pns_id_indices(par, types);
  
  if (curr_exchange->are_optimizing_rate_props() == TRUE)
	  param_num=initialize_rate_prob_params(param_num, par, types);
  

  par[param_num]=get_trs_trv();
  types[param_num++]=TRS_TRV;

  

  if (curr_exchange->is_mol_clock_3()==TRUE)
    {
      mc_start=param_num;
      switch (curr_exchange->get_clock_type() ) {
      case KS_CLOCK: 
	par[mc_start]=save_k_common;
	types[mc_start]=THREE_TREE_COMMON;
	
	par[mc_start+1]=save_k_dupl;
	types[mc_start+1]=THREE_TREE_SPLIT;
	break;
      case KA_CLOCK:
	par[mc_start]=save_k_common;
	types[mc_start]=THREE_TREE_COMMON_KA;
	
	par[mc_start+1]=save_k_dupl;
	types[mc_start+1]=THREE_TREE_SPLIT_KA;
      }


    }

} //End initialize_parameters




//Public MG_94_HKY_model
MG_94_HKY_model::MG_94_HKY_model () 
{
}


MG_94_HKY_model::MG_94_HKY_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree, Genetic_code *ccode) 
{
  curr_code=ccode;
  curr_lin_alg=new  Linear_Algebra(ccode);
  assemble (cexchange, cdata, ctree);
}


void MG_94_HKY_model::num_params_model()
  //Another Polymorphic function that determines the number of parameters to be
  //minimized based on the model of evolution selected
{
  int mol_clock_params=0;
  
  if (curr_exchange->is_mol_clock_3()==TRUE)
    mol_clock_params=2;
  set_num_pns_used();

  if (curr_exchange->fixed_basefreq()==FALSE)
    curr_exchange->set_num_params(4+num_pns_coeff
	+mol_clock_params+num_rate_prob_params());
  else 
    curr_exchange->set_num_params(1+num_pns_coeff
	+mol_clock_params+num_rate_prob_params());
} //End HKY_model::num_params_model



void MG_94_HKY_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
  int i, bfend=0, mc_start, param_num;


  param_num=set_pns_id_indices(par, types);
  if (curr_exchange->are_optimizing_rate_props() == TRUE)
	  param_num=initialize_rate_prob_params(param_num, par, types);
  
 

  par[param_num]=get_trs_trv();
  types[param_num++]=TRS_TRV;
  

  bfend=param_num;
  if (curr_exchange->fixed_basefreq()==FALSE) {
    par[param_num]=pur_pyr_split;
    types[param_num++]=PUR_PYR_SPLIT;
    
    par[param_num]=a_g_split;
    types[param_num++]=A_G_SPLIT;
    
    par[param_num]=c_t_split;
    types[param_num++]=C_T_SPLIT;  
    bfend=4;
  }
 
  
  if (curr_exchange->is_mol_clock_3()==TRUE)
    {
      mc_start=param_num;
      switch (curr_exchange->get_clock_type() ) {
      case KS_CLOCK: 
		par[mc_start]=save_k_common;
		types[mc_start]=THREE_TREE_COMMON;
	
		par[mc_start+1]=save_k_dupl;
		types[mc_start+1]=THREE_TREE_SPLIT;
		break;
      case KA_CLOCK:
		par[mc_start]=save_k_common;
		types[mc_start]=THREE_TREE_COMMON_KA;
	
		par[mc_start+1]=save_k_dupl;
		types[mc_start+1]=THREE_TREE_SPLIT_KA;
      }


    }
} //End initialize_parameters





//Public MG_94_GG_98_model
 MG_94_GG_98_model::MG_94_GG_98_model () 
{
}


MG_94_GG_98_model::MG_94_GG_98_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree, Genetic_code *ccode) 
{
  curr_code=ccode;
  curr_lin_alg=new  Linear_Algebra(ccode);
  assemble (cexchange, cdata, ctree);
}



void MG_94_GG_98_model::num_params_model()
  //Another Polymorphic function that determines the number of parameters to be
  //minimized based on the model of evolution selected
{
  cerr<<"Powell optimization not allowed with Galtier and Gaut 1998 Model!\n";
} //End HKY_model::num_params_model



void MG_94_GG_98_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
 cerr<<"Powell optimization not allowed with Galtier and Gaut 1998 Model!\n";
} //End initialize_parameters




//Public C_00_JC_model
 C_00_JC_model::C_00_JC_model () 
{
}


C_00_JC_model::C_00_JC_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree, Genetic_code *ccode) 
{
  curr_code=ccode;
  curr_lin_alg=new  Linear_Algebra(ccode);
  assemble(cexchange, cdata, ctree);
}



void C_00_JC_model::num_params_model()
  //Another Polymorphic function that determines the number of parameters to be
  //minimized based on the model of evolution selected
{

	set_num_pns_used();
	curr_exchange->set_num_params(num_pns_coeff+num_pi_coeff+
		+num_rate_prob_params());

} //End HKY_model::num_params_model




void C_00_JC_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
  int param_num;

  param_num=set_pns_id_indices(par, types);
  
  if (curr_exchange->are_optimizing_rate_props() == TRUE)
	  param_num=initialize_rate_prob_params(param_num, par, types);
  
  
} //End initialize_parameters



//Public C_00_K2P_model
 C_00_K2P_model::C_00_K2P_model () 
{
}


C_00_K2P_model::C_00_K2P_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree, Genetic_code *ccode) 
{
  curr_code=ccode;
  curr_lin_alg=new  Linear_Algebra(ccode);
  assemble (cexchange, cdata, ctree);
}

void C_00_K2P_model::num_params_model()
  //Another Polymorphic function that determines the number of parameters to be
  //minimized based on the model of evolution selected
{

	set_num_pns_used();
	curr_exchange->set_num_params(1+num_pns_coeff+num_pi_coeff+
		num_rate_prob_params());

} //End HKY_model::num_params_model




void C_00_K2P_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
  int i, param_num=1;
  
  param_num=set_pns_id_indices(par, types);

  if (curr_exchange->are_optimizing_rate_props() == TRUE)
	  param_num=initialize_rate_prob_params(param_num, par, types);
  
 
  par[param_num]=get_trs_trv();
  types[param_num++]=TRS_TRV;

 
} //End initialize_parameters



//Public C_00_HKY_model
C_00_HKY_model::C_00_HKY_model () 
{
}


C_00_HKY_model::C_00_HKY_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree, Genetic_code *ccode) 
{
  curr_code=ccode;
  curr_lin_alg=new  Linear_Algebra(ccode);
  assemble (cexchange, cdata, ctree);
}


void C_00_HKY_model::num_params_model()
  //Another Polymorphic function that determines the number of parameters to be
  //minimized based on the model of evolution selected
{

  set_num_pns_used();

  int basefreq_params;
  if (curr_exchange->fixed_basefreq()==TRUE)
    basefreq_params=0;
  else
    basefreq_params=3;

  curr_exchange->set_num_params(1+basefreq_params+num_pns_coeff+num_pi_coeff+
	  num_rate_prob_params());
} //End HKY_model::num_params_model




void C_00_HKY_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
  int param_num;
  
  param_num=set_pns_id_indices(par, types);

  if (curr_exchange->are_optimizing_rate_props() == TRUE)
	  param_num=initialize_rate_prob_params(param_num, par, types);
  
 

  par[param_num]=get_trs_trv();
  types[param_num++]=TRS_TRV;
 
  if(curr_exchange->fixed_basefreq()==FALSE) {
    par[param_num]=pur_pyr_split;
    types[param_num++]=PUR_PYR_SPLIT;
    
    par[param_num]=a_g_split;
    types[param_num++]=A_G_SPLIT;
    
    par[param_num]=c_t_split;
    types[param_num++]=C_T_SPLIT;  
  }

   
} //End initialize_parameters


//Public C_00_GG_98_model
C_00_GG_98_model::C_00_GG_98_model () 
{
}


C_00_GG_98_model::C_00_GG_98_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree, Genetic_code *ccode) 
{
  curr_code=ccode;
  curr_lin_alg=new  Linear_Algebra(ccode);
  assemble (cexchange, cdata, ctree);
}


void C_00_GG_98_model::num_params_model()
  //Another Polymorphic function that determines the number of parameters to be
  //minimized based on the model of evolution selected
{
 cerr<<"Powell optimization not allowed with Galtier and Gaut 1998 Model!\n";
} //End HKY_model::num_params_model




void C_00_GG_98_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
  cerr<<"Powell optimization not allowed with Galtier and Gaut 1998 Model!\n";
}



//Public LCAP_JC_model
 LCAP_JC_model::LCAP_JC_model () 
{
}


LCAP_JC_model::LCAP_JC_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree, Genetic_code *ccode) 
{
  curr_code=ccode;
  curr_lin_alg=new  Linear_Algebra(ccode);
  assemble (cexchange, cdata, ctree);
}



void LCAP_JC_model::num_params_model()
  //Another Polymorphic function that determines the number of parameters to be
  //minimized based on the model of evolution selected
{
	set_num_aa_props_used();

	curr_exchange->set_num_params(num_aa_prop_coeff+1+num_rate_prob_params());
} //End LCAP_JC_model::num_params_model




void LCAP_JC_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
	int i, param_num;

	param_num=set_aa_prop_id_indices(par, types);

	if (curr_exchange->are_optimizing_rate_props() == TRUE)
	  param_num=initialize_rate_prob_params(param_num, par, types);
  
 
} //End initialize_parameters



//Public LCAP_K2P_model
 LCAP_K2P_model::LCAP_K2P_model () 
{
}


LCAP_K2P_model::LCAP_K2P_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree, Genetic_code *ccode) 
{
  curr_code=ccode;
  curr_lin_alg=new  Linear_Algebra(ccode);
  assemble (cexchange, cdata, ctree);
}



void LCAP_K2P_model::num_params_model()
  //Another Polymorphic function that determines the number of parameters to be
  //minimized based on the model of evolution selected
{
	set_num_aa_props_used();

	curr_exchange->set_num_params(1+num_aa_prop_coeff+
		num_rate_prob_params());
} //End HKY_model::num_params_model




void LCAP_K2P_model::intialize_parameters (double par[], PARAM_TYPE types[])
{

	int i, param_num=0;

	param_num=set_aa_prop_id_indices(par, types);

	if (curr_exchange->are_optimizing_rate_props() == TRUE)
	  param_num=initialize_rate_prob_params(param_num, par, types);
  
 
  
	par[param_num]=get_trs_trv();
	types[param_num]=TRS_TRV;

} //End initialize_parameters



//Public LCAP_HKY_model
 LCAP_HKY_model::LCAP_HKY_model () 
{
}


LCAP_HKY_model::LCAP_HKY_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree, Genetic_code *ccode)
{
  curr_code=ccode;
  curr_lin_alg=new  Linear_Algebra(ccode);
  assemble (cexchange, cdata, ctree);
}



void LCAP_HKY_model::num_params_model()
  //Another Polymorphic function that determines the number of parameters to be
  //minimized based on the model of evolution selected
{
  set_num_aa_props_used();

  if (curr_exchange->fixed_basefreq()==FALSE)
    curr_exchange->set_num_params(4+num_aa_prop_coeff+
		num_rate_prob_params());
  else
     curr_exchange->set_num_params(1+num_aa_prop_coeff+
		num_rate_prob_params());
} //End HKY_model::num_params_model




void LCAP_HKY_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
	int i, param_num=0;

	param_num=set_aa_prop_id_indices(par, types);

	if (curr_exchange->are_optimizing_rate_props() == TRUE)
	  param_num=initialize_rate_prob_params(param_num, par, types);
  
 

	par[param_num]=get_trs_trv();
	types[param_num++]=TRS_TRV;
	if (curr_exchange->fixed_basefreq()==FALSE) {
		par[param_num]=pur_pyr_split;
		types[param_num]=PUR_PYR_SPLIT;
    
		par[param_num]=a_g_split;
		types[param_num++]=A_G_SPLIT;
    
		par[param_num]=c_t_split;
		types[param_num++]=C_T_SPLIT;  
    
	}

	for(i=0; i<curr_exchange->get_num_params(); i++)
		cout<<i<<": "<<par[i]<<endl;

} //End initialize_parameters

#if 0
double LCAP_HKY_model::find_ut(Branch *taxa)
{
	int i, j, k, l, start_codon, end_codon, start[3], end[3];
	double val, val2, freq, diag_sum, diag_sum2, so_far, so_far2, stop_freq;
	
	stop_freq=0;
  
	for (start_codon=0; start_codon<curr_exchange->get_condlike_size(); start_codon++)
    { 
      if (curr_code->is_stop(start_codon)==TRUE) {
		for (j=0; j<3; j++)
			start[j]=curr_code->get_pos_n(start_codon, j);
	    
		stop_freq+=get_basefreq(start[0], 0, taxa)*
		get_basefreq(start[1], 1, taxa)*get_basefreq(start[2], 2, taxa);
      }   
    }
	stop_excess=1/(1-stop_freq);
  

	for(i=0; i<curr_exchange->get_num_branches(); i++)
    {      
		prop_syn=prop_nsyn=neu_prop_syn=neu_prop_nsyn=0;
		diag_sum=diag_sum2=0;
		for (start_codon=0; start_codon<curr_exchange->get_condlike_size(); start_codon++)
		{
			if (curr_code->is_stop(start_codon)!=TRUE)
			{ 
				so_far=so_far2=0;
				for (j=0; j<3; j++)
					start[j]=curr_code->get_pos_n(start_codon, j);
	      
				freq= get_basefreq(start[0], 0, taxa)*
				get_basefreq(start[1], 1, taxa)*
				get_basefreq(start[2], 2, taxa)*stop_excess;
	      
				for(j=0; j<3; j++)
				{
					for (k=0; k<3; k++)
					{
						end_codon=curr_code->one_diff(start_codon, j, k);
		      
						for (l=0; l<3; l++)
							end[l]=curr_code->get_pos_n(end_codon, l);
		      
						if (curr_code->is_stop(end_codon)!=TRUE) {
							val=q_matrix_entry(start[0], end[0], 0, taxa)* 
							  q_matrix_entry(start[1], end[1], 1, taxa)*
							  q_matrix_entry(start[2], end[2], 2, taxa)*
							  prob_nonsyn(start, end, taxa, 0)*stop_excess;
							val2=	q_matrix_entry(start[0], end[0], 0, taxa)* 
							  q_matrix_entry(start[1], end[1], 1, taxa)*
							  q_matrix_entry(start[2], end[2], 2, taxa)*stop_excess;
			
						if (curr_code->is_non_synon(start, end)==FALSE)
						{ 
							prop_syn+=freq*val;
							neu_prop_syn+=freq*val2;
						}
			
						so_far +=val;
						so_far2+=val2;
					}
				}
			}
			diag_sum+=freq*so_far;
			diag_sum2+=freq*so_far2;
	    }    
	}
      
      
    prop_syn/=diag_sum;
    neu_prop_syn/=diag_sum2;
    prop_nsyn=1.0-prop_syn;
    neu_prop_nsyn=1.0-neu_prop_syn;
      
    cout<<"Expected: "<<taxa->expect_subs_site()<<"  New len: "<<((3.0*neu_prop_syn)/prop_syn)*(taxa->expect_subs_site())<<endl;
    if (((3.0*neu_prop_syn)/prop_syn)*(taxa->expect_subs_site()) < MIN_BRLEN) 
		return(((3.0*neu_prop_syn)/prop_syn)*(taxa->expect_subs_site()));
	
     else 
		 return(MIN_BRLEN);
    }
  

}
#endif


//Public LCAP_GG_98_model
 LCAP_GG_98_model::LCAP_GG_98_model () 
{
}


LCAP_GG_98_model::LCAP_GG_98_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree, Genetic_code *ccode) 
{
  curr_code=ccode;
  curr_lin_alg=new  Linear_Algebra(ccode);
  assemble (cexchange, cdata, ctree);
}


void LCAP_GG_98_model::num_params_model()
  //Another Polymorphic function that determines the number of parameters to be
  //minimized based on the model of evolution selected
{
 cerr<<"Powell optimization not allowed with Galtier and Gaut 1998 Model!\n";
} 




void LCAP_GG_98_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
  cerr<<"Powell optimization not allowed with Galtier and Gaut 1998 Model!\n";
} //End initialize_parameters





//Public Matrix_Rates_JC_model
MultiMatrix_JC_model::MultiMatrix_JC_model ()
{
}


MultiMatrix_JC_model::MultiMatrix_JC_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree, Genetic_code *ccode, AA_matrix **the_mats) 
{
  curr_code=ccode;
  curr_lin_alg=new  Linear_Algebra(ccode);
  the_matrices=the_mats;
  assemble (cexchange, cdata, ctree);
}



void MultiMatrix_JC_model::num_params_model()
  //Another Polymorphic function that determines the number of parameters to be
  //minimized based on the model of evolution selected
{
	set_num_pns_used();
	if (curr_exchange->get_num_rates() > 1)
		cerr<<"Error in matrix models: Only one rate will be used\n";
	curr_exchange->set_num_params(num_pns_coeff);

} //End HKY_model::num_params_model




void MultiMatrix_JC_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
	int i, j, param_num=0;

	matcoeff_start=0;
	param_num=set_pns_id_indices(par, types);
 
} //End initialize_parameters


//Public MultiMatrix_K2P_model
MultiMatrix_K2P_model::MultiMatrix_K2P_model () 
{
}


MultiMatrix_K2P_model::MultiMatrix_K2P_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree, Genetic_code *ccode, AA_matrix **the_mats) 
{
  curr_code=ccode;
  curr_lin_alg=new  Linear_Algebra(ccode);
  the_matrices=the_mats;
  assemble (cexchange, cdata, ctree);
}



void MultiMatrix_K2P_model::num_params_model()
  //Another Polymorphic function that determines the number of parameters to be
  //minimized based on the model of evolution selected
{
	set_num_pns_used();
	if (curr_exchange->get_num_rates() > 1)
		cerr<<"Error in matrix models: Only one rate will be used\n";
    curr_exchange->set_num_params(1+num_pns_coeff);
  
} //End K2P_model::num_params_model





void MultiMatrix_K2P_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
	int i, j, param_num=0;

	param_num=set_pns_id_indices(par, types);
	matcoeff_start=0;	


	par[param_num]=get_trs_trv();
	types[param_num++]=TRS_TRV;

    
} //End initialize_parameters




//Public Matrix_Rates_HKY_model
MultiMatrix_HKY_model::MultiMatrix_HKY_model () 
{
}


MultiMatrix_HKY_model::MultiMatrix_HKY_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree, Genetic_code *ccode, AA_matrix **the_mats) 
{
  curr_code=ccode;
  curr_lin_alg=new  Linear_Algebra(ccode);
  the_matrices=the_mats;
  assemble (cexchange, cdata, ctree);
}


void MultiMatrix_HKY_model::num_params_model()
  //Another Polymorphic function that determines the number of parameters to be
  //minimized based on the model of evolution selected
{
	set_num_pns_used();
	if (curr_exchange->get_num_rates() > 1)
		cerr<<"Error in matrix models: Only one rate will be used\n";
	if (curr_exchange->fixed_basefreq()==TRUE)
		curr_exchange->set_num_params(1+num_pns_coeff);
	else
		curr_exchange->set_num_params(4+num_pns_coeff);
} //End HKY_model::num_params_model



void MultiMatrix_HKY_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
	int i, j, param_num;


	param_num=set_pns_id_indices(par, types);
	matcoeff_start=0;	


	par[param_num]=get_trs_trv();
	types[param_num++]=TRS_TRV;

	if (curr_exchange->fixed_basefreq()==FALSE) {
		par[param_num]=pur_pyr_split;
		types[param_num++]=PUR_PYR_SPLIT;
    
		par[param_num]=a_g_split;
		types[param_num++]=A_G_SPLIT;
    
		par[param_num]=c_t_split;
		types[param_num++]=C_T_SPLIT; 
	}
 
 
 

} //End initialize_parameters




