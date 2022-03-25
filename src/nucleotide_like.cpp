//Copyright 1999-2005 Gavin Conant

#include <iostream>
#include <math.h>
#include <iomanip>
#include "nucleotide_like.h"
#include "powell.h"

using namespace::std;

#define A  (1+1/PurPyrfreqs(taxa, site2)*(curr_exchange->get_trs_trv()-1))
 

//Public Nucleotide_model
void Nucleotide_model::describe_results()
{
  describe_nuc_results();
}

double Nucleotide_model::root_freq(int site)
{
  return(get_basefreq(site, curr_tree->find_root()));
} //End Nucleotide_model::root_freq


int Nucleotide_model::rate_param_size()
{
  return(1);
}


//Protected Nucleotide_model
void Nucleotide_model::initialize_arrays()
{
}





void Nucleotide_model::calc_transprobs(Branch *taxa, int rate_num)
{
  int i,j;
  double PurPyrfreqs[4], InvPurPyrfreqs[4];
  
  PurPyrfreqs[0]=PurPyrfreqs[2]=get_basefreq(0, taxa)+get_basefreq(2, taxa);
  PurPyrfreqs[1]=PurPyrfreqs[3]=get_basefreq(1, taxa)+get_basefreq(3, taxa);
  
  for (i=0; i<4; i++)
    InvPurPyrfreqs[i]=1/PurPyrfreqs[i];

  for (i=0;i<4; i++)
      for (j=0; j<4;j++)
	{
	  if (taxa->get_brnlen()>=MIN_BRLEN)
	    {
	      if (i==j)
		taxa->set_trpb(rate_num, i, j, get_basefreq(j, taxa)+
				    get_basefreq(j, taxa)*(InvPurPyrfreqs[j]-1)*
				    exp(-curr_exchange->get_rate(rate_num)*taxa->get_brnlen())+
				    (InvPurPyrfreqs[j]*(PurPyrfreqs[j]-get_basefreq(j, taxa)))*
				    exp(-curr_exchange->get_rate(rate_num)*taxa->get_brnlen()*
					(1+PurPyrfreqs[j]*(get_trs_trv()-1))));
	      else if((i==j+2) || (i==j-2))
		taxa->set_trpb(rate_num, i, j, get_basefreq(j, taxa)+
				    get_basefreq(j, taxa)*(InvPurPyrfreqs[j]-1)*
				    exp(-curr_exchange->get_rate(rate_num)*taxa->get_brnlen())-
				    (InvPurPyrfreqs[j]*get_basefreq(j, taxa))*
				    exp(-curr_exchange->get_rate(rate_num)*taxa->get_brnlen()*
					(1+PurPyrfreqs[j]*(get_trs_trv()-1))));
	      else 
		taxa->set_trpb(rate_num, i, j, get_basefreq(j, taxa)*
				    (1-exp(-curr_exchange->get_rate(rate_num)*taxa->get_brnlen())));
	    }
	  
	  else
	    {
	      if (taxa->get_brnlen()!=0.0)
		taxa->set_brnlen(0.0); 
	      if(i==j)
		taxa->set_trpb(rate_num, i, j, 1.0);
	      else
		taxa->set_trpb(rate_num, i, j, 0.0);	    
	    }
	}
  
  //Sets Transprobs for gaps to 1
  for (i=0; i<5; i++)
    {
      taxa->set_trpb(rate_num, i, 4, 1.0);
      taxa->set_trpb(rate_num, 4, i, 1.0);
    }

}



BOOL Nucleotide_model::change_rate(int rate_num, double *rate_info)
{
  BOOL retval=FALSE;
  if(fabs(curr_exchange->get_rate(rate_num)-rate_info[0])>FLOAT_TOL)
    {
      curr_exchange->set_rate(rate_num, rate_info[0]);
      retval=TRUE;
    }
  return(retval);
}



double Nucleotide_model::prob_nonsyn(int start[3], int end[3], Branch *taxa, int rate_num)
{
  cerr<<"Invalid call to codon function for nucleotide likelihood model\n";
  return(0);
}

void JC_Nucleotide_model::calc_first_derivatives(Branch *taxa, double ut)
{
	calc_first_derivatives(taxa, ut, 0);
}


void JC_Nucleotide_model::calc_second_derivatives(Branch *taxa, double ut)
{
	calc_second_derivatives(taxa,ut, 0);
}

void JC_Nucleotide_model::calc_first_derivatives(Branch *taxa, double ut, int rate)
{
  int site1, site2;
  for (site1=0; site1<curr_exchange->get_condlike_size(); site1++)
    for (site2=0; site2<curr_exchange->get_condlike_size(); site2++) {
      if (site1 == site2)
	taxa->set_trpb_prime(site1, site2, rate, -0.75*exp(-ut*curr_exchange->get_rate(rate)));
      else
	taxa->set_trpb_prime(site1, site2, rate, 0.25*exp(-ut*curr_exchange->get_rate(rate)));
    }
}

void JC_Nucleotide_model::calc_second_derivatives(Branch *taxa, double ut, int rate)
{ 
  int site1, site2;
  
  for (site1=0; site1<curr_exchange->get_condlike_size(); site1++)
    for (site2=0; site2<curr_exchange->get_condlike_size(); site2++) {
      if (site1 == site2)
	taxa->set_trpb_double_prime(site1, site2, rate, 0.75*exp(-ut*curr_exchange->get_rate(rate)));
      else
	taxa->set_trpb_double_prime(site1, site2, rate, -0.25*exp(-ut*curr_exchange->get_rate(rate)));
    }
}



void JC_Nucleotide_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
 
} //End initialize_parameters


double JC_Nucleotide_model::get_obs_trs_trv_ratio()
{
  return(0.50);
}




void K2P_Nucleotide_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
  int i;
 
  par[0]=get_trs_trv();
  types[0]=TRS_TRV;

} //End initialize_parameters

void K2P_Nucleotide_model::calc_first_derivatives(Branch *taxa, double ut)
{
	calc_first_derivatives(taxa, ut, 0);
}


void K2P_Nucleotide_model::calc_second_derivatives(Branch *taxa, double ut)
{
	calc_second_derivatives(taxa,ut, 0);
}


void K2P_Nucleotide_model::calc_first_derivatives(Branch *taxa, double ut, int rate)
{
  int site1, site2;
  double kappa_p_1;
 
  kappa_p_1=get_trs_trv()+1;
  ut*=curr_exchange->get_rate(rate);

  for (site1=0; site1<curr_exchange->get_condlike_size(); site1++)
    for (site2=0; site2<curr_exchange->get_condlike_size(); site2++) {

      if (site1 == site2)
	taxa->set_trpb_prime(site1, site2, rate, -0.25*exp(-ut)-((kappa_p_1/4)*
							   exp(-ut*(kappa_p_1/2))));
      else if (is_transition(site1, site2)==TRUE)
	taxa->set_trpb_prime(site1, site2, rate, -0.25*exp(-ut)+((kappa_p_1/4)
							   *exp(-ut*(kappa_p_1/2))));
      else
	taxa->set_trpb_prime(site1, site2, rate, 0.25*exp(-ut)); 
    }
}



void K2P_Nucleotide_model::calc_second_derivatives(Branch *taxa, double ut, int rate)
{
 int site1, site2;
 double kappa_p_1;
 
  kappa_p_1=get_trs_trv()+1;
  ut*=curr_exchange->get_rate(rate);

 
  for (site1=0; site1<curr_exchange->get_condlike_size(); site1++)
   for (site2=0; site2<curr_exchange->get_condlike_size(); site2++) {
     if (site1 == site2)
       taxa->set_trpb_double_prime(site1, site2, rate, 0.25*exp(-ut)+
				   (((kappa_p_1*kappa_p_1)/8)*
				    exp(-ut*(kappa_p_1/2))));
     else if (is_transition(site1, site2)==TRUE)
       taxa->set_trpb_double_prime(site1, site2, rate, 0.25*exp(-ut)-(((kappa_p_1*kappa_p_1)/8)*
								exp(-ut*(kappa_p_1/2))));
     else
       taxa->set_trpb_double_prime(site1, site2, rate, -0.25*exp(-ut));
   }
}


double K2P_Nucleotide_model::get_obs_trs_trv_ratio()
{
  return(get_trs_trv()/2.0);
}



void HKY_Nucleotide_model::num_params_model() 
{
  if (curr_exchange->fixed_basefreq()==FALSE) 
    curr_exchange->set_num_params(4);
  else
    curr_exchange->set_num_params(1);
}

void HKY_Nucleotide_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
  int i;
  double tbfs[4];

  if (curr_exchange->fixed_basefreq()==FALSE) {
    bf_end=3;


    par[0]=pur_pyr_split;
    types[0]=PUR_PYR_SPLIT;
    
    par[1]=a_g_split;
    types[1]=A_G_SPLIT;
    
    par[2]=c_t_split;
    types[2]=C_T_SPLIT;
  } 
  else
    bf_end=0;

  par[bf_end]=get_trs_trv();
  types[bf_end]=TRS_TRV;
  
} //End initialize_parameters


void HKY_Nucleotide_model::calc_first_derivatives(Branch *taxa, double ut)
{
	calc_first_derivatives(taxa, ut, 0);
}


void HKY_Nucleotide_model::calc_second_derivatives(Branch *taxa, double ut)
{
	calc_second_derivatives(taxa,ut, 0);
}


void HKY_Nucleotide_model::calc_first_derivatives(Branch *taxa, double ut, int rate)
{
  int site1, site2; 
  
  ut*=curr_exchange->get_rate(rate);

  for (site1=0; site1<curr_exchange->get_condlike_size(); site1++)
    for (site2=0; site2<curr_exchange->get_condlike_size(); site2++) {
   
     
      if (site1 == site2)
	taxa->set_trpb_prime(site1, site2, rate, -get_basefreq(site2, taxa)*(1/PurPyrfreqs(taxa, site2)-1)*
			     exp(-ut)-A*((PurPyrfreqs(taxa, site2)-get_basefreq(site2, taxa))/PurPyrfreqs(taxa, site2))*
			     exp(-ut*A));
      else if (is_transition(site1, site2)==TRUE)
	taxa->set_trpb_prime(site1, site2, rate, -get_basefreq(site2, taxa)*(1/PurPyrfreqs(taxa, site2)-1)*
			     exp(-ut)+A*(get_basefreq(site2, taxa)/PurPyrfreqs(taxa, site2))*exp(-ut*A));
      else
	taxa->set_trpb_prime(site1, site2, rate, get_basefreq(site2, taxa)*exp(-ut)); 
    }
}


void HKY_Nucleotide_model::calc_second_derivatives(Branch *taxa, double ut, int rate)
{  
  int site1, site2; 
  ut*=curr_exchange->get_rate(rate);

  
 
  for (site1=0; site1<curr_exchange->get_condlike_size(); site1++)
    for (site2=0; site2<curr_exchange->get_condlike_size(); site2++) {
     
      if (site1 == site2)
	taxa->set_trpb_double_prime(site1, site2, rate, get_basefreq(site2, taxa)*(1/PurPyrfreqs(taxa, site2)-1)
				    *exp(-ut)+A*A*((PurPyrfreqs(taxa, site2)-
						    get_basefreq(site2, taxa))/PurPyrfreqs(taxa, site2))*exp(-ut*A));
      else if (is_transition(site1, site2)==TRUE)
	taxa->set_trpb_double_prime(site1, site2, rate, get_basefreq(site2, taxa)*(1/PurPyrfreqs(taxa, site2)-1)*
				    exp(-ut)-A*A*(get_basefreq(site2, taxa)/PurPyrfreqs(taxa, site2))*exp(-ut*A));
      else
	taxa->set_trpb_double_prime(site1, site2, rate, -get_basefreq(site2, taxa)*exp(-ut));
    } 
}


double HKY_Nucleotide_model::get_obs_trs_trv_ratio()
{
    return(get_trs_trv()*((get_basefreq(0,(*curr_tree)[0])*get_basefreq(2,(*curr_tree)[0])+
			   get_basefreq(1,(*curr_tree)[0])*get_basefreq(3,(*curr_tree)[0]))/
			  ((get_basefreq(0,(*curr_tree)[0])+get_basefreq(2,(*curr_tree)[0]))*
			   (get_basefreq(1,(*curr_tree)[0])+get_basefreq(3,(*curr_tree)[0])))));
  
}




double HKY_Nucleotide_model::PurPyrfreqs(Branch *taxa, int base)
{
  return(get_basefreq(base, taxa)+ get_basefreq(((base+2)%4), taxa));
}




void GG_98_Nucleotide_model::num_params_model()
  //Another Polymorphic function that determines the number of parameters to be
  //minimized based on the model of evolution selected
{
  cerr<<"Powell optimization not allowed with Galtier and Gaut 1998 Model!\n";
} //End HKY_model::num_params_model




void GG_98_Nucleotide_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
  cerr<<"Powell optimization not allowed with Galtier and Gaut 1998 Model!\n";
}




double HKY_Nucleotide_12_3_model::find_lnL_on_tree_pos_diff(Tree *ctree)
{
  int i, j;
  double temp=0;

  if (ctree!=0)
    {
      curr_tree=ctree;
      
      for (i=0; i<curr_exchange->get_num_branches(); i++)
	{
	  calc_transprobs((*curr_tree)[i], 0);
	  calc_transprobs((*curr_tree)[i], 1);
	}
    }

  for (i=0; i<curr_exchange->get_num_localities(); i++)
      {
	if (i+1%3!=0)
	  temp+=-log(prob_w_rate(i, 0));
	else
	  temp+=-log(prob_w_rate(i, 1));	
      }
  return(temp);
}



void HKY_Nucleotide_12_3_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
 int i;
 double tbfs[4];
 
 bf_end=3;
 
 
 par[0]=pur_pyr_split;
 types[0]=PUR_PYR_SPLIT;
 
 par[1]=a_g_split;
 types[1]=A_G_SPLIT;
 
 par[2]=c_t_split;
 types[2]=C_T_SPLIT;
 
 par[bf_end]=get_trs_trv();
 types[bf_end]=TRS_TRV;
 
 par[bf_end+1]=curr_exchange->get_rate(1);
 types[bf_end+1]=THIRD_POS_RATE;
}
