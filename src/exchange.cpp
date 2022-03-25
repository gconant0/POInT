//argument.cpp
//Copyright 1999-2002 Gavin Conant


#include "exchange.h"

#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>

#ifdef _OPEN_MP_VERSION_
#include "omp.h"
#endif


using namespace::std;
//#define DEBUG

Exchange::Exchange()
  //Constructor for Exchange class
{
  int i, j;
  
  curr_groups=0;
  num_taxa=num_sites=num_branches=0;
  strcpy(treefile, "\0");
  strcpy(datafile, "\0");
  strcpy(gencode_file, "\0");


  nonstandard_gencode=FALSE;
  model_codons=FALSE;
  pre_optimize=FALSE;
  branch_freqs=FALSE;
  have_dataset=TRUE;
  use_full_conprobs=FALSE;			 //Stores all site node conditional probs for non-rev model brn opt
  opt_rate_probs=FALSE;
  allow_LCAP_pos_weights=TRUE;
  rooted_tree=FALSE;
  fix_zero_brns=FALSE;
  track_gaps_as_missing = FALSE;
  scale_likelihood=FALSE;
    guess_order=FALSE;
    use_break_pos=FALSE;
  
  codon_freqs=FALSE;
  basefreqs_fixed=FALSE;
  three_tree_mol_clock=FALSE;
  is_kaks1=FALSE;
    all_taxa_used=TRUE;
    taxa_used=0;
  branch_kaks1=0;
  clock_type=KS_CLOCK;              /*Only relevent for three-tree test*/
  curr_genetic_code=UNIVERSAL;
  condlike_size=4;
  num_localities=&Exchange::num_sites;
  site_rate_func=&Exchange::single_rate_site_rates;
  site_rate_num_func=&Exchange::single_rate_site_rate_num;
  site_rate_prob_func=&Exchange::singe_rate_site_rate_prob;
  site_rates=0; 
  site_rate_probs=0;
  dupl_parallel_rate=dupl_fix_rate=0;
  strand_switch_prob=0;
  fix_rate_scale=0.01;
  loss_rate_scale=0.05;
  fix_loss_rate=0.05;
  snp_gains_to_losses=1.0;
  snp_loss_to_switch=1.0;
  
  num_matrices=0;
  matrix_coeff=0;
  matrix_coeff_used=0;
    WGX_depth=2;
    max_breaks=0;

  //The model is set to HKY by default
  //If the trs/trv ratio and/or basefreqs are not entered by the user
  //or in the treefile, this gives the same results as the JC model
  current_model=HKY_NUCLEOTIDE;

	num_rates=1;
  rate_type=SINGLE_RATE;
  rates=new double [num_rates];
  aa_rate_index= new int *[num_rates];
	for(i=0; i<num_rates ; i++)
		aa_rate_index[i]=new int [get_num_aa_props()+1];
 
  num_p_nonsyn=0;
	p_non_syn=0;
 
  for(i=0; i<get_num_aa_props()+1; i++)
	  aa_property_allowed[i]=TRUE;
  num_live_properties=NUM_AA_PROPS;

  set_num_p_nonsyn(1);

  num_patterns=&Exchange::num_p_nonsyn;

  set_rates();
	for(i=0; i<num_rates; i++) {
		for(j=0; j<get_num_aa_props()+1; j++)
			aa_rate_index[i][j]=i;
	}

	
  num_nonsyn_params=&Exchange::standard_model_nonsyn_params;
	num_unallowed_SNP_states=1;
	unallowed_states[0]=SNP_ABSENT;
  
	

  //Note that this is the ratio of events, so setting it to 0.5 corresponds to
  //no difference between transitions and transversions
  //(the kappa parameter is set in the Like_model class.  For 
  //trs/trv=0.5, kappa=1)
  obs_trs_trv=0.50;
  kappa=1.0;

  for(i=0; i<4; i++)
    Exchange::basefreqs[i]=0.25;

  rate_probs = 0;

//I believe the number of open mp threads is avalaible only in parallel sections
//but we need it for memory allocation in the Tree object.  Save it here.
#ifdef _OPEN_MP_VERSION_
#pragma omp parallel
	{
	num_open_mp_threads=omp_get_num_threads();
	}
	//cout<<"Using "<<num_open_mp_threads<<" threads for this run"<<endl;
#endif
} //End Exchange::Exchange (constructor)


Exchange& Exchange::operator=(Exchange &assign_from)
{
  int i, j, *sites;
  double freqs[3][4];
	
	//cout<<"In exchange copy operator\n";

  strcpy(treefile, assign_from.get_treefile());
  strcpy(datafile, assign_from.get_datafile());
  strcpy(gencode_file, assign_from.get_gencode_file());

  set_num_taxa(assign_from.get_num_taxa());
  set_num_sites(assign_from.get_num_sites()); 

  if (num_rates != assign_from.get_num_rates())
    set_num_rates(assign_from.get_num_rates());
	
	if (assign_from.get_num_p_non_syn() != num_p_nonsyn) set_num_p_nonsyn(assign_from.get_num_p_non_syn());

  rate_type=assign_from.get_site_rate_type();
  if (assign_from.have_site_rate_probs() == TRUE) {
	use_prob_site_rates(num_rates);
	for(i=0; i<num_sites; i++)
		for(j=0; j<num_rates; j++)
			set_site_rate_prob(i, j, assign_from.get_site_rate_prob(i, j));
  }

    num_params=assign_from.get_num_params();
    nonstandard_gencode=assign_from.non_standard_code();
    if (assign_from.get_model() != DUPL_ARBITRARY)
        set_model(assign_from.get_model());
    else
        set_model(assign_from.get_model(), assign_from.get_condlike_size());
    set_dataformat(assign_from.get_dataformat());
    set_dataset(assign_from.get_dataset());
    set_have_data(assign_from.have_data());
    allow_LCAP_pos_weights=assign_from.allow_pos_LCAP_weights();
    rooted_tree=assign_from.is_rooted_tree();
    fix_zero_brns=assign_from.zero_len_brns_fixed();
    dupl_fix_rate=assign_from.get_dupl_fix_rate();
    dupl_parallel_rate=assign_from.get_dupl_parallel_rate();
    strand_switch_prob=assign_from.get_strand_switch_prob();
    fix_rate_scale=assign_from.get_fix_rate_scale();
    loss_rate_scale=assign_from.get_loss_rate_scale();
    fix_loss_rate=assign_from.get_fix_loss_rate();
    snp_gains_to_losses=assign_from.get_snp_gains_to_losses();
    snp_loss_to_switch=assign_from.get_snp_loss_to_switch();
    use_full_conprobs=assign_from.full_conprobs();
    scale_likelihood=assign_from.likelihood_is_scaled();
    track_gaps_as_missing=assign_from.use_track_gaps_as_missing();
    WGX_depth=assign_from.get_WGX_depth();
    max_breaks=assign_from.get_max_breaks();
    guess_order=assign_from.use_guessed_order();
    use_break_pos=assign_from.count_break_pos();

  if (assign_from.curr_groups != 0) {
	  if(curr_groups !=0) {
		  if (assign_from.curr_groups->get_num_groups() != curr_groups->get_num_groups())
		  {
				delete curr_groups;
				curr_groups=new Amino_acid_group(assign_from.curr_groups->get_num_groups());
		  }
	  }
	  else
		 curr_groups=new Amino_acid_group(assign_from.curr_groups->get_num_groups()); 

	  (*curr_groups)=(*assign_from.curr_groups);
  
  }
  
  if (assign_from.using_codon_position_rates() == TRUE)
    set_use_codon_position_rates();
  
  if (assign_from.using_arbitrary_site_rates() == TRUE)
    {
      sites=new int [this->*num_localities];
      for (i=0; i< this->*num_localities; i++)
		sites[i]=assign_from.get_site_rate_num(i);
       set_use_arbitrary_site_rates(sites);
       delete[] sites;
    }

  if (assign_from.using_generic_site_rate_probs() == TRUE) {
	set_use_generic_site_rates();

	for(i=0; i<num_rates; i++)
		set_generic_site_rate_prob(i, assign_from.get_rate_prob(i));
  }

  if(assign_from.are_optimizing_rate_props() == TRUE)
	  optimize_rate_probs();

  pre_optimize=assign_from.pre_optimizing();
  
    if (assign_from.use_all_taxa()==FALSE ) {
        if (taxa_used ==0) taxa_used=new BOOL [num_taxa];
        
        all_taxa_used=FALSE;
        
        for(i=0; i<num_taxa; i++) taxa_used[i]=assign_from.is_taxa_used(i);
        
    }
  for (i=0; i<assign_from.get_num_rates(); i++) {
    rates[i]=assign_from.get_rate(i);
		  for(j=0; j<num_p_nonsyn; j++) {
			  p_non_syn[j][i]=assign_from.get_p_non_syn(j, i);
			  p_inter_group[j][i]=assign_from.get_p_inter_group(j, i);
			  aa_properties[j][i][get_prop_index_num(CHEM_COMP)]=assign_from.get_aa_prop_fac(j, i, CHEM_COMP);
			  aa_properties[j][i][get_prop_index_num(POLARITY)]=assign_from.get_aa_prop_fac(j, i, POLARITY);
			  aa_properties[j][i][get_prop_index_num(VOLUME)]=assign_from.get_aa_prop_fac(j, i, VOLUME);
			  aa_properties[j][i][get_prop_index_num(ISO_ELEC)]=assign_from.get_aa_prop_fac(j, i, ISO_ELEC);		
			  aa_properties[j][i][get_prop_index_num(HYDROPATHY)]=assign_from.get_aa_prop_fac(j, i, HYDROPATHY);
			  aa_properties[j][i][get_prop_index_num(SCALING)]=assign_from.get_aa_prop_fac(j, i, SCALING);
		  }
	  
  }
	
	
		for(j=0; j<num_p_nonsyn; j++) {
			p_non_syn_used[j]=assign_from.is_p_non_syn_used(j);
			p_inter_group_used[j]=assign_from.is_p_inter_group_used(j);
			for (i=0; i<num_rates; i++) {
				p_non_syn_fix[j][i]=assign_from.p_non_syn_fixed(j,i);
				p_inter_group_fix[j][i]=assign_from.p_inter_group_fixed(j, i);
		  
			}
			for(i=0; i<get_num_aa_props()+1; i++)
				aa_prop_used[j][i]=assign_from.is_aa_prop_used(j,i);
		}

		for(i=0; i<get_num_aa_props()+1; i++)
			aa_property_allowed[i]=assign_from.aa_prop_allowed(i);

		num_live_properties=assign_from.get_num_live_aa_props();
	
  for(i=0; i<4; i++)
    basefreqs[i]=assign_from.return_basefreq(i);
  
  obs_trs_trv=assign_from.get_obs_trs_trv();
  kappa=assign_from.get_trs_trv();

  if (assign_from.using_codon_basefreqs() == TRUE) {
    for(i=0; i<3; i++)
      for(j=0; j<4; j++)
	freqs[i][j]=assign_from. return_codon_basefreq(i,j);
    set_use_codon_basefreqs(freqs); 
  }
  if (assign_from.fixed_basefreq() == TRUE)
    fix_basefreq();

  saved_lnL=assign_from.get_saved_lnL();

#ifdef _OPEN_MP_VERSION_
	num_open_mp_threads=assign_from.get_num_open_mp_threads();
#endif
 return(assign_from); 
}


int Exchange::get_branch_kaks1(int brn_num)                          
{
  if (branch_kaks1 != 0 && brn_num>=0 && brn_num < num_branches)
    return(branch_kaks1[brn_num]);
  else {
    cerr<<"Invalid branch number in get_branch_kaks1"<<endl;
    return(-1);
  }

}



int Exchange::get_prop_index_num(AA_PROPERTIES prop)
{
  switch(prop)
    {
      case CHEM_COMP:
	return(0);
	break;
    case POLARITY:
      return(1);
      break; 
    case VOLUME:
      return(2);
      break; 
    case ISO_ELEC:
      return(3);
      break;
	case HYDROPATHY:
		return(4);
		break;
    case SCALING:
      return(5);
      break;
    }
}




AA_PROPERTIES Exchange::get_prop_num_n(int n)
{
   switch(n)
    {
    case 0:
      return CHEM_COMP;
      break;
    case 1:
      return POLARITY;
      break; 
    case 2:
      return VOLUME;
      break; 
    case 3:
      return ISO_ELEC;
      break;
	case 4:
		return HYDROPATHY;
		break;
    case 5:
    default:
      return SCALING;
      break;
    }
}


int Exchange::get_num_p_non_syn_used()
{
	int i, num=0;	
	
	for(i=0; i<num_p_nonsyn; i++)
		if(p_non_syn_used[i] == TRUE)
			num++;
	return(num);
}


int Exchange::get_num_p_inter_group_used()
{
	int i, num=0;	
	
	for(i=0; i<num_p_nonsyn; i++)
		if(p_inter_group_used[i] == TRUE)
			num++;
	return(num);	
}
 

int Exchange::get_num_aa_props_used(int prop)
{
	int i, num=0;
	for(i=0; i<num_p_nonsyn; i++)
		if(aa_prop_used[i][prop]==TRUE)
			num++;
	return(num);
}


BOOL Exchange::using_codon_position_rates()
{
	if (rate_type == CODON_RATES)
		return(TRUE);
	else
		return(FALSE);
}


BOOL Exchange::using_arbitrary_site_rates()
{
	if (rate_type == ARBITRARY_SITE_RATES)
		return(TRUE);
	else
		return(FALSE);
}

BOOL Exchange::have_site_rate_probs()
{
	if (rate_type == SITE_SPECFIC_RATE_PROB)
		return(TRUE);
	else
		return(FALSE);
}

BOOL Exchange::using_generic_site_rate_probs()
{
	if (rate_type == GENERIC_RATE_PROB)
		return(TRUE);
	else
		return(FALSE);
}

BOOL Exchange::aa_prop_allowed(AA_PROPERTIES prop)		  
{
	return(aa_property_allowed[get_prop_index_num(prop)]);
}


BOOL Exchange::aa_prop_allowed(int prop)				  
{
	return(aa_property_allowed[prop]);
}


BOOL Exchange::is_p_non_syn_used(int pns_num)			  
{
	return(p_non_syn_used[pns_num]);
}

BOOL Exchange::is_p_inter_group_used(int pns_num)		  
{
	return(p_inter_group_used[pns_num]);
}

BOOL Exchange::is_aa_prop_used(int pns_num, AA_PROPERTIES prop)
{
	return(aa_prop_used[pns_num][get_prop_index_num(prop)]);
}


BOOL Exchange::is_aa_prop_used(int pns_num, int prop)
{
	return(aa_prop_used[pns_num][prop]);
}


BOOL Exchange::is_matrix_coeff_used(int pns_num, int matrix)
{
	if (matrix_coeff_used != 0)
		return(matrix_coeff_used[pns_num][matrix]);
    else
        return(FALSE);
}


BOOL Exchange::is_taxa_used(int taxa_num)
{
    if (all_taxa_used==TRUE) return(TRUE);
    else return(taxa_used[abs(taxa_num)%num_taxa]);
    
}
double Exchange::return_basefreq(int freqnum)
{
 if(freqnum>=0 && freqnum<=3)
   return(Exchange::basefreqs[freqnum]);
 else
   return(0);
}  //End Exchange::return_basefreq



double Exchange::return_codon_basefreq(int codon_pos, int base)
{
  if (((codon_pos>=0) && (codon_pos<3)) && ((base>=0) && (base<4)))
    return(codon_basefreqs[codon_pos][base]);
  else
    return(0);
}


double Exchange::get_p_non_syn()                          
{
	return(p_non_syn[0][0]);
}
  
double Exchange::get_p_non_syn(int pns_num)               
{
	return(p_non_syn[pns_num][0]);
}
  

double Exchange::get_p_non_syn(int pns_num, int rate_num)    
{
	return(p_non_syn[pns_num][rate_num]);
}


double Exchange::get_p_inter_group(int pns_num, int rate_num)          
{
	if (!( (p_inter_group_fix[pns_num][rate_num]==TRUE) && (p_inter_group_same[pns_num][rate_num]==TRUE))) 
		return(p_inter_group[pns_num][rate_num]);
	else
		return(p_non_syn[pns_num][rate_num]);
}


double Exchange::get_p_inter_group(int pns_num)			  
{
	return(p_inter_group[pns_num][0]);
}



double Exchange::get_aa_prop_fac(int pns_num, AA_PROPERTIES prop)
{
	if((prop == SCALING) || (aa_property_allowed[get_prop_index_num(prop)] ==TRUE))
		return(aa_properties[pns_num][0][get_prop_index_num(prop)]);
	else
		return(0.0);
}

double Exchange::get_aa_prop_fac(int pns_num, int prop)
{
	if((get_prop_num_n(prop) == SCALING) || (aa_property_allowed[prop] ==TRUE) )
		return(aa_properties[pns_num][0][prop]);
	else
		return(0.0);
}
 

double Exchange::get_aa_prop_fac(int pns_num, int rate_num, AA_PROPERTIES prop)
{
	if((prop == SCALING) || (aa_property_allowed[get_prop_index_num(prop)] ==TRUE))
		return(aa_properties[pns_num][aa_rate_index[rate_num][get_prop_index_num(prop)]][get_prop_index_num(prop)]);
	else
		return(0.0);
}


double Exchange::get_aa_prop_fac(int pns_num, int rate_num, int prop)
{
	if((prop == SCALING) || (aa_property_allowed[prop] ==TRUE))
		return(aa_properties[pns_num][aa_rate_index[rate_num][prop]][prop]);
	else
		return(0.0);
}


double Exchange::get_matrix_coeff(int coeff_num, int pns_num)	
{
	return(matrix_coeff[pns_num][coeff_num]);
}



double Exchange::get_site_rate_prob(int site, int rate) 
{
	return((this->*site_rate_prob_func)(site, rate));
}


void Exchange::set_num_taxa(int taxa)
{
  num_taxa=taxa;
  num_branches=2*num_taxa-1; 
}  //End Exchange::set_num_taxa
  



void Exchange::set_num_sites(int sites)
{
  num_sites=sites;
  num_codons=(int)(sites/3);
}  //End Exchange::set_num_sites




void Exchange::set_num_rates(int rates_requested)
{
  int i, j, k;

  if (rates_requested!=num_rates)
    {
      delete[] rates;
      for(i=0; i<num_p_nonsyn; i++) {
	    delete[] p_non_syn[i];
		delete[] p_inter_group[i];
		for(j=0; j<num_rates; j++) {
			delete[] aa_rate_index[j];
			delete[] aa_properties[i][j];
		}
		delete[] aa_rate_index;

	  }
	  

      rates=new double[rates_requested];
	  aa_rate_index = new int*[rates_requested];
      for(i=0; i<num_p_nonsyn; i++) {
		p_non_syn[i]=new double[rates_requested];
		p_inter_group[i]=new double[rates_requested];
		aa_properties[i] = new double* [rates_requested];
		for (j=0; j<rates_requested; j++) {
			aa_properties[i][j]=new double [get_num_aa_props()+1];
			aa_rate_index[j] = new int [get_num_aa_props()+1];
			for(k=0; k<get_num_aa_props()+1; k++)
				aa_rate_index[j][k]=j;
		}
	  }
	  
	}
  num_rates=rates_requested;
  set_rates();
}  //End Exchange::set_num_rates


void Exchange::use_prob_site_rates(int rates_requested)
{
	int i, j;

	rate_type = SITE_SPECFIC_RATE_PROB;
	if (rates_requested != num_rates) 
		set_num_rates(rates_requested);

	site_rate_prob_func=&Exchange::site_specific_site_rate_prob;

	if (site_rate_probs == 0) {
		site_rate_probs=new double* [get_num_localities()];
	
		for(i=0; i<get_num_localities(); i++)
			site_rate_probs[i]=new double[num_rates];
	}

	for(i=0; i<get_num_localities(); i++)
		for(j=0; j<num_rates; j++)
			site_rate_probs[i][j]=1.0/num_rates;

}



void Exchange::set_p_non_syn(int pns_num, double pns)
{
  if (p_non_syn_fix[pns_num][0]!=TRUE)
    p_non_syn[pns_num][0]=pns;
 else
    cerr<<"Can't set p_non_syn.  Fixed at: "<<p_non_syn[pns_num][0]<<endl;
}  //End Exchange::set_p_non_syn


void Exchange::set_p_non_syn(int pns_num, int rate, double pns)
{
  if (p_non_syn_fix[pns_num][rate]!=TRUE)
    p_non_syn[pns_num][rate]=pns;
 else
    cerr<<"Can't set p_non_syn.  Fixed at: "<<p_non_syn[pns_num][rate]<<endl;
}  //End Exchange::set_p_non_syn




void Exchange::set_p_inter_group(int pns_num, double pgroup)
{ 
  if (p_inter_group_fix[pns_num][0]!=TRUE)
    Exchange::p_inter_group[pns_num][0]=pgroup;
  else
    cerr<<"Can't set p_inter_group.  Fixed at: "<<p_inter_group[pns_num][0]<<endl;
}  //End Exchange::set_p_ns_extra_group


void Exchange::set_p_inter_group(int pns_num, int rate_num, double pgroup)
{ 
  if (p_inter_group_fix[pns_num][rate_num]!=TRUE)
    Exchange::p_inter_group[pns_num][rate_num]=pgroup;
  else
    cerr<<"Can't set p_inter_group.  Fixed at: "<<p_inter_group[pns_num][rate_num]<<endl;
}  //End Exchange::set_p_ns_extra_group


void Exchange::set_aa_prop_fact(int pns_num, AA_PROPERTIES property, double value)
{
	if (allow_LCAP_pos_weights == TRUE)
		aa_properties[pns_num][0][get_prop_index_num(property)]=value;
	else
		aa_properties[pns_num][0][get_prop_index_num(property)]=-fabs(value);
}


void Exchange::set_aa_prop_fact(int pns_num, int property, double value)
{
  if (allow_LCAP_pos_weights == TRUE)
	aa_properties[pns_num][0][property]=value;
  else
	aa_properties[pns_num][0][property]=-fabs(value);
}

void Exchange::set_aa_prop_fact(int pns_num, int rate_num, AA_PROPERTIES property, double value)
{
  if (allow_LCAP_pos_weights == TRUE)
	aa_properties[pns_num][rate_num][get_prop_index_num(property)]=value;
  else
	aa_properties[pns_num][rate_num][get_prop_index_num(property)]=-fabs(value);
}


void Exchange::set_aa_prop_fact(int pns_num, int rate_num, int property, double value)
{
	if (allow_LCAP_pos_weights == TRUE)
		aa_properties[pns_num][rate_num][property]=value;
	else
		aa_properties[pns_num][rate_num][property]=-fabs(value);

}

void Exchange::set_scaling_initial(int pns_num, double p_ns)
{
  aa_properties[pns_num][0][get_prop_index_num(SCALING)]=log(p_ns)+2.15;
  if ((allow_LCAP_pos_weights==FALSE) && (aa_properties[pns_num][0][get_prop_index_num(SCALING)] >0))
	aa_properties[pns_num][0][get_prop_index_num(SCALING)]=-0.1;
}

void Exchange::set_branch_kaks1(int branch)           
{
  int i;

  if (branch_kaks1 == 0) {
    branch_kaks1=new int[num_branches];
    for(i=0; i<num_branches; i++)
      branch_kaks1[i]=0;
  }

  branch_kaks1[branch]=1; 
  is_kaks1=TRUE;
  
}


void Exchange::fix_p_non_syn(double **p_ns)
{
	int i, j;

	for(j=0; j<num_p_nonsyn; j++) {
		for(i=0; i<num_rates; i++) {
			p_non_syn[j][i]=p_ns[j][i]; 
			p_non_syn_fix[j][i]=TRUE;
		}
	}
}

void Exchange::fix_p_non_syn(int pns_num, double p_ns)
{

  p_non_syn_fix[pns_num][0]=TRUE;
  
  p_non_syn[pns_num][0]=p_ns;    
}

	
void Exchange::fix_p_non_syn(int pns_num, int rate, double p_ns)
{		
	p_non_syn_fix[pns_num][rate]=TRUE;
		
	p_non_syn[pns_num][rate]=p_ns;    
}

void Exchange::fix_p_inter_group(int pns_num, double p_intgrp[])
{
  int i;
 
  p_inter_group_fix[pns_num][0]=TRUE;
  p_inter_group_same[pns_num][0]=FALSE;

  for(i=0; i<num_rates; i++)
    p_inter_group[pns_num][i]=p_intgrp[i];  
}

void Exchange::fix_p_inter_group(int pns_num, int rate_num, double p_intgrp)
{
	int i;
	
	p_inter_group_fix[pns_num][rate_num]=TRUE;
	p_inter_group_same[pns_num][rate_num]=FALSE;
	
	p_inter_group[pns_num][rate_num]=p_intgrp;  
}


void Exchange::fix_p_inter_group(int pns_num, int rate_num)
{
	int i;
	
	p_inter_group_fix[pns_num][rate_num]=TRUE;
	
	p_inter_group_same[pns_num][rate_num]=TRUE;
}


void Exchange::set_num_p_nonsyn(int num, double **pns)
{
  int i, j;

  
  if (num_p_nonsyn != num) {
	if (p_non_syn !=0) 
		{

			if(matrix_coeff != 0) {
				
				for(i=0; i<num_p_nonsyn; i++) {
					delete[] matrix_coeff[i];
					delete[] matrix_coeff_used[i];
				}
				delete[] matrix_coeff;
				delete[] matrix_coeff_used;
			}
  

			delete[] p_non_syn_used;
			delete[] p_inter_group_used;
			

			for(i=0; i<num_p_nonsyn; i++)
				delete[] aa_prop_used[i];
			delete[] aa_prop_used;

			for(i=0; i<num_p_nonsyn; i++) {
				delete[] p_non_syn[i];
				delete[] p_non_syn_fix[i];
				delete[] p_inter_group_fix[i];
				delete[] p_inter_group_same[i];
				delete[] p_inter_group[i];
				for(j=0; j<num_rates; j++)
					delete[] aa_properties[i][j];
				delete[] aa_properties[i];
			}
			delete[] p_non_syn;
			delete[] p_non_syn_fix;
			delete[] p_inter_group_fix;
			delete[] p_inter_group_same;
			delete[] p_inter_group;
			delete[] aa_properties;
		}
		
		num_p_nonsyn=num;

		p_non_syn=new double*[num_p_nonsyn];
		p_inter_group=new double* [num_p_nonsyn];
		aa_properties=new double**[num_p_nonsyn];
		p_non_syn_fix=new BOOL*[num_p_nonsyn];
		p_inter_group_fix=new BOOL* [num_p_nonsyn];
		p_inter_group_same=new BOOL* [num_p_nonsyn];
		p_non_syn_used=new BOOL [num_p_nonsyn];
		p_inter_group_used=new BOOL [num_p_nonsyn];
		aa_prop_used=new BOOL* [num_p_nonsyn];

		for(i=0; i<num_p_nonsyn; i++)
		{
			aa_prop_used[i]=new BOOL [get_num_aa_props()+1];
			for(j=0; j<get_num_aa_props()+1; j++)
				aa_prop_used[i][j]=TRUE;
			p_non_syn_used[i]=TRUE;
			p_inter_group_used[i]=TRUE;
			p_non_syn_fix[i]=new BOOL [num_rates];
			p_inter_group_fix[i]=new BOOL [num_rates];
			p_inter_group_same[i]=new BOOL [num_rates];
			

			p_non_syn[i]=new double [num_rates];
			p_inter_group[i]=new double [num_rates];
			aa_properties[i]=new double* [num_rates];
			for(j=0; j<num_rates; j++) 
				aa_properties[i][j]=new double[get_num_aa_props()+1];
			for(j=0; j<num_rates; j++) {
				p_non_syn[i][j]=pns[i][j];
				p_non_syn_fix[i][j]=FALSE;
				p_inter_group_fix[i][j]=FALSE;
				p_inter_group_same[i][j]=TRUE;
			}
		
		}
  

		if(num_matrices != 0) {
			matrix_coeff=new double* [num_p_nonsyn];
			matrix_coeff_used=new BOOL* [num_p_nonsyn];
			for(i=0; i<num_p_nonsyn; i++){
				matrix_coeff_used[i]=new BOOL[num_matrices];;
				matrix_coeff[i]=new double[num_matrices];
				for(j=0; j<num_matrices; j++) {
					matrix_coeff[i][j]=1.0;
					matrix_coeff_used[i][j]=TRUE;
				}
			}
		}

	}

}


void Exchange::set_num_p_nonsyn(int num)
//Warning: Only use this function if you don't need
//to specific Ka/Ks ratios at the point of calling
{
  int i, j, k;
  double pns, ping, lcapvals[6], *matrix_vals=0;

	//cout<<"CURR PNS: "<<num_p_nonsyn<<" NEW: "<<num<<" RATES "<<num_rates<<endl;
   if(num_p_nonsyn != num) {
		if (p_non_syn !=0) 
		{
			matrix_vals=new double [num_matrices];
			pns=p_non_syn[0][0];
			ping=p_inter_group[0][0];
			for(i=0; i<get_num_aa_props()+1; i++)
				lcapvals[i]=aa_properties[0][0][i];

			for(i=0; i<num_matrices; i++)
				matrix_vals[i]=matrix_coeff[0][i];

			if(matrix_coeff != 0) {
				
				for(i=0; i<num_p_nonsyn; i++) {
					delete[] matrix_coeff[i];
					delete[] matrix_coeff_used[i];
				}
				delete[] matrix_coeff;
				delete[] matrix_coeff_used;
			}


			

			delete[] p_non_syn_used;
			delete[] p_inter_group_used;
			for(i=0; i<num_p_nonsyn; i++)
				delete[] aa_prop_used[i];
			delete[] aa_prop_used;

			for(i=0; i<num_p_nonsyn; i++) {
				delete[] p_non_syn[i];
				delete[] p_non_syn_fix[i];
				delete[] p_inter_group_fix[i];
				delete[] p_inter_group_same[i];
				delete[] p_inter_group[i];
				for(j=0; j<num_rates; j++)
					delete[] aa_properties[i][j];
				delete[] aa_properties[i];
			}
			delete[] p_non_syn;
			delete[] p_non_syn_fix;
			delete[] p_inter_group_fix;
			delete[] p_inter_group_same;
			delete[] p_inter_group;
			delete[] aa_properties;
		}

	  num_p_nonsyn=num;
	  
	  p_non_syn=new double*[num_p_nonsyn];
	  p_inter_group=new double* [num_p_nonsyn];
	  aa_properties=new double**[num_p_nonsyn];
	  p_non_syn_fix=new BOOL* [num_p_nonsyn];
	  p_inter_group_fix=new BOOL* [num_p_nonsyn];
	  p_inter_group_same=new BOOL*[num_p_nonsyn];
	  p_non_syn_used=new BOOL [num_p_nonsyn];
	  p_inter_group_used=new BOOL [num_p_nonsyn];
	  aa_prop_used=new BOOL* [num_p_nonsyn];

		for(i=0; i<num_p_nonsyn; i++)
		{

			aa_prop_used[i]=new BOOL [get_num_aa_props()+1];
			for(j=0; j<get_num_aa_props()+1; j++)
				aa_prop_used[i][j]=TRUE;
			p_non_syn_used[i]=TRUE;
			p_inter_group_used[i]=TRUE;
				
			p_non_syn_fix[i]= new BOOL[num_rates];
			p_inter_group_fix[i]=new BOOL[num_rates];
			p_inter_group_same[i]=new BOOL[num_rates];
			p_non_syn[i]=new double [num_rates];
			p_inter_group[i]=new double [num_rates];
			aa_properties[i]=new double* [num_rates];
	 
	  
			for(j=0; j<num_rates; j++) {
				aa_properties[i][j]=new double[get_num_aa_props()+1];
				p_non_syn[i][j]=pns;
				p_non_syn_fix[i][j]=FALSE;
				p_inter_group_fix[i][j]=FALSE;
				p_inter_group_same[i][j]=TRUE;
				p_inter_group[i][j]=ping;
				for(k=0; k<get_num_aa_props()+1; k++)
					aa_properties[i][j][k]=lcapvals[k];
			}
	  
		}
		if(num_matrices != 0) {
			matrix_coeff=new double* [num_p_nonsyn];
			matrix_coeff_used=new BOOL* [num_p_nonsyn];
			for(i=0; i<num_p_nonsyn; i++){
				matrix_coeff_used[i]=new BOOL[num_matrices];;
				matrix_coeff[i]=new double[num_matrices];
				for(j=0; j<num_matrices; j++) {
					matrix_coeff[i][j]=matrix_vals[j];
					matrix_coeff_used[i][j]=TRUE;
				}
			}
		}
	

	}
  	if(matrix_vals !=0)
		delete[] matrix_vals;

	//cout<<"MADE "<<p_non_syn<<endl;
}


void Exchange::set_num_p_nonsyn(int n_pns, int n_patterns)
{
	num_nonsyn_patterns=n_patterns;
	num_patterns=&Exchange::num_nonsyn_patterns;
	set_num_p_nonsyn(n_pns);
}


void Exchange::set_num_matrices(int num) 
{
	int i, j;

	if (matrix_coeff != 0) {
		for (i=0; i<num_p_nonsyn; i++) {
			delete[] matrix_coeff[i];
			delete[] matrix_coeff_used[i];
		}
		delete[] matrix_coeff;
		delete[] matrix_coeff_used;
	}
	num_matrices=num;

	if (num_matrices != 0) {
		matrix_coeff=new double*[num_p_nonsyn];
		matrix_coeff_used=new BOOL*[num_p_nonsyn];
		for(i=0; i<num_p_nonsyn; i++) {
			matrix_coeff[i]=new double[num_matrices];
			matrix_coeff_used[i]=new BOOL[num_matrices];
			for(j=0; j<num_matrices; j++) {
				matrix_coeff[i][j]=1.0;
				matrix_coeff_used[i][j]=TRUE;
			}
		}
	  }
}


void Exchange::set_model(LKMODEL curmodel)
{
	int i;
	current_model=curmodel;

  switch(current_model)
    {
    case GG_98_NUCLEOTIDE:
    case MG_94_GG_98:
    case C_00_GG_98:
    case LCAP_GG_98:
      set_branch_basefreqs(TRUE);
      break;
    default:
        break;
    }

  switch(current_model)
    {
    case LCAP_JC:
    case LCAP_K2P:
    case LCAP_HKY:
    case LCAP_GG_98:
      set_aa_class_modeling(TRUE);
	  standard_model_nonsyn_params=get_num_aa_props()+1;
      break;
    default:
        break;
    }

  switch(current_model)
    {
    case MG_94_JC:
    case MG_94_K2P:
    case MG_94_HKY:
    case MG_94_GG_98:
		standard_model_nonsyn_params=1;
		break;
    case C_00_JC:
    case C_00_K2P:
    case C_00_HKY:
    case C_00_GG_98:
		standard_model_nonsyn_params=2;
		break;
	case AAMULMAT_JC:
	case AAMULMAT_K2P:
	case AAMULMAT_HKY:
		num_nonsyn_params=&Exchange::num_matrices;
    default:
        break;
	}
  switch(current_model)
    {
    case MG_94_JC:
    case MG_94_K2P:
    case MG_94_HKY:
    case MG_94_GG_98:
    case C_00_JC:
    case C_00_K2P:
    case C_00_HKY:
    case C_00_GG_98:
    case LCAP_JC:
    case LCAP_K2P:
    case LCAP_HKY:
    case LCAP_GG_98:
	case AAMULMAT_JC:
	case AAMULMAT_K2P:
	case AAMULMAT_HKY:
    case BPS_JC:
    case BPS_K2P:
    case BPS_HKY:
      model_codons=TRUE;
      Exchange::num_localities=&Exchange::num_codons;
      condlike_size=64;
      break;
    default:
        break;
    }


  switch(current_model)
  {
	case AAMULMAT_JC:
	case AAMULMAT_K2P:
	case AAMULMAT_HKY:
		set_num_matrices(1);
		break;
    default:
      break;
  }
  switch (current_model)
  {
  case DUPL:
  case DUPL_FIX:
  case DUPL_PARALLEL:
  case DUPL_PARALLEL_2_RATE:
  case DUPL_PARALLEL_2_RATE_NOSTATE:
  case DUPL_PARALLEL_FIX:
  case DUPL_PARALLEL_FIX_SUBF:
  case DUPL_SUBF_3_RATE:
  case DUPL_SUBF_3_RATE_NOSTATE:
  case DUPL_SUBF_3_RATE_ALLSTATES_NOSTATE:
  case DUPL_SUBF_ONLY:
  case DUPL_SUBF_ONLY_NOSTATE:
  case DUPL_NOSTATE:
  case DUPL_FIX_NOSTATE:
  case DUPL_PARALLEL_NOSTATE:
  case DUPL_PARALLEL_FIX_NOSTATE:
  case DUPL_PARALLEL_FIX_SUBF_NOSTATE:
  case DUPL_SUBF_ALL_STATES:
  case DUPL_SUBF_ALL_STATES_NOSTATE:
  case DUPL_2_RATE_NOSUBF:
  case DUPL_2_RATE_NOSUBF_NOSTATE:
  case DUPL_2_RATE_NOSUBF_ALLSTATES_NOSTATE:
  case DUPL_SLOW_LOSS_CONV_FIX:
  case DUPL_SLOW_LOSS_CONV_FIX_NOSTATE:
	  condlike_size=14;
  case  DUPL_ARBITRARY:
	  rooted_tree=TRUE;
          cerr<<"ERROR: Setting invalid condlike size for DUPL_ARBITRARY model\n";
          condlike_size=6;
	  break;
  case SNP_EVOL_MODEL:
	  rooted_tree=TRUE;
	  condlike_size=4;
	  use_full_conprobs=TRUE;
	  scale_likelihood=TRUE;
	  break;
    default:
      break;
  }

}  //End Exchange::set_model

void Exchange::set_model(LKMODEL curmodel, int num_states)
{
    int i;
    current_model=curmodel;
    
    
    if (current_model==  DUPL_ARBITRARY){
        rooted_tree=TRUE;
        condlike_size=num_states;
    }
    else
        cerr<<"ERROR set_model with num_states only valid for DUPL_ARBITRARY model\n";
    

    
}  //End Exchange::set_model



void Exchange::set_use_codon_position_rates()
{
  if (num_rates< 3) 
    cerr<<"Can't use codon position rates with less than three subsitition rates\n";
  else if (model_codons == TRUE) 
    cerr<<"Can't use codon position rates with codon models of evolution\n";
  else { 
    rate_type=CODON_RATES;
    
    site_rate_num_func=&Exchange::codon_positions_site_rate_num;
	site_rate_func=&Exchange::codon_positions_site_rates;
	site_rate_prob_func=&Exchange::codon_position_site_rate_prob;
  }
}


void Exchange::set_use_arbitrary_site_rates(int *initial_rates)
{ 
  int i;

  if (num_rates>= 2) {
    rate_type = ARBITRARY_SITE_RATES;
    
	site_rate_num_func=&Exchange::arbitrary_site_rate_num;
    site_rate_func=&Exchange::arbitrary_site_rates;
	site_rate_prob_func=&Exchange::arbitrary_site_rate_prob;

    site_rates =new int[this->*num_localities];
    for (i=0; i<this->*num_localities; i++)
      site_rates[i]= initial_rates[i];

  }
  else {
    cerr<<"Can't use arbitrary site rates with less than two subsitition rates\n";
  }
  
}


void Exchange::set_use_generic_site_rates()
{
	int i;
	double prob;

	rate_type = GENERIC_RATE_PROB;

	site_rate_prob_func=&Exchange::generic_rate_prob;

	if (rate_probs != 0)
		delete[] rate_probs;

	rate_probs = new double [num_rates];

	prob = 1.0/(double)num_rates;

	for(i=0; i<num_rates; i++)
		rate_probs[i]=prob;

}


void Exchange::set_generic_site_rate_prob(int rate, double prob)
{
	if ((rate >=0) && (rate<num_rates))
		rate_probs[rate]=prob;
	else
		cerr<<"Error: invalid request to set rate probabilitiy\n";
}

void Exchange::set_site_rate(int site, int rate_num)
{
  if ((site>=0) && (site < this->*num_localities) && (rate_num >=0 ) && (rate_num <num_rates))
    site_rates[site]=rate_num;
  else
    cerr<<"Error setting rate "<<rate_num<<" for site "<<site<<endl;
}



void Exchange::set_use_codon_basefreqs(double newfreqs[3][4])
{
  int i,j;
  codon_freqs=TRUE;

  for (i=0; i<3; i++)
    for(j=0;j<4; j++)
      codon_basefreqs[i][j]=newfreqs[i][j];
}


void Exchange::disable_property(AA_PROPERTIES prop)
{
	if (aa_property_allowed[get_prop_index_num(prop)]!=FALSE)
		num_live_properties--;
	aa_property_allowed[get_prop_index_num(prop)]=FALSE;
}


void Exchange::set_p_non_syn_used(int pns_num, BOOL val)
{
	p_non_syn_used[pns_num]=val;
}


void Exchange::set_p_inter_group_used(int pns_num, BOOL val)
{
	p_inter_group_used[pns_num]=val;
}


void Exchange::set_aa_property_used(int pns_num, AA_PROPERTIES prop, BOOL val)
{
	aa_prop_used[pns_num][get_prop_index_num(prop)]=val;
}
  

void Exchange::set_aa_property_used(int pns_num, int prop, BOOL val)
{
	aa_prop_used[pns_num][prop]=val;
}

void Exchange::set_matrix_coeff_used(int pns_num, int matrix_num, BOOL val)
{
	if (num_matrices != 0)
		matrix_coeff_used[pns_num][matrix_num]=val;

}


void Exchange::set_site_rate_prob(int site, int rate, double rate_prob)
{
	if ((0<=site) && (get_num_localities() > site) && (0<=rate) && (num_rates> rate) && 
		(0<=rate_prob) && (rate_prob<=1.0))
			site_rate_probs[site][rate]=rate_prob;
	else
		cerr<<"Error setting site rate prob.\n";

}


void Exchange::use_single_property(AA_PROPERTIES prop)
{
	use_single_property(get_prop_index_num(prop));
}


void Exchange::use_single_property(int prop_num)
{
	int i;
	for(i=1; i<num_rates; i++)
		aa_rate_index[i][prop_num]=0;
}


void Exchange::use_only_hetero_internal()
{
	num_unallowed_SNP_states=3;
	unallowed_states[0]=SNP_ABSENT;
	unallowed_states[1]=TYPE_A;
	unallowed_states[2]=TYPE_B;
}

void Exchange::mark_taxa_unused(int taxa_num)
{
    int i;
    all_taxa_used=FALSE;
    
    if (taxa_used==0) {
        taxa_used=new BOOL[num_taxa];
        
        for(i=0; i<num_taxa; i++) taxa_used[i]=TRUE;
        
    }
    
    taxa_used[taxa_num]=FALSE;
}

Exchange::~Exchange()
{
  int i, j;
  for(j=0; j<num_p_nonsyn; j++){
	  delete[] aa_prop_used[j];
	  
	  for (i=0; i<num_rates; i++) 
		delete[] aa_properties[j][i];

		delete[] aa_properties[j];
	 
	  delete[] p_non_syn_fix[j];
	  delete[] p_inter_group_fix[j];
	  delete[] p_inter_group_same[j];
		delete[] p_non_syn[j];
		delete[] p_inter_group[j];
  }
  delete[] p_non_syn;
  delete[] p_inter_group;
  delete[] aa_properties;
  delete[] rates;	 
  delete[] p_non_syn_used;
  delete[] p_inter_group_used;
  delete[] aa_prop_used;
  delete[] p_non_syn_fix;
  delete[] p_inter_group_fix;
  delete[] p_inter_group_same;
  
  
  for(i=0; i<num_rates; i++)
	delete[] aa_rate_index[i];

  delete[] aa_rate_index;

  if (site_rates != 0)
    delete[] site_rates;

  if(rate_probs != 0)
	  delete[] rate_probs;

  if(site_rate_probs != 0)
  {
	  for(i=0; i<get_num_localities(); i++)
		  delete[] site_rate_probs[i];
	  delete[] site_rate_probs;
  }



  if (num_matrices !=0)
  {
	for(i=0; i<num_p_nonsyn; i++)
	{
		delete[] matrix_coeff[i];
		delete[] matrix_coeff_used[i];
	}
	delete[] matrix_coeff;
	delete[] matrix_coeff_used;

  }
    
    if( taxa_used!=0) delete[] taxa_used;

}



//Private Functions
void Exchange::set_rates()
{
  int i, j, k;

	//cout<<"NR: "<<num_rates<<" NPS: "<<num_p_nonsyn<<" OBJ: "<<p_non_syn<<endl;
	
  for (i=0; i<num_rates; i++)
    {
      rates[i]=1;
      for(j=0; j< num_p_nonsyn; j++) {
		p_non_syn[j][i]=1.0;
		p_inter_group[j][i]=1.0;
	
		//Note that GCC didn't estimate HYDROPATHY initial value
		aa_properties[j][i][get_prop_index_num(CHEM_COMP)]=0.059;
		aa_properties[j][i][get_prop_index_num(POLARITY)]=-0.069;
		aa_properties[j][i][get_prop_index_num(VOLUME)]=-0.223;
		aa_properties[j][i][get_prop_index_num(ISO_ELEC)]=-0.065;
		aa_properties[j][i][get_prop_index_num(HYDROPATHY)]=-0.1;
		aa_properties[j][i][get_prop_index_num(SCALING)]=0.0;
		if (allow_LCAP_pos_weights==FALSE)
			aa_properties[j][i][get_prop_index_num(CHEM_COMP)]=0.0;
    }
  }
}  //End Exchange::set_rates


int Exchange::codon_positions_site_rate_num(int site)
{
	return(site%3);
}
 


int Exchange::arbitrary_site_rate_num(int site)
{
if ((site >= 0) && (site<this->*num_localities))
    return(site_rates[site]);
  else
    cerr<<"Invalid site in call to Exchange::arbitrary_site_rate_num()\n";
    return(-1);

}



int Exchange::get_rate_aa_prop_rate(int rate_num, AA_PROPERTIES prop)
{
	return(get_rate_aa_prop_rate(rate_num, get_prop_index_num(prop)));
}


int Exchange::get_rate_aa_prop_rate(int rate_num, int prop_num)
{
	return(aa_rate_index[rate_num][prop_num]);
}



double Exchange::codon_positions_site_rates(int site)
{
  if (site%3 == 0)
    return(rates[0]);
  else if (site%3 == 1)
    return(rates[1]);
  else
    return(rates[2]);
}

double Exchange::arbitrary_site_rates(int site)
{
  if ((site >= 0) && (site<this->*num_localities))
    return(rates[site_rates[site]]);
  else
    cerr<<"Invalid site in call to Exchange::arbitrary_site_rates()\n";
    return(0.0);
}


		

double Exchange::arbitrary_site_rate_prob(int site, int rate)
{
	if (rate == arbitrary_site_rate_num(site))
		return(1.0);
	else
		return(0.0);
}
  

double Exchange::codon_position_site_rate_prob(int site, int rate)
{
	if (rate == codon_positions_site_rate_num(site))
		return(1.0);
	else
		return(0.0);
}
  
double Exchange::singe_rate_site_rate_prob(int site, int rate)
{
	return(1.0);
}


double Exchange::site_specific_site_rate_prob(int site, int rate)
{
		return(site_rate_probs[site][rate]);
}
  

double Exchange::generic_rate_prob(int site, int rate)
{
	return(rate_probs[rate]);
}












