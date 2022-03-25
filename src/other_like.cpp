//Copyright 2005-2006 Gavin Conant

#include <iostream>
#include <math.h>
#include <iomanip>
#include <fstream>
#include "other_like.h"
#include "powell.h"

#ifdef _OPEN_MP_VERSION_
#include "omp.h"
#endif

using namespace::std;

//#define DO_TIME
//#define MAC_TIME

#ifdef DO_TIME
#include <sys/time.h>
#include <sys/resource.h>

unsigned long RunTime() 
{
	unsigned long retval;
	struct timespec t1;
	
	clock_gettime(CLOCK_MONOTONIC,  &t1);
	
	retval = (t1.tv_sec*1e9 + t1.tv_nsec)/1000000;
	return(retval);
}
#endif

#ifdef MAC_TIME
#include <mach/mach_time.h>


#endif

int cmp(const void *x, const void *y)
{
  double xx = *(double*)x, yy = *(double*)y;
  if (xx < yy) return -1;
  if (xx > yy) return  1;
  return 0;
}

TRANSPOINT_STATES int_to_transpoint_state(int n)
{
	switch(n) {
	case 0:
		return(DUPLICATE);
	case 1:
		return(DUPLICATE_TO_SINGLE);
	default:
		return(SINGLE);
	}
}


int transpoint_state_to_int(TRANSPOINT_STATES state)
{
	switch (state) {
	case DUPLICATE:
		return(0);
	case DUPLICATE_TO_SINGLE:
		return(1);
	case SINGLE:
		return(2);
	}

}


BOOL is_trans_in_transpoint_state(DUPL_LOSS_STATES starting_state, DUPL_LOSS_STATES ending_state, 
								  TRANSPOINT_STATES the_state)
{
	switch (the_state) {
	case DUPLICATE:
		switch (starting_state) {
		case BOTH_PRESENT:
		case BOTH_PRESENT_FIXED:
		case BOTH_1_BIAS:
		case BOTH_2_BIAS:
		case BOTH_FIXED_SUBF:
			switch (ending_state) {
			case BOTH_PRESENT:
			case BOTH_PRESENT_FIXED:
			case BOTH_1_BIAS:
			case BOTH_2_BIAS:
			case BOTH_FIXED_SUBF:
				return((BOOL)TRUE);
			default:
				return((BOOL)FALSE);
			}
			break;
			default:
				return((BOOL)FALSE);
		}
		break;
	case DUPLICATE_TO_SINGLE:
		switch (starting_state) {
		case BOTH_PRESENT:
		case BOTH_PRESENT_FIXED:
		case BOTH_1_BIAS:
		case BOTH_2_BIAS:
		case BOTH_FIXED_SUBF:
			switch(ending_state) {
			case COPY1:
			case COPY2:
			case COPY1_BIAS:
			case COPY2_BIAS:
				return((BOOL)TRUE);
			default:
				return((BOOL)FALSE);
			}
			break;
		default:
			return((BOOL)FALSE);
		}
		break;
	case SINGLE:
		switch(starting_state) {
			case COPY1:
			case COPY2:
			case COPY1_BIAS:
			case COPY2_BIAS:
				switch(ending_state) {
				case COPY1:
				case COPY2:
				case COPY1_BIAS:
				case COPY2_BIAS:
					return((BOOL)TRUE);
				default:
					return((BOOL)FALSE);
				}
			default:
				return((BOOL)FALSE);
			}
		break;

	}


}



double Dupl_Base_model::root_freq(int site)
{
	//We assume that all genes are identical and in state BOTH_PRESENT at the
	//root of the tree
	switch (dupl_to_loss_state(site)) {
	case LOST:
		return(0.0);
		break;
	case COPY1:
		return(0.0);
		break;
	case COPY2:
		return(0.0);
		break;
	case BOTH_PRESENT:
		return(1.0);
		break;
	case GENERIC_SINGLE_COPY:
	case MISSING:
	case BOTH_PRESENT_FIXED:
	case BOTH_1_BIAS:
	case BOTH_2_BIAS:
	case COPY1_OR_BOTH:
	case COPY2_OR_BOTH:
	case COPY1_BIAS:
	case COPY2_BIAS:
	case BOTH_FIXED_SUBF:
		return(0.0);
		break;
	default:
		return(0.0);
		break;
	}

}

void Dupl_Base_model::initialize_arrays()
{
}
 

void Dupl_Base_model::get_site_state_probs(double **&prob_array, int taxa_id)
{
	cerr<<"Error: get_site_state_probs can only be called is the model is Dupl_Subf_All_State_model\n";
}


void Dupl_Base_model::partial_prob_w_rate(int locale, Branch *lsib, int rate_num, Branch *stop_id)
 //The heart of the likelihood calculation: uses the tranisition probablity matrices stored in
  //the tree to calculate the likelihood of part of that tree.  Uses post-order tree transversal.
  //This is modified version for the Dupl_models that allows the tip states to be ambiguous--hence
  //probability is summed across those ambiguous states
{
  int i, j, k;
  double temp_p1, temp_p2;
  Branch *rsib, *parent;
  
  rsib=lsib->get_sibling();
  parent=lsib->get_parent();
  

  
  if (lsib->is_tip()==(BOOL)FALSE)
    {
      if (rsib->is_tip()==(BOOL)FALSE) 
		{
			partial_prob_w_rate(locale, curr_tree->find_left_tip(rsib), rate_num, rsib);
			for(i=0; i<curr_exchange->get_condlike_size(); i++)
			{
				temp_p1=temp_p2=0;
				for(j=0; j<curr_exchange->get_condlike_size(); j++)
				{
#ifdef _OPEN_MP_VERSION_
					temp_p1+=lsib->get_cond_prob_locale(omp_get_thread_num(), j)*lsib->get_trpb(rate_num, i, j);
					temp_p2+=rsib->get_cond_prob_locale(omp_get_thread_num(), j)*rsib->get_trpb(rate_num, i, j);
#else
					temp_p1+=lsib->get_cond_prob(j)*lsib->get_trpb(rate_num, i, j);
					temp_p2+=rsib->get_cond_prob(j)*rsib->get_trpb(rate_num, i, j);
#endif
				} 
#ifdef _OPEN_MP_VERSION_
				parent->set_cond_prob_locale(omp_get_thread_num(), i, temp_p1*temp_p2);
#else
				parent->set_cond_prob(i, temp_p1*temp_p2);
#endif
			}	   
		}
    
      else
		{
			for(i=0; i<curr_exchange->get_condlike_size(); i++)
			{
				temp_p1=0;
				for(j=0; j<curr_exchange->get_condlike_size(); j++) 
#ifdef _OPEN_MP_VERSION_
					temp_p1+=lsib->get_cond_prob_locale(omp_get_thread_num(),j)*lsib->get_trpb(rate_num, i, j);
#else
					temp_p1+=lsib->get_cond_prob(j)*lsib->get_trpb(rate_num, i, j);
#endif
		  
#ifdef _OPEN_MP_VERSION_
				parent->set_cond_prob_locale(omp_get_thread_num(), i, temp_p1*get_tip_prob(locale, rsib, rate_num, i));
#else
				parent->set_cond_prob(i, temp_p1*get_tip_prob(locale, rsib, rate_num, i));
#endif
			}
		}
    }

  else
    {
      if (rsib->is_tip()==(BOOL)FALSE) 
		{
			partial_prob_w_rate(locale, curr_tree->find_left_tip(rsib), rate_num, rsib);
	
			for(i=0; i<curr_exchange->get_condlike_size(); i++)
			{
				temp_p2=0;
				
				for(j=0; j<curr_exchange->get_condlike_size(); j++)
#ifdef _OPEN_MP_VERSION_
					temp_p2+=rsib->get_cond_prob_locale(omp_get_thread_num(), j)*rsib->get_trpb(rate_num, i, j);
#else
					temp_p2+=rsib->get_cond_prob(j)*rsib->get_trpb(rate_num, i, j); 
#endif
#ifdef _OPEN_MP_VERSION_
				parent->set_cond_prob_locale(omp_get_thread_num(), i, get_tip_prob(locale, lsib, rate_num, i)*temp_p2);
#else
				parent->set_cond_prob(i, get_tip_prob(locale, lsib, rate_num, i)*temp_p2); 
#endif
			}   
		}

      else
		  for(i=0; i<curr_exchange->get_condlike_size(); i++) 
#ifdef _OPEN_MP_VERSION_
			parent->set_cond_prob_locale(omp_get_thread_num(), i, get_tip_prob(locale, lsib, rate_num, i)*get_tip_prob(locale, rsib, rate_num, i));
#else
			parent->set_cond_prob(i, get_tip_prob(locale, lsib, rate_num, i)*get_tip_prob(locale, rsib, rate_num, i));
#endif
	  
    }

  if (parent!=stop_id)
      partial_prob_w_rate(locale, parent, rate_num, stop_id);
 
}  //End Like_model::partial_prob_w_rate




void Dupl_Base_model::partial_prob_w_rate_nonhidden(int locale, Branch *lsib, int rate_num, Branch *stop_id, int nonhidden_taxa)
 //The heart of the likelihood calculation: uses the tranisition probablity matrices stored in
  //the tree to calculate the likelihood of part of that tree.  Uses post-order tree transversal.
  //This is modified version for the Dupl_models that allows the tip states to be ambiguous--hence
  //probability is summed across those ambiguous states
{
  int i, j, k;
  double temp_p1, temp_p2;
  Branch *rsib, *parent;
  
  rsib=lsib->get_sibling();
  parent=lsib->get_parent();
  

  
  if (lsib->is_tip()==(BOOL)FALSE)
    {
      if (rsib->is_tip()==(BOOL)FALSE) 
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
		  
				if (rsib->get_taxa_id() != nonhidden_taxa)
					parent->set_cond_prob(i, temp_p1*get_tip_prob(locale, rsib, rate_num, i)); 
				else
					parent->set_cond_prob(i, temp_p1*rsib->get_trpb(rate_num, i, 
						(*curr_data)[rsib->get_taxa_id()][locale]));
			}
		}
    }

  else
    {
      if (rsib->is_tip()==(BOOL)FALSE) 
		{
			partial_prob_w_rate(locale, curr_tree->find_left_tip(rsib), rate_num, rsib);
	
			for(i=0; i<curr_exchange->get_condlike_size(); i++)
			{
				temp_p2=0;
				
				for(j=0; j<curr_exchange->get_condlike_size(); j++)
					temp_p2+=rsib->get_cond_prob(j)*rsib->get_trpb(rate_num, i, j); 
				if (lsib->get_taxa_id() != nonhidden_taxa)
					parent->set_cond_prob(i, get_tip_prob(locale, lsib, rate_num, i)*temp_p2); 
				else
					parent->set_cond_prob(i, temp_p2*lsib->get_trpb(rate_num, i, (*curr_data)[lsib->get_taxa_id()][locale]));
			}   
		}

      else
		  for(i=0; i<curr_exchange->get_condlike_size(); i++) {
			  if ((lsib->get_taxa_id() != nonhidden_taxa) && (rsib->get_taxa_id() != nonhidden_taxa)) 
				parent->set_cond_prob(i, get_tip_prob(locale, lsib, rate_num, i)*get_tip_prob(locale, rsib, rate_num, i)); 
			  else if (lsib->get_taxa_id() == nonhidden_taxa)
				parent->set_cond_prob(i, lsib->get_trpb(rate_num, i, (*curr_data)[lsib->get_taxa_id()][locale])*
					get_tip_prob(locale, rsib, rate_num, i)); 
			  else if (rsib->get_taxa_id() == nonhidden_taxa)
				parent->set_cond_prob(i, rsib->get_trpb(rate_num, i, (*curr_data)[rsib->get_taxa_id()][locale])*
					get_tip_prob(locale, lsib, rate_num, i));
		  }
	  
    }

  if (parent!=stop_id)
      partial_prob_w_rate(locale, parent, rate_num, stop_id);
 
}  //End Like_model::partial_prob_w_rate





double Dupl_Base_model::get_tip_prob(int locale, Branch *taxa, int rate_num, int cond_state)
{	
	int i;
	double ret_val=0;

	for(i=0; i<state_redundancy_val(dupl_to_loss_state((*curr_data)[taxa->get_taxa_id()][locale]));
					i++)
					ret_val+=
						taxa->get_trpb(rate_num, cond_state, 
						loss_state_to_dupl(get_redund_position_n(dupl_to_loss_state((*curr_data)[taxa->get_taxa_id()][locale]), i)));

	return(ret_val);	
}


BOOL Dupl_Base_model::allowed_state(DUPL_LOSS_STATES the_state)
{
	BOOL retval=(BOOL)FALSE;

	switch (the_state) {
	case BOTH_PRESENT:
	case COPY1:
	case COPY2:
	case COPY1_OR_BOTH:
	case COPY2_OR_BOTH:
	case GENERIC_SINGLE_COPY:
	case MISSING:
	case LOST:
		retval=(BOOL)TRUE;
		break;

	}
	return(retval);
}


double Dupl_Base_model::calculate_transpoint_prob(Branch *taxa, TRANSPOINT_STATES state, int locale, int rate_num)
{
	//Note that this function currently does not allow for ambiguous tip duplication states
	int i, j, myid=0;
	double retval=0, new_prob, sib_prob;
	Branch *curr;

#ifdef _OPEN_MP_VERSION_
	myid=omp_get_thread_num();
#endif
	
	for(i=0; i<curr_exchange->get_condlike_size(); i++)
#ifdef _OPEN_MP_VERSION_
		transpoint_condprobs1_thread[myid][i]=0.0;
#else
		transpoint_condprobs1[i]=0.0;
#endif
	curr=taxa;
	if (taxa != curr_tree->find_root()) {
		
		//Make conditional probs for this branch with only the allowed transitions
		if (taxa->is_tip() == (BOOL)TRUE) {
			for(i=0; i<curr_exchange->get_condlike_size(); i++)
			{
				if (is_trans_in_transpoint_state(dupl_to_loss_state(i), dupl_to_loss_state((*curr_data)[taxa->get_taxa_id()][locale]), state) ==(BOOL)TRUE)
#ifdef _OPEN_MP_VERSION_
					transpoint_condprobs1_thread[myid][i] = get_tip_prob(locale, taxa, rate_num, i);
#else
					transpoint_condprobs1[i] = get_tip_prob(locale, taxa, rate_num, i);
#endif
				
			}
		}
		else {
			for(i=0; i<curr_exchange->get_condlike_size(); i++)
			{
				for(j=0; j<curr_exchange->get_condlike_size(); j++)
				{
					if (is_trans_in_transpoint_state(dupl_to_loss_state (i), dupl_to_loss_state (j), state) ==(BOOL)TRUE)
#ifdef _OPEN_MP_VERSION_
						transpoint_condprobs1_thread[myid][i]+=taxa->get_cond_prob_locale(myid, j)*taxa->get_trpb(rate_num, i, j);
#else
						transpoint_condprobs1[i]+=taxa->get_cond_prob(j)*taxa->get_trpb(rate_num, i, j);
#endif
				} 
				
			}
		}

		if (taxa->get_sibling()->is_tip() == (BOOL)TRUE) {
			for(i=0; i<curr_exchange->get_condlike_size(); i++) 
#ifdef _OPEN_MP_VERSION_
				transpoint_condprobs1_thread[myid][i] *=	get_tip_prob(locale, taxa->get_sibling(), rate_num, i);
#else
				transpoint_condprobs1[i] *=	get_tip_prob(locale, taxa->get_sibling(), rate_num, i);
#endif

		}
		else {
			for(i=0; i<curr_exchange->get_condlike_size(); i++) {
				sib_prob=0;
#ifdef _OPEN_MP_VERSION_
				for(j=0; j<curr_exchange->get_condlike_size(); j++) 
					sib_prob+=taxa->get_sibling()->get_trpb(rate_num, i, j)*taxa->get_sibling()->get_cond_prob_locale(myid, j);

				transpoint_condprobs1_thread[myid][i]*=sib_prob;
#else
				for(j=0; j<curr_exchange->get_condlike_size(); j++) 
					sib_prob+=taxa->get_sibling()->get_trpb(rate_num, i, j)*taxa->get_sibling()->get_cond_prob(j);
				transpoint_condprobs1[i]*=sib_prob;
#endif
			}
		}

		//Now walk down the tree with these partial conditional probs
		curr=curr->get_parent();

		while(curr_tree->is_root(curr) == (BOOL)FALSE) {
			for(i=0; i<curr_exchange->get_condlike_size(); i++) {
				new_prob=0;
				
				for(j=0; j<curr_exchange->get_condlike_size(); j++) 
#ifdef _OPEN_MP_VERSION_
					new_prob += transpoint_condprobs1_thread[myid][j]*curr->get_trpb(rate_num, i, j);
#else
						new_prob += transpoint_condprobs1[j]*curr->get_trpb(rate_num, i, j);
#endif
			
				if (curr->get_sibling()->is_tip() == (BOOL)TRUE) 
							sib_prob = get_tip_prob(locale, curr->get_sibling(), rate_num, i);
				else {
					sib_prob=0;
					for(j=0; j<curr_exchange->get_condlike_size(); j++) 
#ifdef _OPEN_MP_VERSION_
						sib_prob+=curr->get_sibling()->get_trpb(rate_num, i, j)*curr->get_sibling()->get_cond_prob_locale(myid, j);
#else
						sib_prob+=curr->get_sibling()->get_trpb(rate_num, i, j)*curr->get_sibling()->get_cond_prob(j);
#endif
				
				}
#ifdef _OPEN_MP_VERSION_
				transpoint_condprobs2_thread[myid][i]=new_prob*sib_prob;
#else
				transpoint_condprobs2[i]=new_prob*sib_prob;
#endif
			}

			for(i=0; i<curr_exchange->get_condlike_size(); i++)
#ifdef _OPEN_MP_VERSION_
				transpoint_condprobs1_thread[myid][i]=transpoint_condprobs2_thread[myid][i];
#else
			transpoint_condprobs1[i]=transpoint_condprobs2[i];
#endif

			curr=curr->get_parent();
		}
	}
	else {
		for(i=0; i<curr_exchange->get_condlike_size(); i++)
#ifdef _OPEN_MP_VERSION_
			transpoint_condprobs1_thread[myid][i]=curr_tree->find_root()->get_cond_prob_locale(myid, i);
#else
			transpoint_condprobs1[i]=curr_tree->find_root()->get_cond_prob(i);
#endif
	}


	retval=0;
	for(i=0; i<curr_exchange->get_condlike_size(); i++) {
		if ((curr_tree->is_root(taxa) == (BOOL)FALSE) || 
			(is_trans_in_transpoint_state(BOTH_PRESENT, dupl_to_loss_state(i), state)==(BOOL)TRUE)) 
#ifdef _OPEN_MP_VERSION_
			retval += curr_tree->find_root()->get_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), i) *
			transpoint_condprobs1_thread[myid][i];
#else
				retval += curr_tree->find_root()->get_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), i) *
					transpoint_condprobs1[i];
#endif

	}

	return(retval);
}



double Dupl_NoState_Base_model::find_appropriate_ln_like()
{
	int i, j, k, num_prob_indices;
#ifdef MAC_TIME
	uint64_t  pretime, pptime, walktime, posttime, preptime, prepaftertime, siteaftertime, makelnltime, makelnlaftertime, sitetotal=0;
#endif
#ifdef DO_TIME
	unsigned long pretime, pptime, walktime, posttime, preptime, prepaftertime, siteaftertime, makelnltime, makelnlaftertime, sitetotal=0;
#endif
	double lnL, sum_prob, site_prob;


	
#ifdef MAC_TIME
	pretime= mach_absolute_time();
#endif
	
#ifdef DO_TIME
	pretime= RunTime();
#endif
	
	if(the_tracks->tracking_current() ==(BOOL)FALSE)
		the_tracks->update_tracking();
	
	calc_pattern_probs();

#ifdef MAC_TIME
	pptime= mach_absolute_time();
	cout<<"Pattern prob calc time: "<<pptime-pretime<<endl;
#endif	
	
#ifdef DO_TIME
	pptime= RunTime();
	cout<<"Pattern prob calc time: "<<pptime-pretime<<endl;
#endif	
	
	num_prob_indices=get_locus_state_index(0);
	site_mask(0);

#ifdef _OPEN_MP_VERSION_
#pragma omp parallel for private (i, j, site_prob)
#endif
	for(i=0; i<pow2[curr_exchange->get_num_taxa()]; i++)  {
		for(j=0; j<num_prob_indices; j++)
#ifdef _OPEN_MP_VERSION_
			sum_probs_thread[omp_get_thread_num()][j] =  pattern_probs[dupl_indices[j]][site_indices[j]][i^transition_mask];
		site_prob = get_ln_sum_prob_sort(sum_probs_thread[omp_get_thread_num()], num_prob_indices);
#else
			sum_probs[j] =  pattern_probs[dupl_indices[j]][site_indices[j]][i^transition_mask];
		site_prob = get_ln_sum_prob_sort(sum_probs, num_prob_indices);
#endif
		
		
		
		last_cumulative_probs[i] = site_prob;
	}

	//cout<<"Site 0 Check prob "<<last_cumulative_probs[0]<<" and 32: "<<last_cumulative_probs[31]<<endl;

	for(i=1; i<the_homologs->get_num_homologs(); i++) {
		num_prob_indices=get_locus_state_index(i);
		//pair_site_mask(i, i-1);
		
		
		site_mask(i);
		build_transition_vector(i);
#ifdef MAC_TIME
		prepaftertime= mach_absolute_time();
		//cout<<"Time for prep: "<<prepaftertime-preptime<<endl;
#endif
		
#ifdef DO_TIME
		prepaftertime= RunTime();
		//cout<<"Time for prep: "<<prepaftertime-preptime<<endl;
#endif

#ifdef _OPEN_MP_VERSION_
#pragma omp parallel for private (j, k, sum_prob, site_prob)
#endif
		for(j=0; j<pow2[curr_exchange->get_num_taxa()]; j++) {
			for(k=0; k<pow2[curr_exchange->get_num_taxa()]; k++) 
#ifdef _OPEN_MP_VERSION_
				sum_probs_thread[omp_get_thread_num()][k]=state_trans_probs[j ^ k]+last_cumulative_probs[k];
			sum_prob=get_ln_sum_prob_sort(sum_probs_thread[omp_get_thread_num()]);
#else
				sum_probs[k]=state_trans_probs[j ^ k]+last_cumulative_probs[k];
			sum_prob=get_ln_sum_prob_sort(sum_probs);
#endif		
			
			
			for(k=0; k<num_prob_indices; k++)
#ifdef _OPEN_MP_VERSION_
				sum_probs_thread[omp_get_thread_num()][k]= pattern_probs[dupl_indices[k]][site_indices[k]][j^transition_mask];
			site_prob=get_ln_sum_prob_sort(sum_probs_thread[omp_get_thread_num()], num_prob_indices);
#else
				sum_probs[k]= pattern_probs[dupl_indices[k]][site_indices[k]][j^transition_mask];
			site_prob=get_ln_sum_prob_sort(sum_probs, num_prob_indices);
			
			
#endif		
			
			cumulative_probs[j] = site_prob + sum_prob;
		}
		//cout<<"Site "<<i<<"Check prob "<<cumulative_probs[0]<<endl;
#ifdef MAC_TIME
		siteaftertime= mach_absolute_time();
		sitetotal+=siteaftertime-prepaftertime;
		//cout<<"Time for one site calc: "<<i<<": "<<siteaftertime-prepaftertime<<endl;
#endif
		
#ifdef DO_TIME
		siteaftertime= RunTime();
		sitetotal+=siteaftertime-prepaftertime;
		//cout<<"Time for one site calc: "<<i<<": "<<siteaftertime-prepaftertime<<endl;
#endif
		
#ifdef _OPEN_MP_VERSION_
#pragma omp parallel for private (j)
#endif
		for(j=0; j<pow2[curr_exchange->get_num_taxa()]; j++) 
			last_cumulative_probs[j]=cumulative_probs[j];

	}
	

		lnL=get_ln_sum_prob(cumulative_probs);
	

#ifdef MAC_TIME
	posttime=mach_absolute_time();
	cout<<"Runtime for likelihood calc: "<<posttime-pretime<<endl;
	cout<<"Time in main loop: "<<sitetotal<<endl;
#endif	
	
#ifdef DO_TIME
	posttime=RunTime();
	cout<<"Runtime for likelihood calc: "<<posttime-pretime<<endl;
	cout<<"Time in main loop: "<<sitetotal<<endl;
#endif	
	
	return(lnL);
}

double Dupl_NoState_Base_model::get_post_prob(int site, int pattern)
{
	if ((pattern >=0) &&(pattern < pow2[curr_exchange->get_num_taxa()])) 
		return(post_probs[site][pattern]);
	else
		return(post_probs[site][0]);

}


double Dupl_NoState_Base_model::get_site_prob(int site, DUPL_LOSS_STATES state)
{
	if (site_probs !=0) {
		return(site_probs[site][loss_state_to_dupl(state)]);
	}
	else 
		return(0.0);
}


void Dupl_NoState_Base_model::print_tracking_probs(char *outfile)
{
	int i, j, k;
	ofstream fout;

	get_gene_conditional_probs();

	fout.open(outfile);
	if (!fout.fail()) {
		for(i=0; i<curr_exchange->get_num_taxa(); i++)
			fout<<(*the_tracks)[i].get_genome()->get_name()<<"\t";

		for(i=0; i<curr_exchange->get_num_taxa(); i++)
			fout<<(*the_tracks)[i].get_genome()->get_name()<<"\t";


		for(j=0; j<pow2[curr_exchange->get_num_taxa()]; j++) {
			fout<<"#";
			for(k=curr_exchange->get_num_taxa()-1; k>=0; k--) {
				if (pow2[k] & j)
					fout<<"1";
				else
					fout<<"0";
			}
			if (j!= pow2[curr_exchange->get_num_taxa()]-1)
				fout<<"\t";
			else
				fout<<"\n";
		}

		for (i=0; i<the_homologs->get_num_homologs(); i++) {
			for(j=0; j<curr_exchange->get_num_taxa(); j++) {
				if ((*the_tracks)[j].get_gene_track(i, 0)->my_locus != 0)
						fout<<(*the_genomes)[j][(*the_tracks)[j].get_gene_track(i, 0)->
						my_locus->get_contig((*the_tracks)[j].get_gene_track(i, 0)->index_num)]
					[(*the_tracks)[j].get_gene_track(i, 0)->
						my_locus->get_gene((*the_tracks)[j].get_gene_track(i, 0)->index_num)].get_name()<<"\t";
				else
					fout<<"NONE\t";
			}
			for(j=0; j<curr_exchange->get_num_taxa(); j++) {
				if ((*the_tracks)[j].get_gene_track(i, 1)->my_locus != 0)
					fout<<(*the_genomes)[j][(*the_tracks)[j].get_gene_track(i, 1)->
					my_locus->get_contig((*the_tracks)[j].get_gene_track(i, 1)->index_num)]
					[(*the_tracks)[j].get_gene_track(i, 1)->
					my_locus->get_gene((*the_tracks)[j].get_gene_track(i, 1)->index_num)].get_name()<<"\t";
				else
					fout<<"NONE\t";
			}
			for(j=0; j<pow2[curr_exchange->get_num_taxa()]; j++) {
				if (j != pow2[curr_exchange->get_num_taxa()]-1) 
					fout<<post_probs[i][j]<<"\t";
				else
					fout<<post_probs[i][j]<<"\n";
			}
		}
	}	
}
	
void Dupl_NoState_Base_model::get_gene_conditional_probs()
//Calculates the probabilities of various trackings for each "pillar" given those to the left and right
{
	int i, j, k, num_prob_indices, last_mask=0, this_mask;
	double sum, site_prob, sum_prob, *temp_probs;
	
	post_probs=new double*[the_homologs->get_num_homologs()];
	left_cond_probs=new double*[the_homologs->get_num_homologs()];
	right_cond_probs=new double*[the_homologs->get_num_homologs()];
	
	temp_probs=new double[pow2[curr_exchange->get_num_taxa()]];
	
	for(i=0; i<the_homologs->get_num_homologs(); i++) {
		post_probs[i]=new double[pow2[curr_exchange->get_num_taxa()]];
		left_cond_probs[i]=new double[pow2[curr_exchange->get_num_taxa()]];
		right_cond_probs[i]=new double[pow2[curr_exchange->get_num_taxa()]];
	}
	

	the_tracks->update_tracking();
	recalculate_transprobs();
	calc_pattern_probs();

	num_prob_indices=get_locus_state_index(0);
	site_mask(0);

	for(i=0; i<pow2[curr_exchange->get_num_taxa()]; i++)  {
		for(j=0; j<num_prob_indices; j++)
			sum_probs[j] = pattern_probs[dupl_indices[j]][site_indices[j]][i^transition_mask];
		site_prob = get_ln_sum_prob_sort(sum_probs, num_prob_indices);
		left_cond_probs[0][i] = site_prob;
	}

	num_prob_indices=get_locus_state_index(the_homologs->get_num_homologs()-1);
	site_mask(the_homologs->get_num_homologs()-1);

	for(i=0; i<pow2[curr_exchange->get_num_taxa()]; i++) {
		for(j=0; j<num_prob_indices; j++)
			sum_probs[j] = pattern_probs[dupl_indices[j]][site_indices[j]][i^transition_mask];
		site_prob = get_ln_sum_prob_sort(sum_probs, num_prob_indices);
		right_cond_probs[the_homologs->get_num_homologs()-1][i] = site_prob; 
		;
	}

	for(i=1; i<the_homologs->get_num_homologs(); i++) {
		num_prob_indices=get_locus_state_index(i);
		site_mask(i);
		//pair_site_mask(i, i-1);
		build_transition_vector(i, (BOOL)TRUE);
#ifdef _OPEN_MP_VERSION_
#pragma omp parallel for private (j, k, sum_prob, site_prob)
#endif	
		for(j=0; j<pow2[curr_exchange->get_num_taxa()]; j++) {
		
			left_cond_probs[i][j]=0;
			for(k=0; k<pow2[curr_exchange->get_num_taxa()]; k++) {
#ifdef _OPEN_MP_VERSION_
				sum_probs_thread[omp_get_thread_num()][k]=state_trans_probs[j ^ k]+left_cond_probs[i-1][k];
#else
				sum_probs[k]=state_trans_probs[j ^ k]+left_cond_probs[i-1][k];
#endif
			}
#ifdef _OPEN_MP_VERSION_
			sum_prob=get_ln_sum_prob_sort(sum_probs_thread[omp_get_thread_num()]);
#else
			sum_prob=get_ln_sum_prob_sort(sum_probs);
#endif

			for(k=0; k<num_prob_indices; k++)
#ifdef _OPEN_MP_VERSION_
				sum_probs_thread[omp_get_thread_num()][k] = pattern_probs[dupl_indices[k]][site_indices[k]][j^transition_mask];
#else
				sum_probs[k] = pattern_probs[dupl_indices[k]][site_indices[k]][j^transition_mask];
#endif
#ifdef _OPEN_MP_VERSION_
			site_prob = get_ln_sum_prob_sort(sum_probs_thread[omp_get_thread_num()], num_prob_indices);
#else
			site_prob = get_ln_sum_prob_sort(sum_probs, num_prob_indices);
#endif
			left_cond_probs[i][j] = site_prob+sum_prob;
		}

	}

	for(i=the_homologs->get_num_homologs()-2; i>0; i--) {
	num_prob_indices=get_locus_state_index(i);
		//pair_site_mask(i, i+1);
		site_mask(i);
		build_transition_vector(i, (BOOL)FALSE);
#ifdef _OPEN_MP_VERSION_
#pragma omp parallel for private (j, k, sum_prob, site_prob)
#endif			
		for(j=0; j<pow2[curr_exchange->get_num_taxa()]; j++) {
			right_cond_probs[i][j]=0;
			for(k=0; k<pow2[curr_exchange->get_num_taxa()]; k++) {
#ifdef _OPEN_MP_VERSION_
				sum_probs_thread[omp_get_thread_num()][k]=state_trans_probs[j ^ k]+right_cond_probs[i+1][k];
#else
				sum_probs[k]=state_trans_probs[j ^ k]+right_cond_probs[i+1][k];
#endif
			}
#ifdef _OPEN_MP_VERSION_
			sum_prob = get_ln_sum_prob_sort(sum_probs_thread[omp_get_thread_num()]);
#else
			sum_prob = get_ln_sum_prob_sort(sum_probs);
#endif

			for(k=0; k<num_prob_indices; k++) 
#ifdef _OPEN_MP_VERSION_
				sum_probs_thread[omp_get_thread_num()][k] = pattern_probs[dupl_indices[k]][site_indices[k]][j^transition_mask];
#else
				sum_probs[k] = pattern_probs[dupl_indices[k]][site_indices[k]][j^transition_mask];	
#endif

#ifdef _OPEN_MP_VERSION_
			site_prob=get_ln_sum_prob_sort(sum_probs_thread[omp_get_thread_num()], num_prob_indices);
#else
			site_prob=get_ln_sum_prob_sort(sum_probs, num_prob_indices);
#endif

			right_cond_probs[i][j] = site_prob + sum_prob;
		}

	}
	
	
	num_prob_indices=get_locus_state_index(0);
	//pair_site_mask(0, 1);
	site_mask(0);
	build_transition_vector(0,(BOOL)FALSE);
	
	for(j=0; j<pow2[curr_exchange->get_num_taxa()]; j++) {
		
		post_probs[0][j]=0;
		for(k=0; k<pow2[curr_exchange->get_num_taxa()]; k++) {
			sum_probs[k]=state_trans_probs[j ^ k]+right_cond_probs[1][k];
		}
		sum_prob = get_ln_sum_prob_sort(sum_probs);

		for(k=0; k<num_prob_indices; k++)
			sum_probs[k] = pattern_probs[dupl_indices[k]][site_indices[k]][j^transition_mask];
		site_prob = get_ln_sum_prob_sort(sum_probs, num_prob_indices);
		post_probs[0][j] = site_prob + sum_prob;
	}

	num_prob_indices=get_locus_state_index(the_homologs->get_num_homologs()-1);
	//pair_site_mask(the_homologs->get_num_homologs()-1, the_homologs->get_num_homologs()-2);
	site_mask(the_homologs->get_num_homologs()-1);
	build_transition_vector(the_homologs->get_num_homologs()-1, (BOOL)TRUE);
	
	for(j=0; j<pow2[curr_exchange->get_num_taxa()]; j++) {
		
		post_probs[the_homologs->get_num_homologs()-1][j]=0;
		for(k=0; k<pow2[curr_exchange->get_num_taxa()]; k++) {
			sum_probs[k]=state_trans_probs[j ^ k]+left_cond_probs[the_homologs->get_num_homologs()-2][k];
		}
		sum_prob = get_ln_sum_prob_sort(sum_probs);

		for(k=0; k<num_prob_indices; k++)
			sum_probs[k] = pattern_probs[dupl_indices[k]][site_indices[k]][j^transition_mask];
		site_prob = get_ln_sum_prob_sort(sum_probs, num_prob_indices);
		post_probs[the_homologs->get_num_homologs()-1][j] = site_prob + sum_prob; 
	}

	for(i=1; i<the_homologs->get_num_homologs()-1; i++) {
		num_prob_indices=get_locus_state_index(i);
		//pair_site_mask(i, i-1);
		site_mask(i);
		build_transition_vector(i, (BOOL)TRUE);
#ifdef _OPEN_MP_VERSION_
#pragma omp parallel for private (j, k, sum_prob, site_prob)
#endif				
		for(j=0; j<pow2[curr_exchange->get_num_taxa()]; j++) {
		
			post_probs[i][j]=0;
			for(k=0; k<pow2[curr_exchange->get_num_taxa()]; k++) {
#ifdef _OPEN_MP_VERSION_
				sum_probs_thread[omp_get_thread_num()][k]=state_trans_probs[j ^ k]+left_cond_probs[i-1][k];
#else
				sum_probs[k]=state_trans_probs[j ^ k]+left_cond_probs[i-1][k];
#endif
			}
#ifdef _OPEN_MP_VERSION_
			sum_prob = get_ln_sum_prob_sort(sum_probs_thread[omp_get_thread_num()]);
#else
			sum_prob = get_ln_sum_prob_sort(sum_probs);
#endif

			for(k=0; k<num_prob_indices; k++)
#ifdef _OPEN_MP_VERSION_
				sum_probs_thread[omp_get_thread_num()][k] = pattern_probs[dupl_indices[k]][site_indices[k]][j^transition_mask];
#else
				sum_probs[k] = pattern_probs[dupl_indices[k]][site_indices[k]][j^transition_mask];
#endif
#ifdef _OPEN_MP_VERSION_
			site_prob = get_ln_sum_prob_sort(sum_probs_thread[omp_get_thread_num()], num_prob_indices);
#else
			site_prob = get_ln_sum_prob_sort(sum_probs, num_prob_indices);
#endif
			
			post_probs[i][j] = site_prob + sum_prob;
		}

	
		//pair_site_mask(i, i+1);
		site_mask(i);
		build_transition_vector(i, (BOOL)FALSE);	
		for(j=0; j<pow2[curr_exchange->get_num_taxa()]; j++) {	
			for(k=0; k<pow2[curr_exchange->get_num_taxa()]; k++) {
				sum_probs_r[k]=state_trans_probs[j ^ k]+right_cond_probs[i+1][k];
			}
			post_probs[i][j] += get_ln_sum_prob_sort(sum_probs_r);
	
		}

		//pair_site_mask(i, i-1);
		/*site_mask(i);
		this_mask = (last_mask ^ transition_mask);
		if(this_mask != 0) {
			for(j=0; j<pow2[curr_exchange->get_num_taxa()]; j++) {
				temp_probs[(j ^ this_mask)]=post_probs[i][j];
			}
			for(j=0; j<pow2[curr_exchange->get_num_taxa()]; j++) {
				post_probs[i][j]=temp_probs[j];
			}
		}
		last_mask=this_mask;*/

	}

	//pair_site_mask(the_homologs->get_num_homologs()-1, the_homologs->get_num_homologs()-2);
	/*site_mask(the_homologs->get_num_homologs()-1);
	this_mask = (last_mask ^ transition_mask);
	if(this_mask != 0) {
		for(j=0; j<pow2[curr_exchange->get_num_taxa()]; j++) {
			temp_probs[(j ^ this_mask)]=post_probs[the_homologs->get_num_homologs()-1][j];
		}
			for(j=0; j<pow2[curr_exchange->get_num_taxa()]; j++) {
				post_probs[the_homologs->get_num_homologs()-1][j]=temp_probs[j];
			}
	}*/
	
	//Normalize the probs:
	for(i=0; i<the_homologs->get_num_homologs(); i++) {
		sum=get_ln_sum_prob(post_probs[i]);
		for(j=0; j<pow2[curr_exchange->get_num_taxa()]; j++)
			post_probs[i][j]=exp(post_probs[i][j]-sum);
	}

	
	delete[] temp_probs;

}


int Dupl_NoState_Base_model::get_locus_num_dupls(int locus)
{
	int i, num_dupls=0;

	for(i=0; i<curr_exchange->get_num_taxa(); i++) {
		if ((*the_tracks)[i].get_gene_track(locus, 0)->my_locus !=0)
			if ((*the_tracks)[i].get_gene_track(locus, 0)->my_locus->has_duplicate() == (BOOL)TRUE)
				num_dupls++;
	}
	return(num_dupls);
}



int Dupl_NoState_Base_model::get_locus_state_index(int locus)
{
	int i, j, num_ambig, num_dupls = 0, shift;

	if (curr_exchange->use_track_gaps_as_missing() == (BOOL)FALSE) {
		for(i=0; i<curr_exchange->get_num_taxa(); i++) {
			dupl_states[i]=1;
			if ((*the_tracks)[i].get_gene_track(locus, 0)->my_locus !=0) {
				if ((*the_tracks)[i].get_gene_track(locus, 0)->my_locus->has_duplicate() == (BOOL)TRUE) {
					dupl_states[i]=2;
					num_dupls++;
				}
			}
		}
			site_indices[0]=get_dupl_state_index();
			dupl_indices[0]= num_dupls;
			return(1);
	}
	else {
		num_ambig = 0;
		shift=0;
		for(i=0; i<curr_exchange->get_num_taxa(); i++) {
			is_ambig[i]=(BOOL)FALSE;
			if ((*the_tracks)[i].get_gene_track(locus, 0)->my_locus !=0) {
				if ((*the_tracks)[i].get_gene_track(locus, 0)->my_locus->has_duplicate() == (BOOL)FALSE) {
					//Could have a gap on track 1
					if (((*the_tracks)[i].get_gene_track(locus, 1)->last == 0) && 
						((*the_tracks)[i].get_gene_track(locus, 1)->next == 0)) {
						num_ambig++;
						is_ambig[i]=(BOOL)TRUE;
					}
					
				}
			}
			else {
				//Could have a gap on track 0
					if (((*the_tracks)[i].get_gene_track(locus, 0)->last == 0) && 
						((*the_tracks)[i].get_gene_track(locus, 0)->next == 0)) {
						num_ambig++;
						is_ambig[i]=(BOOL)TRUE;
					}
			}
		}


		for(i=0; i< pow2[num_ambig]; i++) {
		//loop over all possible assignments of ambiguous duplicates
			num_dupls=0;
			for(j=0; j<curr_exchange->get_num_taxa(); j++) {
				if (is_ambig[j] == (BOOL)FALSE) {
					shift++;
					dupl_states[j]=1;
					if ((*the_tracks)[j].get_gene_track(locus, 0)->my_locus !=0) {
						if ((*the_tracks)[j].get_gene_track(locus, 0)->my_locus->has_duplicate() == (BOOL)TRUE) {
							dupl_states[j]=2;
							num_dupls++;
						}
					}
				}
				else {
					if ((pow2[j]>>shift & i) != 0) {
						dupl_states[j] = 2;
						num_dupls++;
					}
					else 
						dupl_states[j] = 1;
				}
			}

			site_indices[i]=get_dupl_state_index();
			dupl_indices[i]= num_dupls;


		}

		return(pow2[num_ambig]);

	}
}



int Dupl_NoState_Base_model::get_dupl_state_index()
{
	int i, num_dupl=0, retval=0;

	for(i=0; i<curr_exchange->get_num_taxa(); i++)
		if (dupl_states[i] == 2)
			num_dupl++;

		if ((num_dupl==0) || (num_dupl==curr_exchange->get_num_taxa()))
			retval=0;
		else 
			recurse_states(num_dupl, retval, 0);
			
		return(retval);

}



void Dupl_NoState_Base_model::recurse_states(int num_dupls, int &so_far, int pos)
{
	int i, num_ahead, total;

	i=pos;
	while(dupl_states[i] != 2) i++;
	
	//This is the number of patterns with num_dupls duplications appearing *before* all those
	num_ahead=n_choose_k((curr_exchange->get_num_taxa()-pos)-(i-pos), num_dupls);
	total=n_choose_k((curr_exchange->get_num_taxa()-pos), num_dupls);
	so_far+=total-num_ahead;
	
	if (num_dupls > 1)
		recurse_states(num_dupls-1, so_far, i+1);

}



void Dupl_NoState_Base_model::get_states(int index_num, int num_dupls, int pos)
{
	int i, num_ahead, total, num_pos;
	
	if (num_dupls == 0)
	{
		for (i=pos; i<curr_exchange->get_num_taxa(); i++)
			dupl_states[i]=1;
	}
	else if (num_dupls == curr_exchange->get_num_taxa())
	{
		for (i=pos; i<curr_exchange->get_num_taxa(); i++)
			dupl_states[i]=2;
	}
	else {
		num_ahead=n_choose_k((curr_exchange->get_num_taxa()-1-pos), num_dupls);
		total=n_choose_k((curr_exchange->get_num_taxa()-pos), num_dupls);
		num_pos=total-num_ahead;

		if ((index_num < num_pos) || (num_pos == 0)) {
			dupl_states[pos]=2;
			num_dupls--;
		}
		else {

			dupl_states[pos]=1;
			index_num-=num_pos;
		}

		if (pos < curr_exchange->get_num_taxa()-1)
			get_states(index_num, num_dupls, pos+1);
	}


}



void Dupl_NoState_Base_model::allocate_state_model (Exchange *cexchange, Tree *ctree, Clade *cgenomes, WGD_Data *chomologs)
{
	int i, j;

	if (cexchange->get_strand_switch_prob() == 0.0)
		cexchange->set_strand_switch_prob(0.01);
	curr_exchange=cexchange;
#ifdef _OPEN_MP_VERSION_
	curr_data = new Sequence_dataset(cexchange->get_num_taxa(), cexchange->get_num_open_mp_threads(), DUPL_STATUS);
#else
	curr_data = new Sequence_dataset(cexchange->get_num_taxa(), 1, DUPL_STATUS);
#endif
	
	//cout<<"Curr data has size "<<(*curr_data)[0].Sequence_size()<<" sites"<<" : Threads by exchange: "<<cexchange->get_num_open_mp_threads()<<endl;
	
	the_genomes=cgenomes;
	the_homologs=chomologs;
	the_tracks = new WGD_Tracks(chomologs, cgenomes);
	site_probs=0;
	do_transpoint_probs=(BOOL)FALSE;

	pow2[0]=1;
	for(i=1; i<32; i++)
		pow2[i] =pow2[i-1]*2;
	
	dupl_states=new int[cexchange->get_num_taxa()];



	state_trans_probs=new double [pow2[cexchange->get_num_taxa()]];
	cumulative_probs=new double [pow2[cexchange->get_num_taxa()]];
	last_cumulative_probs=new double [pow2[cexchange->get_num_taxa()]];
	sum_probs=new double [pow2[cexchange->get_num_taxa()]];
	
#ifdef _OPEN_MP_VERSION_
	sum_probs_thread=new double * [cexchange->get_num_open_mp_threads()];
	for (i=0; i<cexchange->get_num_open_mp_threads(); i++) 
		sum_probs_thread[i] = new double [pow2[cexchange->get_num_taxa()]];
#endif
	
	sum_probs_r=new double [pow2[cexchange->get_num_taxa()]];
	sum_probs_used=new BOOL [pow2[cexchange->get_num_taxa()]];

#ifdef _OPEN_MP_VERSION_
	sum_probs_used_thread=new BOOL * [cexchange->get_num_open_mp_threads()];
	for (i=0; i<cexchange->get_num_open_mp_threads(); i++) 
		sum_probs_used_thread[i]=new BOOL [pow2[cexchange->get_num_taxa()]];
#endif
	
	site_indices = new int [pow2[cexchange->get_num_taxa()]];
	dupl_indices = new int [pow2[cexchange->get_num_taxa()]];
	is_ambig = new BOOL [cexchange->get_num_taxa()];

	post_probs=0;
	branch_state_probs=0;
	
	
	//First index is number of duplicates 0->n
	pattern_probs=new double **[cexchange->get_num_taxa()+1];

	for(i=0; i<=cexchange->get_num_taxa(); i++) {
	//Next index is number of arrangements of i duplicates in n taxa
		pattern_probs[i] = new double * [n_choose_k(cexchange->get_num_taxa(), i)];
		for(j=0; j<n_choose_k(cexchange->get_num_taxa(), i); j++)
			//Final index is 2^n
			//Note that we're assuming we can fit 2^n states array into a 32-bit pointer
			pattern_probs[i][j]=new double [(int)pow2[cexchange->get_num_taxa()]];
	}

	
	build_masks();
	assemble(cexchange, curr_data, ctree);
}


void Dupl_NoState_Base_model::will_calculate_branch_transpoint_probs()
{
	int i, j, k, l;

	do_transpoint_probs=(BOOL)TRUE;
	branch_transpoint_probs=new double ****[curr_exchange->get_num_taxa()+1];
#ifdef _OPEN_MP_VERSION_
	transpoint_condprobs1_thread = new double* [curr_exchange->get_num_open_mp_threads()];
	transpoint_condprobs2_thread = new double* [curr_exchange->get_num_open_mp_threads()];
	for(i=0; i<curr_exchange->get_num_open_mp_threads(); i++) {
		transpoint_condprobs1_thread[i] = new double [curr_exchange->get_condlike_size()];
		transpoint_condprobs2_thread[i] = new double [curr_exchange->get_condlike_size()];
	}
#endif
	
	transpoint_condprobs1 = new double[curr_exchange->get_condlike_size()];
	transpoint_condprobs2 = new double[curr_exchange->get_condlike_size()];



	for(i=0; i<=curr_exchange->get_num_taxa(); i++) {
		branch_transpoint_probs[i]=new double *** [n_choose_k(curr_exchange->get_num_taxa(), i)];
		for(j=0; j<n_choose_k(curr_exchange->get_num_taxa(), i); j++) {
				branch_transpoint_probs[i][j]=new double ** [pow2[curr_exchange->get_num_taxa()]];
			for(k=0; k<pow2[curr_exchange->get_num_taxa()]; k++) {
				branch_transpoint_probs[i][j][k] = new double *[curr_exchange->get_num_branches()];
					for(l=0; l<curr_exchange->get_num_branches(); l++)
						branch_transpoint_probs[i][j][k][l]=new double [NUM_TRANSPOINT_STATES];
			}
		}	
	}

}




void Dupl_NoState_Base_model::calc_transpoint_branch_probs()
{
	int i, j, k, l, m, num_prob_indices;
#ifdef MAC_TIME
	uint64_t  pretime, pptime, walktime, posttime, preptime, prepaftertime, siteaftertime, makelnltime, makelnlaftertime, sitetotal=0;
#endif
#ifdef DO_TIME
	unsigned long pretime, pptime, walktime, posttime, preptime, prepaftertime, siteaftertime, makelnltime, makelnlaftertime, sitetotal=0;
#endif
	
	double *transpoint_state_probs, sum_prob, site_prob, full_lnL, prop_prob;
	
	
	branch_state_probs=new double **[the_homologs->get_num_homologs()];
	for(i=0; i<the_homologs->get_num_homologs(); i++) {
		branch_state_probs[i]=new double *[curr_exchange->get_num_branches()];
		for(j=0; j<curr_exchange->get_num_branches(); j++)
			branch_state_probs[i][j]=new double[NUM_TRANSPOINT_STATES];
	}
	
	transpoint_state_probs = new double  [pow2[curr_exchange->get_num_taxa()]];
	
	will_calculate_branch_transpoint_probs();
	full_lnL=find_appropriate_ln_like();

#ifdef MAC_TIME
	pretime= mach_absolute_time();
	//cout<<"Time for prep: "<<prepaftertime-preptime<<endl;
#endif
	
#ifdef DO_TIME
	pretime= RunTime();
	//cout<<"Time for prep: "<<prepaftertime-preptime<<endl;
#endif	
	
	get_gene_conditional_probs();


#ifdef MAC_TIME
	posttime= mach_absolute_time();
	cout<<"Time for get_gene_conditional_probs: "<<posttime-pretime<<endl;
#endif
	
#ifdef DO_TIME
	posttime= RunTime();
	cout<<"Time for get_gene_conditional_probs: "<<posttime-pretime<<endl;
#endif	
	
	for(m=0; m<curr_exchange->get_num_branches(); m++) {
		for(l=0; l<NUM_TRANSPOINT_STATES; l++) {
			num_prob_indices=get_locus_state_index(0);
			//pair_site_mask(0, 1);
			site_mask(0);
			build_transition_vector(0, (BOOL)FALSE);
			for(j=0; j<pow2[curr_exchange->get_num_taxa()]; j++) {
				
				for(k=0; k<pow2[curr_exchange->get_num_taxa()]; k++) {
					sum_probs[k]=state_trans_probs[j ^ k]+right_cond_probs[1][k];
				}
				sum_prob = get_ln_sum_prob_sort(sum_probs);
				
				for(k=0; k<num_prob_indices; k++) {
					if (branch_transpoint_probs[dupl_indices[k]][site_indices[k]][j^transition_mask][m][l] > 0)
						sum_probs[k] = log(branch_transpoint_probs[dupl_indices[k]][site_indices[k]][j^transition_mask][m][l]);
					else
						sum_probs[k] = 1.0;
				}
				site_prob = get_ln_sum_prob_zero_sort(sum_probs, num_prob_indices);
				
				
				if (site_prob != 1.0)
					transpoint_state_probs[j] = site_prob + sum_prob;    
				else
					transpoint_state_probs[j] = 1.0;
			}
			
			site_prob=get_ln_sum_prob_zero(transpoint_state_probs, pow2[curr_exchange->get_num_taxa()]);
			if (site_prob != 1.0) 
				branch_state_probs[0][m][l] = exp(site_prob - full_lnL);
			else 
				branch_state_probs[0][m][l]=0.0;
		}
	}
#ifdef MAC_TIME
	pretime= mach_absolute_time();
	//cout<<"Time for prep: "<<prepaftertime-preptime<<endl;
#endif
	
#ifdef DO_TIME
	pretime= RunTime();
	//cout<<"Time for prep: "<<prepaftertime-preptime<<endl;
#endif	
	
	
	for(i=1; i<the_homologs->get_num_homologs()-1; i++) {
		for(m=0; m<curr_exchange->get_num_branches(); m++) {
			for(l=0; l<NUM_TRANSPOINT_STATES; l++) {
				num_prob_indices=get_locus_state_index(i);
				site_mask(i);
				build_transition_vector(i, (BOOL)TRUE);
#ifdef _OPEN_MP_VERSION_
#pragma omp parallel for private (j, k, sum_prob, site_prob)
#endif
				for(j=0; j<pow2[curr_exchange->get_num_taxa()]; j++) {
					for(k=0; k<pow2[curr_exchange->get_num_taxa()]; k++) {
#ifdef _OPEN_MP_VERSION_
						sum_probs_thread[omp_get_thread_num()][k]=state_trans_probs[j ^ k]+left_cond_probs[i-1][k];
#else
						sum_probs[k]=state_trans_probs[j ^ k]+left_cond_probs[i-1][k];
#endif
					}
#ifdef _OPEN_MP_VERSION_
					sum_prob = get_ln_sum_prob_sort(sum_probs_thread[omp_get_thread_num()]);
#else
					sum_prob = get_ln_sum_prob_sort(sum_probs);
#endif
					
					for(k=0; k<num_prob_indices; k++) {
						if (branch_transpoint_probs[dupl_indices[k]][site_indices[k]][j^transition_mask][m][l] > 0)
#ifdef _OPEN_MP_VERSION_
							sum_probs_thread[omp_get_thread_num()][k]=log(branch_transpoint_probs[dupl_indices[k]][site_indices[k]][j^transition_mask][m][l]);
#else
							sum_probs[k] = log(branch_transpoint_probs[dupl_indices[k]][site_indices[k]][j^transition_mask][m][l]);
#endif
						else
#ifdef _OPEN_MP_VERSION_
							sum_probs_thread[omp_get_thread_num()][k]= 1.0;
#else		
							sum_probs[k] =1.0;
#endif
					}

#ifdef _OPEN_MP_VERSION_
					site_prob = get_ln_sum_prob_zero_sort(sum_probs_thread[omp_get_thread_num()], num_prob_indices);
#else
					site_prob = get_ln_sum_prob_zero_sort(sum_probs, num_prob_indices);
#endif
					if (site_prob != 1.0)
						transpoint_state_probs[j] = site_prob + sum_prob;
					else
						transpoint_state_probs[j] = 1.0;
				}
				
				
				//pair_site_mask(i, i+1);
				site_mask(i);
				build_transition_vector(i, (BOOL)FALSE);	
#ifdef _OPEN_MP_VERSION_
#pragma omp parallel for private (j, k, site_prob)
#endif
				for(j=0; j<pow2[curr_exchange->get_num_taxa()]; j++) {	
					for(k=0; k<pow2[curr_exchange->get_num_taxa()]; k++) {
#ifdef _OPEN_MP_VERSION_
						sum_probs_thread[omp_get_thread_num()][k]=state_trans_probs[j ^ k]+right_cond_probs[i+1][k];
#else
						sum_probs_r[k]=state_trans_probs[j ^ k]+right_cond_probs[i+1][k];
#endif
					}
					if (transpoint_state_probs[j] != 1.0)
#ifdef _OPEN_MP_VERSION_
						transpoint_state_probs[j] += get_ln_sum_prob_sort(sum_probs_thread[omp_get_thread_num()]);
#else
						transpoint_state_probs[j] += get_ln_sum_prob_sort(sum_probs_r);
#endif
					
				}
				
				site_prob=get_ln_sum_prob_zero(transpoint_state_probs, pow2[curr_exchange->get_num_taxa()]);
				if (site_prob != 1.0) 
					branch_state_probs[i][m][l] = exp(site_prob - full_lnL);
				else
					branch_state_probs[i][m][l]=0.0;
			}
		}
	}
	
#ifdef MAC_TIME
	posttime= mach_absolute_time();
	cout<<"Time for main homolog loop in calc_transpoint_probs: "<<posttime-pretime<<endl;
#endif
	
#ifdef DO_TIME
	posttime= RunTime();
	cout<<"Time for main homolog loop in calc_transpoint_probs: "<<posttime-pretime<<endl;
#endif		
	
	num_prob_indices=get_locus_state_index(the_homologs->get_num_homologs()-1);
	site_mask(the_homologs->get_num_homologs()-1);
	build_transition_vector(the_homologs->get_num_homologs()-1, (BOOL)TRUE);
	for(m=0; m<curr_exchange->get_num_branches(); m++) {
		for(l=0; l<NUM_TRANSPOINT_STATES; l++) {
			for(j=0; j<pow2[curr_exchange->get_num_taxa()]; j++) {
				
				post_probs[the_homologs->get_num_homologs()-1][j]=0;
				for(k=0; k<pow2[curr_exchange->get_num_taxa()]; k++) {
					sum_probs[k]=state_trans_probs[j ^ k]+left_cond_probs[the_homologs->get_num_homologs()-2][k];
				}
				sum_prob = get_ln_sum_prob_sort(sum_probs);
				
				for(k=0; k<num_prob_indices; k++) {
					if (branch_transpoint_probs[dupl_indices[k]][site_indices[k]][j^transition_mask][m][l] > 0)
						sum_probs[k] = log(branch_transpoint_probs[dupl_indices[k]][site_indices[k]][j^transition_mask][m][l]);
					else
						sum_probs[k] = 1.0;
				}
				
				site_prob = get_ln_sum_prob_zero_sort(sum_probs, num_prob_indices);
				if (site_prob != 1.0)
					transpoint_state_probs[j] = site_prob + sum_prob; 
				else
					transpoint_state_probs[j] = 1.0;
			}
			
			site_prob=get_ln_sum_prob_zero(transpoint_state_probs, pow2[curr_exchange->get_num_taxa()]);
			if (site_prob != 1.0) 
				branch_state_probs[the_homologs->get_num_homologs()-1][m][l]= exp(site_prob - full_lnL);
			else
				branch_state_probs[the_homologs->get_num_homologs()-1][m][l]=0.0;
		}
	}
	delete[] transpoint_state_probs;
}  //End calc_transpoint_branch_probs





void Dupl_NoState_Base_model::print_transpoint_branch_probs(char *filename)
{
	int i, m,l;
	ofstream fout;

	calc_transpoint_branch_probs();
	
	
	fout.open(filename);
	fout<<"Site\t";
	for(m=0; m<curr_exchange->get_num_branches(); m++) {
			fout<<(*curr_tree)[m]->get_name()<<"\t";
			fout<<(*curr_tree)[m]->expect_subs_site()<<"\t";
			fout<<(*curr_tree)[m]->get_dist_to_root()<<"\t";
	}
	fout<<endl;
	fout<<"Site\t";
	for(m=0; m<curr_exchange->get_num_branches(); m++) 
	{
			fout<<"DUPL\tDUPL->SING\tSING\t";
	}
	fout<<endl;

	for(i=0; i<the_homologs->get_num_homologs(); i++) {
		fout<<i<<"\t";
		for(m=0; m<curr_exchange->get_num_branches(); m++) {
			for(l=0; l<NUM_TRANSPOINT_STATES; l++) {
				fout<<branch_state_probs[i][m][l];
				if ((m==curr_exchange->get_num_branches()-1) && (l==NUM_TRANSPOINT_STATES-1))
					fout<<"\n";
				else
					fout<<"\t";

			}
		}
	}


	fout.close();	
}
	



void Dupl_NoState_Base_model::calc_pattern_probs()
{
	int i, j, k, l, m, old_index;
	

	for(i=0; i<=curr_exchange->get_num_taxa(); i++) {
		for(j=0; j<n_choose_k(curr_exchange->get_num_taxa(), i); j++) {
				get_states(j, i, 0);
#ifdef _OPEN_MP_VERSION_
#pragma omp parallel for private(k,l,m, old_index)
#endif
			for(k=0; k<pow2[curr_exchange->get_num_taxa()]; k++) {
				if ((k & masks[i][j]) == 0) {
					//First instance of this redundancy class
					setup_data(i, dupl_states, k);
#ifdef _OPEN_MP_VERSION_
					pattern_probs[i][j][k]=log(prob_w_rate(omp_get_thread_num(),0));
#else
					pattern_probs[i][j][k]=log(prob_w_rate(0,0));
#endif
						//cout<<"Pattern prob for "<<i<<", "<<j<<", "<<k<<" : "<<pattern_probs[i][j][k]<<endl;
					if (do_transpoint_probs == (BOOL)TRUE) {
						for(l=0; l<curr_exchange->get_num_branches(); l++) {
							for (m=0; m<NUM_TRANSPOINT_STATES; m++) 
#ifdef _OPEN_MP_VERSION_	
								branch_transpoint_probs[i][j][k][l][m] = 
								calculate_transpoint_prob((*curr_tree)[l], int_to_transpoint_state(m), omp_get_thread_num(), 0);
#else
								branch_transpoint_probs[i][j][k][l][m] = 
								calculate_transpoint_prob((*curr_tree)[l], int_to_transpoint_state(m), 0, 0);
#endif
						}

					}

				}
									
			}
					
		}
		
	}
	
	
	for(i=0; i<=curr_exchange->get_num_taxa(); i++) {
		for(j=0; j<n_choose_k(curr_exchange->get_num_taxa(), i); j++) {
#ifdef _OPEN_MP_VERSION_
#pragma omp parallel for private(k,l,m, old_index)
#endif
			for(k=0; k<pow2[curr_exchange->get_num_taxa()]; k++) {
				if ((k & masks[i][j]) != 0) {
				//Look up previously determined prob
						old_index = k & (~masks[i][j]);   //0 out bits from duplicated taxa and look up
						pattern_probs[i][j][k]=pattern_probs[i][j][old_index];
						
						if (do_transpoint_probs == (BOOL)TRUE) {
							for(l=0; l<curr_exchange->get_num_branches(); l++) {
								
								for (m=0; m<NUM_TRANSPOINT_STATES; m++) 
									branch_transpoint_probs[i][j][k][l][m] = branch_transpoint_probs[i][j][old_index][l][m];
								
							}
							
						}
					
				}
				
			}
		}
	}
	
	

}



void Dupl_NoState_Base_model::build_masks()
{
	int i, j, k;

	//First index is number of duplicates 0->n
	masks=new int *[curr_exchange->get_num_taxa()+1];

	for(i=0; i<=curr_exchange->get_num_taxa(); i++) 
	//Next index is number of arrangements of i duplicates in n taxa
		masks[i] = new int [n_choose_k(curr_exchange->get_num_taxa(), i)];
		

	for(i=0; i<=curr_exchange->get_num_taxa(); i++) {
		for(j=0; j<n_choose_k(curr_exchange->get_num_taxa(), i); j++) {
			get_states(j, i, 0);
			
			masks[i][j]=0;
			for(k=1; k<curr_exchange->get_num_taxa(); k++) {
				if (dupl_states[k] == 2)
					masks[i][j] = masks[i][j] | pow2[k];
			}
		}
	}
	
}



void Dupl_NoState_Base_model::setup_data(int num_dupls, int *dupl_positions, int track_index)
{
	int i, site=0;
#ifdef _OPEN_MP_VERSION_
	site=omp_get_thread_num();
#endif
	
	for(i=0; i<curr_exchange->get_num_taxa(); i++) {
		if (dupl_positions[i] == 2) 
			(*curr_data)[i].Assign_site(site, loss_state_to_dupl(BOTH_PRESENT));
		else {
			if ((pow2[i] & track_index) == 0)
				(*curr_data)[i].Assign_site(site, loss_state_to_dupl(COPY1));
			else
				(*curr_data)[i].Assign_site(site, loss_state_to_dupl(COPY2));

		}

	}
}

void Dupl_NoState_Base_model::site_mask(int site_num)
{
	int i;

	transition_mask=0;
	for(i=0; i<curr_exchange->get_num_taxa(); i++) {
		if ((*the_tracks)[i].get_gene_track(site_num, 0)->my_locus == 0)
			//Gene must be on track 2
			transition_mask = transition_mask | pow2[i];
	}
	
}


void Dupl_NoState_Base_model::pair_site_mask(int this_site, int last_site)
{
	int i;

	transition_mask=0;
	
	for(i=0; i<curr_exchange->get_num_taxa(); i++) {
		if ((*the_tracks)[i].get_gene_track(last_site, 0)->my_locus != 0) {
			if (((*the_tracks)[i].get_gene_track(this_site, 0)->my_locus ==0) &&
				((*the_tracks)[i].get_gene_track(last_site, 1)->my_locus == 0))
				//If the only gene at this position is on the opposite track to the 
				//only gene at the last position, add a 1 to our mask at that position
				transition_mask = transition_mask | pow2[i];
		}
		else {
			if ((*the_tracks)[i].get_gene_track(this_site, 1)->my_locus ==0)
				//If the only gene at this position is on the opposite track to the 
				//only gene at the last position, add a 1 to our mask at that position
				transition_mask = transition_mask | pow2[i];
		}

	}
}


void Dupl_NoState_Base_model::build_transition_vector(int locus)
{
	build_transition_vector(locus, (BOOL)TRUE);
}


void Dupl_NoState_Base_model::build_transition_vector(int locus, BOOL look_back)
{
	int i, j;

#ifdef _OPEN_MP_VERSION_
#pragma omp parallel for private(i, j)
#endif
	for(i=0; i<pow2[curr_exchange->get_num_taxa()]; i++) {
		state_trans_probs[i]=1.0;
		for(j=0; j<curr_exchange->get_num_taxa(); j++) {
			if (look_back == (BOOL)TRUE) {
				if ((*the_tracks)[j].has_back_link(locus, 0) == (BOOL)FALSE) {
					if ((*the_tracks)[j].has_back_link(locus, 1) == (BOOL)FALSE) {
						state_trans_probs[i]*=0.5;
					}
					else {
						//if (((i ^ (curr_pattern ^ transition_mask)) & pow2[j]) == 0)
							if ((i & pow2[j]) == 0)
					
								state_trans_probs[i]*=
									pow(1.0-curr_exchange->get_strand_switch_prob(), 
										(*the_tracks)[j].get_dist_to_last(locus, 1));
							else
								state_trans_probs[i]*=1.0-
									pow(1.0-curr_exchange->get_strand_switch_prob(), 
										(*the_tracks)[j].get_dist_to_last(locus, 1));
					}

				}
				else { //Track 0 has link
					if ((*the_tracks)[j].has_back_link(locus, 1) == (BOOL)TRUE) {
						if ((*the_tracks)[j].get_dist_to_last(locus, 0) < 
							(*the_tracks)[j].get_dist_to_last(locus, 1)) {  //Track 0 is closer
							//if (((i ^ (curr_pattern ^ transition_mask)) & pow2[j]) == 0)
								if ((i & pow2[j]) == 0)
									state_trans_probs[i]*=
									pow(1.0-curr_exchange->get_strand_switch_prob(), 
										(*the_tracks)[j].get_dist_to_last(locus, 0));
							else
								state_trans_probs[i]*=1.0-
									pow(1.0-curr_exchange->get_strand_switch_prob(), 
										(*the_tracks)[j].get_dist_to_last(locus, 0));
						}
						else {  //Track 1 is closer
							//if (((i ^ (curr_pattern ^ transition_mask)) & pow2[j]) == 0)
							if ((i & pow2[j]) == 0)
					
								state_trans_probs[i]*=
									pow(1.0-curr_exchange->get_strand_switch_prob(), 
										(*the_tracks)[j].get_dist_to_last(locus, 1));
							else
								state_trans_probs[i]*=1.0-
									pow(1.0-curr_exchange->get_strand_switch_prob(), 
										(*the_tracks)[j].get_dist_to_last(locus, 1));
						}
					}
					else {  //Track 1 has no link
						//if (((i ^ (curr_pattern ^ transition_mask)) & pow2[j]) == 0)
						if ((i & pow2[j]) == 0)					
							state_trans_probs[i]*=
								pow(1.0-curr_exchange->get_strand_switch_prob(), 
									(*the_tracks)[j].get_dist_to_last(locus, 0));
						else
							state_trans_probs[i]*=1.0-
								pow(1.0-curr_exchange->get_strand_switch_prob(), 
									(*the_tracks)[j].get_dist_to_last(locus, 0));
					

					}
				}
				

			}
			else {
				if ((*the_tracks)[j].has_back_link(locus+1, 0)==(BOOL)FALSE) {
					if ((*the_tracks)[j].has_back_link(locus+1, 1) == (BOOL)FALSE) {
						state_trans_probs[i]*=0.5;
					}
					else {
						//if (((i ^ (curr_pattern ^ transition_mask)) & pow2[j]) == 0)
							if ((i & pow2[j]) == 0)
					
							state_trans_probs[i]*=
								pow(1.0-curr_exchange->get_strand_switch_prob(), 
									(*the_tracks)[j].get_dist_to_last(locus+1, 1));
						else	
							state_trans_probs[i]*=1.0-
								pow(1.0-curr_exchange->get_strand_switch_prob(), 
									(*the_tracks)[j].get_dist_to_last(locus+1, 1));
					}
				}
				else {  //Track 0 has link
					if ((*the_tracks)[j].has_back_link(locus+1, 1) == (BOOL)TRUE) {
						if ((*the_tracks)[j].get_dist_to_last(locus+1, 0) < 
							(*the_tracks)[j].get_dist_to_last(locus+1, 1)) {  //Track 0 is closer
								//if (((i ^ (curr_pattern ^ transition_mask)) & pow2[j]) == 0)
								if ((i & pow2[j]) == 0)
					
									state_trans_probs[i]*=pow(1.0-curr_exchange->get_strand_switch_prob(), 
										(*the_tracks)[j].get_dist_to_last(locus+1, 0));
								else
									state_trans_probs[i]*=1.0-pow(1.0-curr_exchange->get_strand_switch_prob(), 
										(*the_tracks)[j].get_dist_to_last(locus+1, 0));
						}
						else {//Track 1 is closer
								//if (((i ^ (curr_pattern ^ transition_mask)) & pow2[j]) == 0)
								if ((i & pow2[j]) == 0)
					
									state_trans_probs[i]*=pow(1.0-curr_exchange->get_strand_switch_prob(), 
										(*the_tracks)[j].get_dist_to_last(locus+1, 1));
								else
									state_trans_probs[i]*=1.0-pow(1.0-curr_exchange->get_strand_switch_prob(), 
										(*the_tracks)[j].get_dist_to_last(locus+1, 1));
						}
					}
					else {  //Track has no link
						//if (((i ^ (curr_pattern ^ transition_mask)) & pow2[j]) == 0)
						if ((i & pow2[j]) == 0)
					
							state_trans_probs[i]*=pow(1.0-curr_exchange->get_strand_switch_prob(), 
									(*the_tracks)[j].get_dist_to_last(locus+1, 0));
						else
							state_trans_probs[i]*=1.0-pow(1.0-curr_exchange->get_strand_switch_prob(), 
									(*the_tracks)[j].get_dist_to_last(locus+1, 0));
					}
				}

			}
		}
		state_trans_probs[i]=log(state_trans_probs[i]);
	}
}



double Dupl_NoState_Base_model::get_ln_sum_prob(double *vals)
{
	return(get_ln_sum_prob(vals, pow2[curr_exchange->get_num_taxa()]));
}

double Dupl_NoState_Base_model::get_ln_sum_prob_sort(double *vals)
{
	return(get_ln_sum_prob_sort(vals, pow2[curr_exchange->get_num_taxa()]));
}


double Dupl_NoState_Base_model::get_ln_sum_prob(double *vals, int num)
{
	int i, my_id;
	double smallest, cumulative;

#ifdef _OPEN_MP_VERSION_
	my_id=omp_get_thread_num();
#endif
	
	for(i=0; i<num; i++) {
#ifdef _OPEN_MP_VERSION_
		sum_probs_used_thread[my_id][i]=(BOOL)FALSE;
#else
		sum_probs_used[i]=(BOOL)FALSE;
#endif
	}
	
	cumulative=find_smallest(vals, num);

	for(i=1; i<num; i++)
	{
		smallest=find_smallest(vals, num);
		cumulative=logadd(cumulative, smallest);
		
	}

	return(cumulative);
}


double Dupl_NoState_Base_model::get_ln_sum_prob_sort(double *vals, int num)
{
	int i;
	double  cumulative;
	
	qsort(vals, num, sizeof(double), cmp);
	cumulative=vals[0];
	for(i=1; i<num; i++) cumulative=logadd(cumulative, vals[i]);
	
	
	return(cumulative);
}



double Dupl_NoState_Base_model::find_smallest(double *vals, int num)
{
	int i, start, index, my_id;
	double retval;
	
	
	start=0;
#ifdef _OPEN_MP_VERSION_
	my_id=omp_get_thread_num();
	while(sum_probs_used_thread[my_id][start] ==(BOOL)TRUE) start++;
#else
	while(sum_probs_used[start] ==(BOOL)TRUE) start++;
#endif
	index=start;
	retval=vals[start];
	for (i=start+1; i<num; i++)
	{
#ifdef _OPEN_MP_VERSION_
		if ((sum_probs_used_thread[my_id][i] ==(BOOL)FALSE) && (vals[i] < retval)) {
#else		
		if ((sum_probs_used[i] ==(BOOL)FALSE) && (vals[i] < retval)) {
#endif
			index=i;
			retval=vals[i];
		}
	}
#ifdef _OPEN_MP_VERSION_
	sum_probs_used_thread[my_id][index]=(BOOL)TRUE;	
#else
	sum_probs_used[index]=(BOOL)TRUE;
#endif
	return(retval);
}


double Dupl_NoState_Base_model::get_ln_sum_prob_zero(double *vals, int num)
{
	int i;
	double smallest, cumulative;

	for(i=0; i<num; i++)
#ifdef _OPEN_MP_VERSION_
		sum_probs_used_thread[omp_get_thread_num()][i]=(BOOL)FALSE;
#else
		sum_probs_used[i]=(BOOL)FALSE;
#endif

	cumulative=find_smallest(vals, num);

	for(i=1; i<num; i++)
	{
		smallest=find_smallest(vals, num);

		if (smallest != 1.0)                          //HACK--1.0 indicates 0 probability
		{
			if (cumulative != 1.0)
				cumulative=logadd(cumulative, smallest);
			else
				cumulative=smallest;
		}
		
	}

	return(cumulative);
}

double Dupl_NoState_Base_model::get_ln_sum_prob_zero_sort(double *vals, int num)
{
	int i;
	double cumulative;
	
	qsort(vals, num, sizeof(double), cmp);
	cumulative=vals[0];
	
	for(i=1; i<num; i++) {
		
		if (vals[i] != 1.0)                          //HACK--1.0 indicates 0 probability
		{
			if (cumulative != 1.0)
				cumulative=logadd(cumulative, vals[i]);
			else
				cumulative=vals[i];
		}
		
	}
	
	return(cumulative);
}
	



double Dupl_NoState_Base_model::find_smallest_zero(double *vals)
{
	int i, start, index;
	double retval;
	
	start=0;
	while(sum_probs_used[start] ==(BOOL)TRUE) start++;
	
	index=start;
	retval=vals[start];
	for (i=start+1; i<pow2[curr_exchange->get_num_taxa()]; i++)
	{
		if ((sum_probs_used[i] ==(BOOL)FALSE) && ((vals[i] < retval) || (vals[i] == 1.0))) {
			index=i;
			retval=vals[i];
		}
	}

	sum_probs_used[index]=(BOOL)TRUE;
	return(retval);
}



void Dupl_NoState_Base_model::get_site_state_probs(int taxa_id)
{
	int i, j, k, l, m, num_prob_indices, last_state;
	double sum, sum_prob, partial_prob, *left_trans_probs, *right_trans_probs;


	get_gene_conditional_probs();
	site_probs=new double * [the_homologs->get_num_homologs()];
	left_trans_probs=new double[pow2[curr_exchange->get_num_taxa()]];
	right_trans_probs=new double[pow2[curr_exchange->get_num_taxa()]];

	
	for(i=0; i<the_homologs->get_num_homologs(); i++)
		site_probs[i] = new double [curr_exchange->get_condlike_size()];

	for(i=0; i<the_homologs->get_num_homologs(); i++)
		for(j=0; j<curr_exchange->get_condlike_size(); j++)
			site_probs[i][j]=0.0;

	for(i=0; i<the_homologs->get_num_homologs(); i++) {
		//This function returns the number of valid probability indices for this locus and their 
		//indicies (dupl_indices and site_indicies
		num_prob_indices=get_locus_state_index(i);
		site_mask(i);
		build_transition_vector(i, (BOOL)TRUE);
	
		for(l=0; l<pow2[curr_exchange->get_num_taxa()]; l++) 
			left_trans_probs[l]=state_trans_probs[l];

		build_transition_vector(i, (BOOL)FALSE);
		for(l=0; l<pow2[curr_exchange->get_num_taxa()]; l++) 
			right_trans_probs[l]=state_trans_probs[l];


		sum=0.0;
		last_state=-1;
		for(j=0; j<num_prob_indices; j++) {
			get_states(site_indices[j], dupl_indices[j], 0);
			if (last_state != dupl_states[taxa_id]) {
				last_state=dupl_states[taxa_id];
				
				for(k=0; k<pow2[curr_exchange->get_num_taxa()]; k++) {
					setup_data(dupl_indices[j], dupl_states, k);
				
					if (i!=0) {
						for(l=0; l<pow2[curr_exchange->get_num_taxa()]; l++) {
							sum_probs[l]=left_trans_probs[k ^ l]+left_cond_probs[i-1][l];
						}
					
						sum_prob = get_ln_sum_prob(sum_probs);
					}
					else sum_prob = 0;

					if (i != the_homologs->get_num_homologs()-1) {
						for(l=0; l<pow2[curr_exchange->get_num_taxa()]; l++) {
							sum_probs[l]=right_trans_probs[k ^ l]+right_cond_probs[i+1][l];
						}
						sum_prob = sum_prob + get_ln_sum_prob(sum_probs);
					}
					
					if (dupl_to_loss_state((*curr_data)[taxa_id][0]) == BOTH_PRESENT) {
						for(l=0; l<curr_exchange->get_condlike_size(); l++) {
							switch (dupl_to_loss_state(l)) {
							case BOTH_PRESENT:
							case BOTH_1_BIAS:
							case BOTH_2_BIAS:
							case BOTH_PRESENT_FIXED:
							case BOTH_FIXED_SUBF:
							if (allowed_state(dupl_to_loss_state(l)) == (BOOL)TRUE) {
									(*curr_data)[taxa_id].Assign_site(0, l);
									partial_prob=0;

									partial_prob_w_rate_nonhidden(0, curr_tree->get_leftmost_tip(), 0, 
										curr_tree->find_root(), taxa_id); 
									//CHANGE LOGADD APPROACH!!!!!
									for (m=0; m<curr_exchange->get_condlike_size();m++)
										partial_prob+=root_freq(loss_state_to_dupl(BOTH_PRESENT))*
										curr_tree->find_root()->get_trpb(0, loss_state_to_dupl(BOTH_PRESENT), m)*
										curr_tree->find_root()->get_cond_prob(m);
										
									if (site_probs[i][l] == 0.0)
										site_probs[i][l]=log(partial_prob)+sum_prob;
									else
										site_probs[i][l]=logadd(site_probs[i][l], log(partial_prob)+sum_prob);
										
									if (sum == 0.0)
										sum=log(partial_prob)+sum_prob;
									else
										sum=logadd(sum, log(partial_prob)+sum_prob);
								}
								break;
							}
						}//...for (l=0; l<curr_exchange->get_condlike_size...
					}
					else if (dupl_to_loss_state((*curr_data)[taxa_id][0]) == COPY1)  {
						for(l=0; l<curr_exchange->get_condlike_size(); l++) {
							switch (dupl_to_loss_state(l)) {
							case COPY1:
							case COPY1_BIAS:
								if (allowed_state(dupl_to_loss_state(l)) == (BOOL)TRUE) {
									(*curr_data)[taxa_id].Assign_site(0, l);
									partial_prob=0;

									partial_prob_w_rate_nonhidden(0, curr_tree->get_leftmost_tip(), 0, 
										curr_tree->find_root(), taxa_id); 
									//CHANGE LOGADD APPROACH!!!!!
									for (m=0; m<curr_exchange->get_condlike_size();m++)
										partial_prob+=root_freq(loss_state_to_dupl(BOTH_PRESENT))*
										curr_tree->find_root()->get_trpb(0, loss_state_to_dupl(BOTH_PRESENT), m)*
										curr_tree->find_root()->get_cond_prob(m);
										
									if (site_probs[i][l] == 0.0)
										site_probs[i][l]=log(partial_prob)+sum_prob;
									else
										site_probs[i][l]=logadd(site_probs[i][l], log(partial_prob)+sum_prob);
										
									if (sum == 0.0)
										sum=log(partial_prob)+sum_prob;
									else
										sum=logadd(sum, log(partial_prob)+sum_prob);
								}
								break;
							}
						}//...for (l=0; l<curr_exchange->get_condlike_size...
					}
					else if (dupl_to_loss_state((*curr_data)[taxa_id][0]) == COPY2)  {
						for(l=0; l<curr_exchange->get_condlike_size(); l++) {
							switch (dupl_to_loss_state(l)) {
							case COPY2:
							case COPY2_BIAS:
								if (allowed_state(dupl_to_loss_state(l)) == (BOOL)TRUE) {
									(*curr_data)[taxa_id].Assign_site(0, l);
									partial_prob=0;

									partial_prob_w_rate_nonhidden(0, curr_tree->get_leftmost_tip(), 0, 
										curr_tree->find_root(), taxa_id); 
									//CHANGE LOGADD APPROACH!!!!!
									for (m=0; m<curr_exchange->get_condlike_size();m++)
										partial_prob+=root_freq(loss_state_to_dupl(BOTH_PRESENT))*
										curr_tree->find_root()->get_trpb(0, loss_state_to_dupl(BOTH_PRESENT), m)*
										curr_tree->find_root()->get_cond_prob(m);
										
									if (site_probs[i][l] == 0.0)
										site_probs[i][l]=log(partial_prob)+sum_prob;
									else
										site_probs[i][l]=logadd(site_probs[i][l], log(partial_prob)+sum_prob);
										
									if (sum == 0.0)
										sum=log(partial_prob)+sum_prob;
									else
										sum=logadd(sum, log(partial_prob)+sum_prob);
								}
								break;
							}
						}//...for (l=0; l<curr_exchange->get_condlike_size...
					}
				}//...for(k=0; k<pow2[num_taxa]
			}//...if last_state != dupl_states...
		}  //..for(j=0; j<num_prob_indices

		for(j=0; j<curr_exchange->get_condlike_size(); j++) {
			if (site_probs[i][j] != 0.0)
				site_probs[i][j] = exp(site_probs[i][j] - sum);
			//site_probs[i][j]=site_probs[i][j];
		}
		
	} //..for(i=0; i<num_localities

	delete[] right_trans_probs;
	delete[] left_trans_probs;
	


}



Dupl_NoState_Base_model::~Dupl_NoState_Base_model()
{
	int i, j, k, l;

	for(i=0; i<=curr_exchange->get_num_taxa(); i++) {
	//Next index is number of arrangements of i duplicates in n taxa
		for(j=0; j<n_choose_k(curr_exchange->get_num_taxa(), i); j++)			
			delete[] pattern_probs[i][j];

		delete[] pattern_probs[i];
	}

	for(i=0; i<curr_exchange->get_num_taxa()+1; i++) {
		delete[] masks[i];
	}

	if (post_probs != 0) {
		for(i=0; i<the_homologs->get_num_homologs(); i++)
			delete[] post_probs[i];
		delete[] post_probs;

	}

	if(branch_state_probs !=0) {
		for(i=0; i<the_homologs->get_num_homologs(); i++) {
			for(j=0; j<curr_exchange->get_num_branches(); j++)
				delete[] branch_state_probs[i][j];
			delete[] branch_state_probs[i];
		}
		delete[] branch_state_probs;
		
	}

	if(do_transpoint_probs == (BOOL)TRUE) {
		delete[] transpoint_condprobs1;
		delete[] transpoint_condprobs2;
		
#ifdef _OPEN_MP_VERSION_
		for (i=0; i<curr_exchange->get_num_open_mp_threads(); i++) {
			delete[] transpoint_condprobs1_thread[i];
			delete[] transpoint_condprobs2_thread[i];
		}
		delete[] transpoint_condprobs1_thread;
		delete[] transpoint_condprobs2_thread;
#endif
		for(i=0; i<=curr_exchange->get_num_taxa(); i++) {
			for(j=0; j<n_choose_k(curr_exchange->get_num_taxa(), i); j++) {
				for(k=0; k<pow2[curr_exchange->get_num_taxa()]; k++) {
					for(l=0; l< curr_exchange->get_num_branches(); l++)
						delete[] branch_transpoint_probs[i][j][k][l];
					delete[] branch_transpoint_probs[i][j][k];	
				}
				delete[] branch_transpoint_probs[i][j];
			}
			delete[] branch_transpoint_probs[i];
		}
		delete[] branch_transpoint_probs;
	}

	delete[] pattern_probs;
	delete[] masks;
	delete[] dupl_states;
	delete[] state_trans_probs;
	delete[] cumulative_probs;
	delete[] last_cumulative_probs;
	delete[] sum_probs;
	
#ifdef _OPEN_MP_VERSION_
	for (i=0; i<curr_exchange->get_num_open_mp_threads(); i++) delete[] sum_probs_thread[i];
	delete[] sum_probs_thread;
#endif
	
	
	delete[] sum_probs_r;
	delete[] sum_probs_used;
	
#ifdef _OPEN_MP_VERSION_
	for (i=0; i<curr_exchange->get_num_open_mp_threads(); i++) delete[] sum_probs_used_thread[i];
	delete[] sum_probs_used_thread;
#endif
	
	delete[] site_indices;
	delete[] dupl_indices;
	delete[] is_ambig;

	if (left_cond_probs !=0) {
		for(i=0; i<the_homologs->get_num_homologs(); i++) {
			delete[] left_cond_probs[i];
			delete[] right_cond_probs[i];
		}

		delete[] left_cond_probs;
		delete[] right_cond_probs;
	}

/*	if(site_probs !=0) {
		for(i=0; i<the_homologs->get_num_homologs(); i++)
			delete[] site_probs[i];
		delete[] site_probs;
	}*/
}



void Dupl_Null_Base_model::calc_transprobs(Branch *taxa, int rate_num)
{
	int i, j;

		//Set transprobs so that unused states in this model are never entered
		for (i=0;i<curr_exchange->get_condlike_size(); i++){
			if (dupl_to_loss_state(i) != LOST) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(LOST), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(LOST), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(LOST),loss_state_to_dupl(LOST), 1.0);

			if (dupl_to_loss_state(i) != BOTH_PRESENT_FIXED) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(BOTH_PRESENT_FIXED), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED),
					loss_state_to_dupl(BOTH_PRESENT_FIXED), 1.0);

			if (dupl_to_loss_state(i) != BOTH_1_BIAS) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(BOTH_1_BIAS), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS),
					loss_state_to_dupl(BOTH_1_BIAS), 1.0);

			if (dupl_to_loss_state(i) != BOTH_2_BIAS) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(BOTH_2_BIAS), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS),
					loss_state_to_dupl(BOTH_2_BIAS), 1.0);

			if (dupl_to_loss_state(i) != GENERIC_SINGLE_COPY) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(GENERIC_SINGLE_COPY), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(GENERIC_SINGLE_COPY), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(GENERIC_SINGLE_COPY),
					loss_state_to_dupl(GENERIC_SINGLE_COPY), 1.0);

			if (dupl_to_loss_state(i) != MISSING) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(MISSING), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(MISSING), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(MISSING),
					loss_state_to_dupl(MISSING), 1.0);

			if (dupl_to_loss_state(i) != COPY1_OR_BOTH) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(COPY1_OR_BOTH), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1_OR_BOTH), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1_OR_BOTH),
					loss_state_to_dupl(COPY1_OR_BOTH), 1.0);

			if (dupl_to_loss_state(i) != COPY2_OR_BOTH) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(COPY2_OR_BOTH), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2_OR_BOTH), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2_OR_BOTH),
					loss_state_to_dupl(COPY2_OR_BOTH), 1.0);
			
			if (dupl_to_loss_state(i) != BOTH_FIXED_SUBF) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(BOTH_FIXED_SUBF), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_FIXED_SUBF), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_FIXED_SUBF),loss_state_to_dupl(BOTH_FIXED_SUBF), 1.0);

			if (dupl_to_loss_state(i) != COPY1_BIAS) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(COPY1_BIAS), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1_BIAS), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1_BIAS),loss_state_to_dupl(COPY1_BIAS), 1.0);

			if (dupl_to_loss_state(i) != COPY2_BIAS) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(COPY2_BIAS), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2_BIAS), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2_BIAS),loss_state_to_dupl(COPY2_BIAS), 1.0);


			

		}

		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(COPY1), 1.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(COPY2), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(COPY2), 1.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(COPY1), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(COPY1),
				0.5-0.5*exp(-2.0*taxa->get_brnlen()));
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(COPY2),
				0.5-0.5*exp(-2.0*taxa->get_brnlen()));
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_PRESENT), 
			exp(-2.0*taxa->get_brnlen()));


}



int Dupl_Null_Base_model::state_redundancy_val(DUPL_LOSS_STATES the_state)
{
	//Tells us if we need to sum over internal model states for a given input tip state

	if (observable_dupl_state(the_state)==(BOOL)TRUE) {
		switch (the_state) {
		case COPY1:
		case COPY2:
		case LOST:
			return(1);
			break;
		case BOTH_PRESENT:
			return(1);
			break;
		case MISSING:
			return(3);
			break;
		case GENERIC_SINGLE_COPY:
			return(2);
			break;
		case COPY1_OR_BOTH:
			return(2);
			break;
		case COPY2_OR_BOTH:
			return(2);
			break;
		}
	}
	else {
		cerr<<"ERROR: Tried to get Tip redundancy for non-observable internal model state\n";
		return(1);
	}

}


DUPL_LOSS_STATES Dupl_Null_Base_model::get_redund_position_n(DUPL_LOSS_STATES the_state, int value)
{
	//Indicates order of the states to be summed through for ambiguous observed states

	switch(the_state) {
	
	//If it is not an ambiguous state, just return it
	case COPY1:
	case COPY2:
	case BOTH_PRESENT:
	case LOST:
		return(the_state);
		break;

	//Otherwise, give the possible states in an arbitrary order
	case MISSING:
		if (value == 0)
			return(BOTH_PRESENT);
		else if (value == 1)
			return(COPY1);
		else
			return(COPY2);
		break;
	case GENERIC_SINGLE_COPY:
		if (value == 0)
			return(COPY1);
		else
			return(COPY2);
		break;
	case COPY1_OR_BOTH:
		if (value == 0)
			return(COPY1);
		else
			return(BOTH_PRESENT);
		break;
	case COPY2_OR_BOTH:
		if (value == 0)
			return(COPY2);
		else
			return(BOTH_PRESENT);
		break;

	}

}


double Dupl_Null_Base_model::find_ut (Branch *taxa)
{
	return(taxa->expect_subs_site()/2.0);
}


void Dupl_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
	int i, brn_cnt=0;

	brn_start=0;
	for(i=0; i<curr_exchange->get_num_branches(); i++) {
		if ((curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE) || ((*curr_tree)[i]->get_brnlen() != 0)) {
				par[brn_cnt]=(*curr_tree)[i]->get_brnlen();
				types[brn_cnt]=BRANCH;
				if (curr_exchange->zero_len_brns_fixed() ==(BOOL)TRUE)
					brn_index[brn_cnt]=i;
				brn_cnt++;
		}
	}
}

void Dupl_model::describe_results()
{

}


void Dupl_model::num_params_model() 
{
	int i, cnt_zero=0;

	if (curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE)
		curr_exchange->set_num_params(curr_exchange->get_num_branches());
	else {
		for(i=0; i<curr_exchange->get_num_branches(); i++)
		{
			if ((*curr_tree)[i]->expect_subs_site() ==0)
				cnt_zero++;
		}
		brn_index=new int[curr_exchange->get_num_branches()-cnt_zero];
		curr_exchange->set_num_params(curr_exchange->get_num_branches()-cnt_zero);
	}
}






void Dupl_NoState_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
	int i, brn_cnt=0;

	par[0]=log(curr_exchange->get_strand_switch_prob());
	types[0]=STRAND_SWITCH_PROB;

	brn_start=1;
	for(i=0; i<curr_exchange->get_num_branches(); i++) {
		if ((curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE) || ((*curr_tree)[i]->get_brnlen() != 0)) {
				par[brn_cnt+brn_start]=(*curr_tree)[i]->get_brnlen();
				types[brn_cnt+brn_start]=BRANCH;
				if (curr_exchange->zero_len_brns_fixed() ==(BOOL)TRUE)
					brn_index[brn_cnt]=i;
				brn_cnt++;
		}
	}
}



void Dupl_NoState_model::num_params_model() 
{
	int i, cnt_zero=0;

	if (curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE)
		curr_exchange->set_num_params(curr_exchange->get_num_branches()+1);
	else {
		for(i=0; i<curr_exchange->get_num_branches(); i++)
		{
			if ((*curr_tree)[i]->expect_subs_site() ==0)
				cnt_zero++;
		}
		brn_index=new int[curr_exchange->get_num_branches()-cnt_zero];
		curr_exchange->set_num_params(curr_exchange->get_num_branches()-cnt_zero+1);
	}
}




void Dupl_Fix_Base_model::calc_transprobs(Branch *taxa, int rate_num)
{
		int i, j;

		//Set transprobs so that unused states in this model are never entered
		for (i=0;i<curr_exchange->get_condlike_size(); i++){
			if (dupl_to_loss_state(i) != LOST) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(LOST), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(LOST), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(LOST),loss_state_to_dupl(LOST), 1.0);

			if (dupl_to_loss_state(i) != BOTH_1_BIAS) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(BOTH_1_BIAS), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS),
					loss_state_to_dupl(BOTH_1_BIAS), 1.0);

			if (dupl_to_loss_state(i) != BOTH_2_BIAS) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(BOTH_2_BIAS), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS),
					loss_state_to_dupl(BOTH_2_BIAS), 1.0);

			if (dupl_to_loss_state(i) != GENERIC_SINGLE_COPY) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(GENERIC_SINGLE_COPY), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(GENERIC_SINGLE_COPY), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(GENERIC_SINGLE_COPY),
					loss_state_to_dupl(GENERIC_SINGLE_COPY), 1.0);

			if (dupl_to_loss_state(i) != MISSING) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(MISSING), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(MISSING), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(MISSING),
					loss_state_to_dupl(MISSING), 1.0);

			if (dupl_to_loss_state(i) != COPY1_OR_BOTH) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(COPY1_OR_BOTH), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1_OR_BOTH), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1_OR_BOTH),
					loss_state_to_dupl(COPY1_OR_BOTH), 1.0);

			if (dupl_to_loss_state(i) != COPY2_OR_BOTH) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(COPY2_OR_BOTH), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2_OR_BOTH), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2_OR_BOTH),
					loss_state_to_dupl(COPY2_OR_BOTH), 1.0);

				if (dupl_to_loss_state(i) != BOTH_FIXED_SUBF) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(BOTH_FIXED_SUBF), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_FIXED_SUBF), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_FIXED_SUBF),loss_state_to_dupl(BOTH_FIXED_SUBF), 1.0);

			if (dupl_to_loss_state(i) != COPY1_BIAS) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(COPY1_BIAS), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1_BIAS), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1_BIAS),loss_state_to_dupl(COPY1_BIAS), 1.0);

			if (dupl_to_loss_state(i) != COPY2_BIAS) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(COPY2_BIAS), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2_BIAS), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2_BIAS),loss_state_to_dupl(COPY2_BIAS), 1.0);
			

		}

		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(COPY1), 1.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(COPY2), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_PRESENT_FIXED), 0.0);

		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(COPY2), 1.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(COPY1), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_PRESENT_FIXED), 0.0);

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(COPY2), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(COPY1), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(BOTH_PRESENT_FIXED), 1.0);

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(COPY1),
				(1.0-exp(-1.0*taxa->get_brnlen()*(2.0+curr_exchange->get_dupl_fix_rate()) ) )/ 
				(2.0+curr_exchange->get_dupl_fix_rate()) );
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(COPY2),
				(1.0-exp(-1.0*taxa->get_brnlen()*(2.0+curr_exchange->get_dupl_fix_rate()) ) )/ 
				(2.0+curr_exchange->get_dupl_fix_rate()) );
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_PRESENT_FIXED),
				(curr_exchange->get_dupl_fix_rate()-curr_exchange->get_dupl_fix_rate()*
				exp(-1.0*taxa->get_brnlen()*(2.0+curr_exchange->get_dupl_fix_rate()) ) )/
				(2.0+curr_exchange->get_dupl_fix_rate()) );
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_PRESENT), 
			exp(-1.0*taxa->get_brnlen()*(2.0+curr_exchange->get_dupl_fix_rate())));



}


int Dupl_Fix_Base_model::state_redundancy_val(DUPL_LOSS_STATES the_state)
{
		//Tells us if we need to sum over internal model states for a given input tip state

	if (observable_dupl_state(the_state)==(BOOL)TRUE) {
		switch (the_state) {
		case COPY1:
		case COPY2:
		case LOST:
			return(1);
			break;
		case BOTH_PRESENT:
			return(2);
			break;
		case MISSING:
			return(4);
			break;
		case GENERIC_SINGLE_COPY:
			return(2);
			break;
		case COPY1_OR_BOTH:
		case COPY2_OR_BOTH:
			return(3);
			break;
		}
	}
	else {
		cerr<<"ERROR: Tried to get Tip redundancy for non-observable internal model state\n";
		return(1);
	}

}



DUPL_LOSS_STATES Dupl_Fix_Base_model::get_redund_position_n(DUPL_LOSS_STATES the_state, int value)
{
		//Indicates order of the states to be summed through for ambiguous observed states

	switch(the_state) {
	
	//If it is not an ambiguous state, just return it
	case COPY1:
	case COPY2:
	case LOST:
		return(the_state);
		break;

	//Otherwise, give the possible states in an arbitrary order
	case BOTH_PRESENT:
		if (value == 0)
			return(BOTH_PRESENT);
		else
			return(BOTH_PRESENT_FIXED);
		break;
	case MISSING:
		if (value == 0)
			return(BOTH_PRESENT);
		else if (value == 1)
			return(COPY1);
		else if (value == 2)
			return(COPY2);
		else
			return(BOTH_PRESENT_FIXED);
		break;
	case GENERIC_SINGLE_COPY:
		if (value == 0)
			return(COPY1);
		else
			return(COPY2);
		break;
	case COPY1_OR_BOTH:
		if (value == 0)
			return(COPY1);
		else if (value == 1)
			return(BOTH_PRESENT);
		else
			return(BOTH_PRESENT_FIXED);
		break;
	case COPY2_OR_BOTH:
		if (value == 0)
			return(COPY2);
		else if (value == 1)
			return(BOTH_PRESENT);
		else
			return(BOTH_PRESENT_FIXED);
		break;

	}


}


BOOL Dupl_Fix_Base_model::allowed_state(DUPL_LOSS_STATES the_state)
{
	BOOL retval=(BOOL)FALSE;

	switch (the_state) {
	case BOTH_PRESENT:
	case BOTH_PRESENT_FIXED:
	case COPY1:
	case COPY2:
	case COPY1_OR_BOTH:
	case COPY2_OR_BOTH:
	case GENERIC_SINGLE_COPY:
	case MISSING:
	case LOST:
		retval=(BOOL)TRUE;
		break;

	}
	return(retval);
}



double Dupl_Fix_Base_model::find_ut (Branch *taxa)
{
	return(taxa->expect_subs_site()/(2.0+curr_exchange->get_dupl_fix_rate()));
}



Dupl_Fix_model::Dupl_Fix_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree)
{
	if (cexchange->get_dupl_fix_rate() == 0)
		cexchange->set_dupl_fix_rate(0.1); 

	assemble (cexchange, cdata, ctree);
}



void Dupl_Fix_model::describe_results()
{
	cout<<"Instantaneous relative rate of duplicate fixation (multiply by two for comparison to branch lengths): "
		<<curr_exchange->get_dupl_fix_rate()<<"\n";
}




void Dupl_Fix_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
	int i, brn_cnt=0;

	brn_start=0;
	for(i=0; i<curr_exchange->get_num_branches(); i++) {
		if ((curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE) || ((*curr_tree)[i]->get_brnlen() != 0)) {
				par[brn_cnt]=(*curr_tree)[i]->get_brnlen();
				types[brn_cnt]=BRANCH;
				if (curr_exchange->zero_len_brns_fixed() ==(BOOL)TRUE)
					brn_index[brn_cnt]=i;
				brn_cnt++;
		}
	}

	par[curr_exchange->get_num_params()-1]=curr_exchange->get_dupl_fix_rate();
	types[curr_exchange->get_num_params()-1]=DUPL_FIX_RATE;
}
    



void Dupl_Fix_model::num_params_model()
{
	int i, cnt_zero=0;

	if (curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE)
		curr_exchange->set_num_params(curr_exchange->get_num_branches()+1);
	else {
		for(i=0; i<curr_exchange->get_num_branches(); i++)
		{
			if ((*curr_tree)[i]->expect_subs_site() ==0)
				cnt_zero++;
		}
		brn_index=new int[curr_exchange->get_num_branches()-cnt_zero];
		curr_exchange->set_num_params(curr_exchange->get_num_branches()-cnt_zero+1);
	}

}


Dupl_Fix_NoState_model::Dupl_Fix_NoState_model (Exchange *cexchange, Tree *ctree, Clade *cgenomes, WGD_Data *chomologs) 
{
	allocate_state_model(cexchange, ctree, cgenomes, chomologs);

	if (cexchange->get_dupl_fix_rate() == 0)
		cexchange->set_dupl_fix_rate(0.1);	

	left_cond_probs=right_cond_probs=0;

	assemble (cexchange, curr_data, ctree);
}



void Dupl_Fix_NoState_model::describe_results()
{
	cout<<"Instantaneous relative rate of duplicate fixation (multiply by two for comparison to branch lengths): "
		<<curr_exchange->get_dupl_fix_rate()<<"\n";
	cout<<"Probability of strand definition switching between two single copy genes: "
		<<curr_exchange->get_strand_switch_prob()<<endl;
}




void Dupl_Fix_NoState_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
	int i, brn_cnt=0;

	brn_start=0;
	for(i=0; i<curr_exchange->get_num_branches(); i++) {
		if ((curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE) || ((*curr_tree)[i]->get_brnlen() != 0)) {
				par[brn_cnt]=(*curr_tree)[i]->get_brnlen();
				types[brn_cnt]=BRANCH;
				if (curr_exchange->zero_len_brns_fixed() ==(BOOL)TRUE)
					brn_index[brn_cnt]=i;
				brn_cnt++;
		}
	}

	par[curr_exchange->get_num_params()-2]=curr_exchange->get_dupl_fix_rate();
	types[curr_exchange->get_num_params()-2]=DUPL_FIX_RATE;

	par[curr_exchange->get_num_params()-1]=log(curr_exchange->get_strand_switch_prob());
	types[curr_exchange->get_num_params()-1]=STRAND_SWITCH_PROB;
}
    



void Dupl_Fix_NoState_model::num_params_model()
{
	int i, cnt_zero=0;

	if (curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE)
		curr_exchange->set_num_params(curr_exchange->get_num_branches()+2);
	else {
		for(i=0; i<curr_exchange->get_num_branches(); i++)
		{
			if ((*curr_tree)[i]->expect_subs_site() ==0)
				cnt_zero++;
		}
		brn_index=new int[curr_exchange->get_num_branches()-cnt_zero];
		curr_exchange->set_num_params(curr_exchange->get_num_branches()-cnt_zero+2);
	}

}


void Dupl_Parallel_Base_model::calc_transprobs(Branch *taxa, int rate_num)
{
		int i, j;
		double one_plus_beta, one_plus_two_beta;

		//Set transprobs so that unused states in this model are never entered
		for (i=0;i<curr_exchange->get_condlike_size(); i++){
			if (dupl_to_loss_state(i) != LOST) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(LOST), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(LOST), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(LOST),loss_state_to_dupl(LOST), 1.0);

			if (dupl_to_loss_state(i) != BOTH_PRESENT_FIXED) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(BOTH_PRESENT_FIXED), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED),
					loss_state_to_dupl(BOTH_PRESENT_FIXED), 1.0);

			if (dupl_to_loss_state(i) != GENERIC_SINGLE_COPY) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(GENERIC_SINGLE_COPY), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(GENERIC_SINGLE_COPY), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(GENERIC_SINGLE_COPY),
					loss_state_to_dupl(GENERIC_SINGLE_COPY), 1.0);

			if (dupl_to_loss_state(i) != MISSING) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(MISSING), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(MISSING), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(MISSING),
					loss_state_to_dupl(MISSING), 1.0);

			if (dupl_to_loss_state(i) != COPY1_OR_BOTH) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(COPY1_OR_BOTH), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1_OR_BOTH), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1_OR_BOTH),
					loss_state_to_dupl(COPY1_OR_BOTH), 1.0);

			if (dupl_to_loss_state(i) != COPY2_OR_BOTH) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(COPY2_OR_BOTH), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2_OR_BOTH), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2_OR_BOTH),
					loss_state_to_dupl(COPY2_OR_BOTH), 1.0);

				if (dupl_to_loss_state(i) != BOTH_FIXED_SUBF) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(BOTH_FIXED_SUBF), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_FIXED_SUBF), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_FIXED_SUBF),loss_state_to_dupl(BOTH_FIXED_SUBF), 1.0);

			if (dupl_to_loss_state(i) != COPY1_BIAS) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(COPY1_BIAS), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1_BIAS), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1_BIAS),loss_state_to_dupl(COPY1_BIAS), 1.0);

			if (dupl_to_loss_state(i) != COPY2_BIAS) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(COPY2_BIAS), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2_BIAS), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2_BIAS),loss_state_to_dupl(COPY2_BIAS), 1.0);
			

		}

		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(COPY1), 1.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(COPY2), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_2_BIAS), 0.0);

		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(COPY2), 1.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(COPY1), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_2_BIAS), 0.0);

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_2_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(COPY1), 1.0-exp(-1.0*taxa->get_brnlen()));
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(COPY2), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_1_BIAS), exp(-1.0*taxa->get_brnlen()));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(COPY1), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(COPY2), 1.0-exp(-1.0*taxa->get_brnlen()));
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_2_BIAS), exp(-1.0*taxa->get_brnlen()));

		
		one_plus_beta=1.0+curr_exchange->get_dupl_parallel_rate();
		one_plus_two_beta=1.0+2.0*curr_exchange->get_dupl_parallel_rate();


		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(COPY1),
				 (one_plus_two_beta - exp(-2.0*taxa->get_brnlen()*one_plus_beta) -
				 2.0*curr_exchange->get_dupl_parallel_rate()*exp(-1.0*taxa->get_brnlen()) )/ 
				 (2*one_plus_two_beta));
		
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(COPY2),
				 (one_plus_two_beta - exp(-2.0*taxa->get_brnlen()*one_plus_beta) -
				 2.0*curr_exchange->get_dupl_parallel_rate()*exp(-1.0*taxa->get_brnlen()) )/ 
				 (2*one_plus_two_beta));
		
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_1_BIAS),
				(curr_exchange->get_dupl_parallel_rate()/one_plus_two_beta)*
				( exp(-1.0*taxa->get_brnlen()) - exp(-2.0*taxa->get_brnlen()*one_plus_beta) ));
		
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_2_BIAS),
				(curr_exchange->get_dupl_parallel_rate()/one_plus_two_beta)*
				( exp(-1.0*taxa->get_brnlen()) - exp(-2.0*taxa->get_brnlen()*one_plus_beta) ));
		
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_PRESENT), 
			exp(-2.0*taxa->get_brnlen()*one_plus_beta));



}


int Dupl_Parallel_Base_model::state_redundancy_val(DUPL_LOSS_STATES the_state)
{
		//Tells us if we need to sum over internal model states for a given input tip state

	if (observable_dupl_state(the_state)==(BOOL)TRUE) {
		switch (the_state) {
		case COPY1:
		case COPY2:
		case LOST:
			return(1);
			break;
		case BOTH_PRESENT:
			return(3);
			break;
		case MISSING:
			return(5);
			break;
		case GENERIC_SINGLE_COPY:
			return(2);
			break;
		case COPY1_OR_BOTH:
		case COPY2_OR_BOTH:
			return(4);
			break;
		}
	}
	else {
		cerr<<"ERROR: Tried to get Tip redundancy for non-observable internal model state\n";
		return(1);
	}

}



DUPL_LOSS_STATES Dupl_Parallel_Base_model::get_redund_position_n(DUPL_LOSS_STATES the_state, int value)
{
		//Indicates order of the states to be summed through for ambiguous observed states

	switch(the_state) {
	
	//If it is not an ambiguous state, just return it
	case COPY1:
	case COPY2:
	case LOST:
		return(the_state);
		break;

	//Otherwise, give the possible states in an arbitrary order
	case BOTH_PRESENT:
		if (value == 0)
			return(BOTH_PRESENT);
		else if (value == 1)
			return(BOTH_1_BIAS);
		else
			return(BOTH_2_BIAS);
		break;
	case MISSING:
		if (value == 0)
			return(BOTH_PRESENT);
		else if (value == 1)
			return(COPY1);
		else if (value == 2)
			return(COPY2);
		else if (value == 3)
			return(BOTH_1_BIAS);
		else
			return(BOTH_2_BIAS);
		break;
	case GENERIC_SINGLE_COPY:
		if (value == 0)
			return(COPY1);
		else
			return(COPY2);
		break;
	case COPY1_OR_BOTH:
		if (value == 0)
			return(COPY1);
		else if (value == 1)
			return(BOTH_PRESENT);
		else if (value == 2)
			return(BOTH_1_BIAS);
		else 
			return(BOTH_2_BIAS);
		break;
	case COPY2_OR_BOTH:
		if (value == 0)
			return(COPY2);
		else if (value == 1)
			return(BOTH_PRESENT);
		else if (value == 2)
			return(BOTH_1_BIAS);
		else 
			return(BOTH_2_BIAS);
		break;

	}


}



BOOL Dupl_Parallel_Base_model::allowed_state(DUPL_LOSS_STATES the_state)
{
	BOOL retval=(BOOL)FALSE;

	switch (the_state) {
	case BOTH_PRESENT:
	case BOTH_1_BIAS:
	case BOTH_2_BIAS:
	case COPY1:
	case COPY2:
	case COPY1_OR_BOTH:
	case COPY2_OR_BOTH:
	case GENERIC_SINGLE_COPY:
	case MISSING:
	case LOST:
		retval=(BOOL)TRUE;
		break;

	}
	return(retval);
}

double Dupl_Parallel_Base_model::find_ut (Branch *taxa)
{
	return(taxa->expect_subs_site()/(2.0+2.0*curr_exchange->get_dupl_parallel_rate()));
}



void Dupl_Parallel_2_Rate_Base_model::calc_transprobs(Branch *taxa, int rate_num)
{
		int i, j;
		double a,b,d;

		//Set transprobs so that unused states in this model are never entered
		for (i=0;i<curr_exchange->get_condlike_size(); i++){
			if (dupl_to_loss_state(i) != LOST) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(LOST), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(LOST), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(LOST),loss_state_to_dupl(LOST), 1.0);

			if (dupl_to_loss_state(i) != BOTH_PRESENT_FIXED) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(BOTH_PRESENT_FIXED), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED),
					loss_state_to_dupl(BOTH_PRESENT_FIXED), 1.0);

			if (dupl_to_loss_state(i) != GENERIC_SINGLE_COPY) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(GENERIC_SINGLE_COPY), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(GENERIC_SINGLE_COPY), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(GENERIC_SINGLE_COPY),
					loss_state_to_dupl(GENERIC_SINGLE_COPY), 1.0);

			if (dupl_to_loss_state(i) != MISSING) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(MISSING), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(MISSING), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(MISSING),
					loss_state_to_dupl(MISSING), 1.0);

			if (dupl_to_loss_state(i) != COPY1_OR_BOTH) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(COPY1_OR_BOTH), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1_OR_BOTH), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1_OR_BOTH),
					loss_state_to_dupl(COPY1_OR_BOTH), 1.0);

			if (dupl_to_loss_state(i) != COPY2_OR_BOTH) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(COPY2_OR_BOTH), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2_OR_BOTH), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2_OR_BOTH),
					loss_state_to_dupl(COPY2_OR_BOTH), 1.0);

				if (dupl_to_loss_state(i) != BOTH_FIXED_SUBF) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(BOTH_FIXED_SUBF), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_FIXED_SUBF), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_FIXED_SUBF),loss_state_to_dupl(BOTH_FIXED_SUBF), 1.0);

			if (dupl_to_loss_state(i) != COPY1_BIAS) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(COPY1_BIAS), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1_BIAS), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1_BIAS),loss_state_to_dupl(COPY1_BIAS), 1.0);

			if (dupl_to_loss_state(i) != COPY2_BIAS) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(COPY2_BIAS), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2_BIAS), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2_BIAS),loss_state_to_dupl(COPY2_BIAS), 1.0);
			

		}


		a= taxa->get_brnlen();
		b=curr_exchange->get_dupl_parallel_rate();
		d=curr_exchange->get_loss_rate_scale();


		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(COPY1), 1.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(COPY2), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_2_BIAS), 0.0);

		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(COPY2), 1.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(COPY1), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_2_BIAS), 0.0);

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_2_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(COPY1), 1.0-exp(-1.0*a*d));
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(COPY2), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_1_BIAS), exp(-1.0*a*d));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(COPY1), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(COPY2), 1.0-exp(-1.0*a*d));
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_2_BIAS), exp(-1.0*a*d));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(COPY1), 
			(2 + 2*b - (2*b)/exp(a*d) + (-2 + d)/exp(2*a*(1 + b)) - d)/(4 + 4*b - 2*d));
		
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(COPY2),
			(2 + 2*b - (2*b)/exp(a*d) + (-2 + d)/exp(2*a*(1 + b)) - d)/(4 + 4*b - 2*d));
		
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_1_BIAS),
				-(((exp(-2*a*(1 + b)) - exp(-(a*d)))*b)/(2 + 2*b - d)));
		
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_2_BIAS),
			-(((exp(-2*a*(1 + b)) - exp(-(a*d)))*b)/(2 + 2*b - d)));
		
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_PRESENT), 
			exp(-2*a*(1 + b)));

}




Dupl_Parallel_model::Dupl_Parallel_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree)
{
	if(cexchange->get_dupl_parallel_rate() ==0)
		cexchange->set_dupl_parallel_rate(0.1); 
	assemble (cexchange, cdata, ctree);
}

void Dupl_Parallel_model::describe_results()
{
	cout<<"Instantaneous relative rate of duplicate transition to parallel loss state: "
		<<curr_exchange->get_dupl_parallel_rate()<<"\n";
}



void Dupl_Parallel_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
	int i, brn_cnt=0;

	brn_start=0;
	for(i=0; i<curr_exchange->get_num_branches(); i++) {
		if ((curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE) || ((*curr_tree)[i]->get_brnlen() != 0)) {
				par[brn_cnt]=(*curr_tree)[i]->get_brnlen();
				types[brn_cnt]=BRANCH;
				if (curr_exchange->zero_len_brns_fixed() ==(BOOL)TRUE)
					brn_index[brn_cnt]=i;
				brn_cnt++;
		}
	}


	par[curr_exchange->get_num_params()-1]=curr_exchange->get_dupl_parallel_rate();
	types[curr_exchange->get_num_params()-1]=DUPL_PARALLEL_RATE;
}
    

void Dupl_Parallel_model::num_params_model() 
{
	int i, cnt_zero=0;

	if (curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE)
		curr_exchange->set_num_params(curr_exchange->get_num_branches()+1);
	else {
		for(i=0; i<curr_exchange->get_num_branches(); i++)
		{
			if ((*curr_tree)[i]->expect_subs_site() ==0)
				cnt_zero++;
		}
		brn_index=new int[curr_exchange->get_num_branches()-cnt_zero];
		curr_exchange->set_num_params(curr_exchange->get_num_branches()-cnt_zero+1);
	}
}




Dupl_Parallel_NoState_model::Dupl_Parallel_NoState_model (Exchange *cexchange, Tree *ctree, Clade *cgenomes, WGD_Data *chomologs) 
{
	allocate_state_model(cexchange, ctree, cgenomes, chomologs);
	
	if(cexchange->get_dupl_parallel_rate() ==0)
		cexchange->set_dupl_parallel_rate(0.1);

	left_cond_probs=right_cond_probs=0;

	assemble (cexchange, curr_data, ctree);
}

void Dupl_Parallel_NoState_model::describe_results()
{
	cout<<"Instantaneous relative rate of duplicate transition to parallel loss state: "
		<<curr_exchange->get_dupl_parallel_rate()<<"\n";
	cout<<"Probability of strand definition switching between two single copy genes: "
		<<curr_exchange->get_strand_switch_prob()<<endl;
}



void Dupl_Parallel_NoState_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
	int i, brn_cnt=0;

	brn_start=0;
	for(i=0; i<curr_exchange->get_num_branches(); i++) {
		if ((curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE) || ((*curr_tree)[i]->get_brnlen() != 0)) {
				par[brn_cnt]=(*curr_tree)[i]->get_brnlen();
				types[brn_cnt]=BRANCH;
				if (curr_exchange->zero_len_brns_fixed() ==(BOOL)TRUE)
					brn_index[brn_cnt]=i;
				brn_cnt++;
		}
	}


	par[curr_exchange->get_num_params()-2]=curr_exchange->get_dupl_parallel_rate();
	types[curr_exchange->get_num_params()-2]=DUPL_PARALLEL_RATE;

	par[curr_exchange->get_num_params()-1]=log(curr_exchange->get_strand_switch_prob());
	types[curr_exchange->get_num_params()-1]=STRAND_SWITCH_PROB;
}
    

void Dupl_Parallel_NoState_model::num_params_model() 
{
	int i, cnt_zero=0;

	if (curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE)
		curr_exchange->set_num_params(curr_exchange->get_num_branches()+2);
	else {
		for(i=0; i<curr_exchange->get_num_branches(); i++)
		{
			if ((*curr_tree)[i]->expect_subs_site() ==0)
				cnt_zero++;
		}
		brn_index=new int[curr_exchange->get_num_branches()-cnt_zero];
		curr_exchange->set_num_params(curr_exchange->get_num_branches()-cnt_zero+2);
	}
}

Dupl_Parallel_2_Rate_model::Dupl_Parallel_2_Rate_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree)
{
	if(cexchange->get_dupl_parallel_rate() ==0)
		cexchange->set_dupl_parallel_rate(0.1); 
	if(cexchange->get_loss_rate_scale() ==0)
		cexchange->set_loss_rate_scale(1.0);
	assemble (cexchange, cdata, ctree);
}

void Dupl_Parallel_2_Rate_model::describe_results()
{
	cout<<"Five state model with convergening states and two rates of loss\n";
	cout<<"Instantaneous relative rate of duplicate transition to parallel loss state: "
		<<curr_exchange->get_dupl_parallel_rate()<<"\n";
	cout<<"Relative rate of duplicate loss after convergence: "<<curr_exchange->get_loss_rate_scale()<<endl;
}



void Dupl_Parallel_2_Rate_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
	int i, brn_cnt=0;

	brn_start=0;
	for(i=0; i<curr_exchange->get_num_branches(); i++) {
		if ((curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE) || ((*curr_tree)[i]->get_brnlen() != 0)) {
				par[brn_cnt]=(*curr_tree)[i]->get_brnlen();
				types[brn_cnt]=BRANCH;
				if (curr_exchange->zero_len_brns_fixed() ==(BOOL)TRUE)
					brn_index[brn_cnt]=i;
				brn_cnt++;
		}
	}


	par[curr_exchange->get_num_params()-2]=curr_exchange->get_dupl_parallel_rate();
	types[curr_exchange->get_num_params()-2]=DUPL_PARALLEL_RATE;

	par[curr_exchange->get_num_params()-1]=curr_exchange->get_loss_rate_scale();
	types[curr_exchange->get_num_params()-1]=CON_REL_LOSS_RATE;
}
    

void Dupl_Parallel_2_Rate_model::num_params_model() 
{
	int i, cnt_zero=0;

	if (curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE)
		curr_exchange->set_num_params(curr_exchange->get_num_branches()+1);
	else {
		for(i=0; i<curr_exchange->get_num_branches(); i++)
		{
			if ((*curr_tree)[i]->expect_subs_site() ==0)
				cnt_zero++;
		}
		brn_index=new int[curr_exchange->get_num_branches()-cnt_zero];
		curr_exchange->set_num_params(curr_exchange->get_num_branches()-cnt_zero+2);
	}
}




Dupl_Parallel_2_Rate_NoState_model::Dupl_Parallel_2_Rate_NoState_model (Exchange *cexchange, Tree *ctree, Clade *cgenomes, WGD_Data *chomologs) 
{
	allocate_state_model(cexchange, ctree, cgenomes, chomologs);
	
	if(cexchange->get_dupl_parallel_rate() ==0)
		cexchange->set_dupl_parallel_rate(0.1);
	if(cexchange->get_loss_rate_scale() == 0)
		cexchange->set_loss_rate_scale(1.0);

	left_cond_probs=right_cond_probs=0;

	assemble (cexchange, curr_data, ctree);
}

void Dupl_Parallel_2_Rate_NoState_model::describe_results()
{
	cout<<"Five state model with convergening duplicates, two rates of gene loss and inferred tracking\n";
	cout<<"Instantaneous relative rate of duplicate transition to parallel loss state: "
		<<curr_exchange->get_dupl_parallel_rate()<<"\n";
	cout<<"Relative rate of gene loss after convergence: "<<curr_exchange->get_loss_rate_scale()<<endl;
	cout<<"Probability of strand definition switching between two single copy genes: "
		<<curr_exchange->get_strand_switch_prob()<<endl;
}



void Dupl_Parallel_2_Rate_NoState_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
	int i, brn_cnt=0;

	brn_start=0;
	for(i=0; i<curr_exchange->get_num_branches(); i++) {
		if ((curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE) || ((*curr_tree)[i]->get_brnlen() != 0)) {
				par[brn_cnt]=(*curr_tree)[i]->get_brnlen();
				types[brn_cnt]=BRANCH;
				if (curr_exchange->zero_len_brns_fixed() ==(BOOL)TRUE)
					brn_index[brn_cnt]=i;
				brn_cnt++;
		}
	}


	par[curr_exchange->get_num_params()-3]=curr_exchange->get_dupl_parallel_rate();
	types[curr_exchange->get_num_params()-3]=DUPL_PARALLEL_RATE;

	par[curr_exchange->get_num_params()-1]=curr_exchange->get_loss_rate_scale();
	types[curr_exchange->get_num_params()-1]=CON_REL_LOSS_RATE;

	par[curr_exchange->get_num_params()-2]=log(curr_exchange->get_strand_switch_prob());
	types[curr_exchange->get_num_params()-2]=STRAND_SWITCH_PROB;
}
    

void Dupl_Parallel_2_Rate_NoState_model::num_params_model() 
{
	int i, cnt_zero=0;

	if (curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE)
		curr_exchange->set_num_params(curr_exchange->get_num_branches()+2);
	else {
		for(i=0; i<curr_exchange->get_num_branches(); i++)
		{
			if ((*curr_tree)[i]->expect_subs_site() ==0)
				cnt_zero++;
		}
		brn_index=new int[curr_exchange->get_num_branches()-cnt_zero];
		curr_exchange->set_num_params(curr_exchange->get_num_branches()-cnt_zero+3);
	}
}


void Dupl_Fix_Parallel_Base_model::set_null_transprobs(Branch *taxa, int rate_num)
{
	int i;
	
	//Set transprobs so that unused states in this model are never entered
		for (i=0;i<curr_exchange->get_condlike_size(); i++){
			if (dupl_to_loss_state(i) != LOST) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(LOST), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(LOST), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(LOST),loss_state_to_dupl(LOST), 1.0);

					if (dupl_to_loss_state(i) != GENERIC_SINGLE_COPY) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(GENERIC_SINGLE_COPY), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(GENERIC_SINGLE_COPY), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(GENERIC_SINGLE_COPY),
					loss_state_to_dupl(GENERIC_SINGLE_COPY), 1.0);

			if (dupl_to_loss_state(i) != MISSING) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(MISSING), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(MISSING), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(MISSING),
					loss_state_to_dupl(MISSING), 1.0);

			if (dupl_to_loss_state(i) != COPY1_OR_BOTH) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(COPY1_OR_BOTH), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1_OR_BOTH), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1_OR_BOTH),
					loss_state_to_dupl(COPY1_OR_BOTH), 1.0);

			if (dupl_to_loss_state(i) != COPY2_OR_BOTH) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(COPY2_OR_BOTH), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2_OR_BOTH), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2_OR_BOTH),
					loss_state_to_dupl(COPY2_OR_BOTH), 1.0);
			
			if (dupl_to_loss_state(i) != BOTH_FIXED_SUBF) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(BOTH_FIXED_SUBF), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_FIXED_SUBF), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_FIXED_SUBF),loss_state_to_dupl(BOTH_FIXED_SUBF), 1.0);

			if (dupl_to_loss_state(i) != COPY1_BIAS) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(COPY1_BIAS), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1_BIAS), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1_BIAS),loss_state_to_dupl(COPY1_BIAS), 1.0);

			if (dupl_to_loss_state(i) != COPY2_BIAS) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(COPY2_BIAS), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2_BIAS), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2_BIAS),loss_state_to_dupl(COPY2_BIAS), 1.0);

		}

}

void Dupl_Fix_Parallel_Base_model::calc_transprobs(Branch *taxa, int rate_num)
{
		int i, j;
		double two_plus_two_beta_gamma, one_plus_two_beta_gamma,
			one_plus_beta, one_plus_beta_gamma;

		set_null_transprobs(taxa, rate_num);

		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(COPY1), 1.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(COPY2), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_PRESENT_FIXED), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_2_BIAS), 0.0);

		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(COPY2), 1.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(COPY1), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_PRESENT_FIXED), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_2_BIAS), 0.0);

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_PRESENT_FIXED), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_2_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(COPY1), 1.0-exp(-1.0*taxa->get_brnlen()));
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(COPY2), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_1_BIAS), exp(-1.0*taxa->get_brnlen()));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_PRESENT_FIXED), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(COPY1), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(COPY2), 1.0-exp(-1.0*taxa->get_brnlen()));
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_2_BIAS), exp(-1.0*taxa->get_brnlen()));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(COPY2), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(COPY1), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(BOTH_2_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(BOTH_PRESENT_FIXED), 1.0);


		two_plus_two_beta_gamma=2.0+curr_exchange->get_dupl_fix_rate()+2.0*curr_exchange->get_dupl_parallel_rate();
		
		one_plus_two_beta_gamma=1.0+curr_exchange->get_dupl_fix_rate()+2.0*curr_exchange->get_dupl_parallel_rate();
	
		one_plus_beta=1.0+curr_exchange->get_dupl_parallel_rate();
	
		one_plus_beta_gamma=1.0+curr_exchange->get_dupl_fix_rate()+curr_exchange->get_dupl_parallel_rate();
	

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(COPY1),
				 (one_plus_beta*one_plus_two_beta_gamma-
				 one_plus_beta_gamma*exp(-1.0*taxa->get_brnlen()*two_plus_two_beta_gamma)-
				 curr_exchange->get_dupl_parallel_rate()*two_plus_two_beta_gamma*exp(-1.0*taxa->get_brnlen()))
				 /(one_plus_two_beta_gamma*two_plus_two_beta_gamma));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(COPY2),
				(one_plus_beta*one_plus_two_beta_gamma-
				 one_plus_beta_gamma*exp(-1.0*taxa->get_brnlen()*two_plus_two_beta_gamma)-
				 curr_exchange->get_dupl_parallel_rate()*two_plus_two_beta_gamma*exp(-1.0*taxa->get_brnlen()))
				 /(one_plus_two_beta_gamma*two_plus_two_beta_gamma));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_1_BIAS),
			(curr_exchange->get_dupl_parallel_rate()/one_plus_two_beta_gamma)*
			(exp(-1.0*taxa->get_brnlen())-exp(-1.0*taxa->get_brnlen()*two_plus_two_beta_gamma)));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_2_BIAS),
			(curr_exchange->get_dupl_parallel_rate()/one_plus_two_beta_gamma)*
			(exp(-1.0*taxa->get_brnlen())-exp(-1.0*taxa->get_brnlen()*two_plus_two_beta_gamma)));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_PRESENT_FIXED),
				 (curr_exchange->get_dupl_fix_rate()/two_plus_two_beta_gamma)*
				 (1.0-exp(-1.0*taxa->get_brnlen()*two_plus_two_beta_gamma)));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_PRESENT), 
			exp(-1.0*taxa->get_brnlen()*two_plus_two_beta_gamma));



}


int Dupl_Fix_Parallel_Base_model::state_redundancy_val(DUPL_LOSS_STATES the_state)
{
		//Tells us if we need to sum over internal model states for a given input tip state

	if (observable_dupl_state(the_state)==(BOOL)TRUE) {
		switch (the_state) {
		case COPY1:
		case COPY2:
		case LOST:
			return(1);
			break;
		case BOTH_PRESENT:
			return(4);
			break;
		case MISSING:
			return(6);
			break;
		case GENERIC_SINGLE_COPY:
			return(2);
			break;
		case COPY1_OR_BOTH:
		case COPY2_OR_BOTH:
			return(5);
			break;
		}
	}
	else {
		cerr<<"ERROR: Tried to get Tip redundancy for non-observable internal model state\n";
		return(1);
	}

}



DUPL_LOSS_STATES Dupl_Fix_Parallel_Base_model::get_redund_position_n(DUPL_LOSS_STATES the_state, int value)
{
		//Indicates order of the states to be summed through for ambiguous observed states

	switch(the_state) {
	
	//If it is not an ambiguous state, just return it
	case COPY1:
	case COPY2:
	case LOST:
		return(the_state);
		break;

	//Otherwise, give the possible states in an arbitrary order
	case BOTH_PRESENT:
		if (value == 0)
			return(BOTH_PRESENT);
		else if (value == 1)
			return(BOTH_1_BIAS);
		else if (value == 2)
			return(BOTH_2_BIAS);
		else
			return(BOTH_PRESENT_FIXED);
		break;
	case MISSING:
		if (value == 0)
			return(BOTH_PRESENT);
		else if (value == 1)
			return(COPY1);
		else if (value == 2)
			return(COPY2);
		else if (value == 3)
			return(BOTH_1_BIAS);
		else if (value == 4)
			return(BOTH_2_BIAS);
		else 
			return(BOTH_PRESENT_FIXED);
		break;
	case GENERIC_SINGLE_COPY:
		if (value == 0)
			return(COPY1);
		else
			return(COPY2);
		break;
	case COPY1_OR_BOTH:
		if (value == 0)
			return(COPY1);
		else if (value == 1)
			return(BOTH_PRESENT);
		else if (value == 2)
			return(BOTH_PRESENT_FIXED);
		else if (value == 3)
			return(BOTH_1_BIAS);
		else 
			return(BOTH_2_BIAS);
		break;
	case COPY2_OR_BOTH:
		if (value == 0)
			return(COPY2);
		else if (value == 1)
			return(BOTH_PRESENT);
		else if (value == 2)
			return(BOTH_PRESENT_FIXED);
		else if (value == 3)
			return(BOTH_1_BIAS);
		else 
			return(BOTH_2_BIAS);
		break;


	}


}


BOOL Dupl_Fix_Parallel_Base_model::allowed_state(DUPL_LOSS_STATES the_state)
{
	BOOL retval=(BOOL)FALSE;

	switch (the_state) {
	case BOTH_PRESENT:
	case BOTH_1_BIAS:
	case BOTH_2_BIAS:
	case BOTH_PRESENT_FIXED:
	case COPY1:
	case COPY2:
	case COPY1_OR_BOTH:
	case COPY2_OR_BOTH:
	case GENERIC_SINGLE_COPY:
	case MISSING:
	case LOST:
		retval=(BOOL)TRUE;
		break;

	}
	return(retval);
}


double Dupl_Fix_Parallel_Base_model::find_ut (Branch *taxa)
{
	return(taxa->expect_subs_site()/(2.0+2.0*curr_exchange->get_dupl_parallel_rate()+curr_exchange->get_dupl_fix_rate()));
}




void Dupl_Fix_Parallel_SubF_Base_model::calc_transprobs(Branch *taxa, int rate_num)
{
		int i, j;
		double two_plus_two_beta_gamma, one_plus_gamma,
			one_plus_beta, one_plus_two_beta, one_plus_beta_gamma, one_plus_two_beta_gamma, d_p_e;

		set_null_transprobs(taxa, rate_num);

		two_plus_two_beta_gamma=2.0+curr_exchange->get_dupl_fix_rate()+2.0*curr_exchange->get_dupl_parallel_rate();
		
		one_plus_two_beta=1.0+2.0*curr_exchange->get_dupl_parallel_rate();
	
		one_plus_beta=1.0+curr_exchange->get_dupl_parallel_rate();
	
		one_plus_gamma=1.0+curr_exchange->get_dupl_fix_rate();
		
		one_plus_beta_gamma=1.0+curr_exchange->get_dupl_parallel_rate()+curr_exchange->get_dupl_fix_rate();
	
		one_plus_two_beta_gamma=1.0+2.0*curr_exchange->get_dupl_parallel_rate()+curr_exchange->get_dupl_fix_rate();


		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(COPY1), 1.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(COPY2), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_PRESENT_FIXED), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_2_BIAS), 0.0);

		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(COPY2), 1.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(COPY1), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_PRESENT_FIXED), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_2_BIAS), 0.0);

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_PRESENT_FIXED), 
			curr_exchange->get_dupl_fix_rate()*(1.0-exp(-1.0*taxa->get_brnlen()*one_plus_gamma))/one_plus_gamma);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_2_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(COPY1), 
			(1.0-exp(-1.0*taxa->get_brnlen()*one_plus_gamma))/one_plus_gamma);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(COPY2), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_1_BIAS), 
				exp(-1.0*taxa->get_brnlen()*one_plus_gamma));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_PRESENT_FIXED),
			curr_exchange->get_dupl_fix_rate()*(1.0-exp(-1.0*taxa->get_brnlen()*one_plus_gamma))/one_plus_gamma);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(COPY1), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(COPY2), 
			(1.0-exp(-1.0*taxa->get_brnlen()*one_plus_gamma))/one_plus_gamma);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_2_BIAS), 
				exp(-1.0*taxa->get_brnlen()*one_plus_gamma));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(COPY2), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(COPY1), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(BOTH_2_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(BOTH_PRESENT_FIXED), 1.0);


		
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(COPY1),(one_plus_two_beta*one_plus_beta_gamma-
			one_plus_beta*one_plus_gamma*exp(-1.0*taxa->get_brnlen()*two_plus_two_beta_gamma)-
			curr_exchange->get_dupl_parallel_rate()*two_plus_two_beta_gamma*exp(-1.0*one_plus_gamma*taxa->get_brnlen()))/
			(one_plus_two_beta*one_plus_gamma*two_plus_two_beta_gamma)  );

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(COPY2),(one_plus_two_beta*one_plus_beta_gamma-
			one_plus_beta*one_plus_gamma*exp(-1.0*taxa->get_brnlen()*two_plus_two_beta_gamma)-
			curr_exchange->get_dupl_parallel_rate()*two_plus_two_beta_gamma*exp(-1.0*one_plus_gamma*taxa->get_brnlen()))/
			(one_plus_two_beta*one_plus_gamma*two_plus_two_beta_gamma) 	);

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_1_BIAS),curr_exchange->get_dupl_parallel_rate()*
			(exp(-1.0*one_plus_gamma*taxa->get_brnlen())-exp(-1.0*taxa->get_brnlen()*two_plus_two_beta_gamma)) / (one_plus_two_beta) );

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_2_BIAS), curr_exchange->get_dupl_parallel_rate()*
			(exp(-1.0*one_plus_gamma*taxa->get_brnlen())-exp(-1.0*taxa->get_brnlen()*two_plus_two_beta_gamma)) / (one_plus_two_beta)	);

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_PRESENT_FIXED), 
			curr_exchange->get_dupl_fix_rate()*
			(one_plus_two_beta_gamma*one_plus_two_beta-one_plus_gamma*exp(-1.0*taxa->get_brnlen()*two_plus_two_beta_gamma)-
			2.0*curr_exchange->get_dupl_parallel_rate()*two_plus_two_beta_gamma*exp(-1.0*taxa->get_brnlen()*one_plus_gamma))/ 
			(one_plus_two_beta*one_plus_gamma*two_plus_two_beta_gamma)	 );

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_PRESENT), 
			exp(-1.0*two_plus_two_beta_gamma*taxa->get_brnlen()));



}



void Dupl_2_Rate_NoSubF_Base_model::calc_transprobs(Branch *taxa, int rate_num)
{
		int i, j;
		double a, b, g, d;

		set_null_transprobs(taxa, rate_num);


		g=curr_exchange->get_dupl_fix_rate();
		b=curr_exchange->get_dupl_parallel_rate();
		a=taxa->get_brnlen();
		d=curr_exchange->get_loss_rate_scale();
	
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(COPY1), 1.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(COPY2), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_PRESENT_FIXED), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_2_BIAS), 0.0);

		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(COPY2), 1.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(COPY1), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_PRESENT_FIXED), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_2_BIAS), 0.0);

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_PRESENT_FIXED), 
			0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_2_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(COPY1), 
			1.0 - exp(-1.0*a*d));
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(COPY2), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_1_BIAS), 
			exp(-1.0*a*d));



		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_PRESENT_FIXED),
			0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(COPY1), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(COPY2), 
			1.0 - exp(-1.0*a*d));
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_2_BIAS), 
			exp(-1.0*a*d));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(COPY2), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(COPY1), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(BOTH_2_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(BOTH_PRESENT_FIXED), 1.0);


		
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(COPY1),			
			-((b*(2 + 2*b + g) - exp(a*d)*(1 + b)*(2 + 2*b + g - d) + (2 + 2*b + g - (1 + b)*d)/exp(a*(2 + 2*b + g - d)))/
        (exp(a*d)*(2 + 2*b + g)*(2 + 2*b + g - d))));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(COPY2),
			-((b*(2 + 2*b + g) - exp(a*d)*(1 + b)*(2 + 2*b + g - d) + (2 + 2*b + g - (1 + b)*d)/exp(a*(2 + 2*b + g - d)))/
        (exp(a*d)*(2 + 2*b + g)*(2 + 2*b + g - d))));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_1_BIAS),
			(b - b/exp(a*(2 + 2*b + g - d)))/(exp(a*d)*(2 + 2*b + g - d)));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_2_BIAS),
			(b - b/exp(a*(2 + 2*b + g - d)))/(exp(a*d)*(2 + 2*b + g - d)));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_PRESENT_FIXED),
			(g - g/exp(a*(2 + 2*b + g)))/(2 + 2*b + g));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_PRESENT), 
			exp(-(a*(2 + 2*b + g))));



}



void Dupl_2_Rate_NoSubF_AllStates_Base_model::calc_transprobs(Branch *taxa, int rate_num)
{
		int i, j;
		double a, b, g, d;

		//Set transprobs so that unused states in this model are never entered
		for (i=0;i<curr_exchange->get_condlike_size(); i++){
			if (dupl_to_loss_state(i) != LOST) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(LOST), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(LOST), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(LOST),loss_state_to_dupl(LOST), 1.0);

					if (dupl_to_loss_state(i) != GENERIC_SINGLE_COPY) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(GENERIC_SINGLE_COPY), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(GENERIC_SINGLE_COPY), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(GENERIC_SINGLE_COPY),
					loss_state_to_dupl(GENERIC_SINGLE_COPY), 1.0);

			if (dupl_to_loss_state(i) != MISSING) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(MISSING), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(MISSING), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(MISSING),
					loss_state_to_dupl(MISSING), 1.0);

			if (dupl_to_loss_state(i) != COPY1_OR_BOTH) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(COPY1_OR_BOTH), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1_OR_BOTH), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1_OR_BOTH),
					loss_state_to_dupl(COPY1_OR_BOTH), 1.0);

			if (dupl_to_loss_state(i) != COPY2_OR_BOTH) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(COPY2_OR_BOTH), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2_OR_BOTH), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2_OR_BOTH),
					loss_state_to_dupl(COPY2_OR_BOTH), 1.0);

				if (dupl_to_loss_state(i) != BOTH_FIXED_SUBF) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(BOTH_FIXED_SUBF), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_FIXED_SUBF), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_FIXED_SUBF),loss_state_to_dupl(BOTH_FIXED_SUBF), 1.0);

			if (dupl_to_loss_state(i) != COPY1_BIAS) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(COPY1_BIAS), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1_BIAS), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1_BIAS),loss_state_to_dupl(COPY1_BIAS), 1.0);

			if (dupl_to_loss_state(i) != COPY2_BIAS) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(COPY2_BIAS), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2_BIAS), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2_BIAS),loss_state_to_dupl(COPY2_BIAS), 1.0);
			

		}


		g=curr_exchange->get_dupl_fix_rate();
		b=curr_exchange->get_dupl_parallel_rate();
		a=taxa->get_brnlen();
		d=curr_exchange->get_loss_rate_scale();
	
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(COPY1), 1.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(COPY2), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(COPY1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(COPY2_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_PRESENT_FIXED), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_2_BIAS), 0.0);

		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(COPY2), 1.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(COPY1), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(COPY1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(COPY2_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_PRESENT_FIXED), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_2_BIAS), 0.0);

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_PRESENT_FIXED), 
			0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_2_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(COPY1_BIAS), 
			1.0 - exp(-1.0*a*d));
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(COPY2), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(COPY1), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(COPY2_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_1_BIAS), 
			exp(-1.0*a*d));



		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_PRESENT_FIXED),
			0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(COPY1), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(COPY2_BIAS), 
			1.0 - exp(-1.0*a*d));
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(COPY2), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(COPY1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_2_BIAS), 
			exp(-1.0*a*d));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(COPY2), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(COPY1), 0.0);	
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(COPY1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(COPY2_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(BOTH_2_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(BOTH_PRESENT_FIXED), 1.0);


		
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(COPY1),
			(1 - exp(-(a*(2 + 2*b + g))))/(2 + 2*b + g));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(COPY1_BIAS),
			(b*(-2 - 2*b - g + exp(a*d)*(2 + 2*b + g - d) + d/exp(a*(2 + 2*b + g - d))))/
			(exp(a*d)*(2 + 2*b + g)*(2 + 2*b + g - d)));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(COPY2),
			(1 - exp(-(a*(2 + 2*b + g))))/(2 + 2*b + g));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(COPY2_BIAS),
			(b*(-2 - 2*b - g + exp(a*d)*(2 + 2*b + g - d) + d/exp(a*(2 + 2*b + g - d))))/
			(exp(a*d)*(2 + 2*b + g)*(2 + 2*b + g - d)));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_1_BIAS),
		(b - b/exp(a*(2 + 2*b + g - d)))/(exp(a*d)*(2 + 2*b + g - d)));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_2_BIAS),
			(b - b/exp(a*(2 + 2*b + g - d)))/(exp(a*d)*(2 + 2*b + g - d)));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_PRESENT_FIXED),
			(g - g/exp(a*(2 + 2*b + g)))/(2 + 2*b + g));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_PRESENT), 
			exp(-(a*(2 + 2*b + g))));



}


BOOL Dupl_2_Rate_NoSubF_AllStates_Base_model::allowed_state(DUPL_LOSS_STATES the_state)
{
	BOOL retval=(BOOL)FALSE;

	switch (the_state) {
	case BOTH_PRESENT:
	case BOTH_1_BIAS:
	case BOTH_2_BIAS:
	case BOTH_PRESENT_FIXED:
	case COPY1:
	case COPY2:
	case COPY1_BIAS:
	case COPY2_BIAS:
	case COPY1_OR_BOTH:
	case COPY2_OR_BOTH:
	case GENERIC_SINGLE_COPY:
	case MISSING:
	case LOST:
		retval=(BOOL)TRUE;
		break;

	}
	return(retval);
}


void Dupl_SubF_3_Rate_Base_model::calc_transprobs(Branch *taxa, int rate_num)
{
		int i, j;
		double a, b, g, d, e;

		set_null_transprobs(taxa, rate_num);


		g=curr_exchange->get_dupl_fix_rate();
		b=curr_exchange->get_dupl_parallel_rate();
		a=taxa->get_brnlen();
		d=curr_exchange->get_loss_rate_scale();
		e=curr_exchange->get_fix_rate_scale();

		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(COPY1), 1.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(COPY2), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_PRESENT_FIXED), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_2_BIAS), 0.0);

		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(COPY2), 1.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(COPY1), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_PRESENT_FIXED), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_2_BIAS), 0.0);

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_PRESENT_FIXED), 
			e*(1.0 - exp(-1.0*a*(d+e)))/(d+e));
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_2_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(COPY1), 
			d*(1.0 - exp(-1.0*a*(d+e)))/(d+e));
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(COPY2), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_1_BIAS), 
			exp(-1.0*a*(d+e)));



		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_PRESENT_FIXED),
			e*(1.0 - exp(-1.0*a*(d+e)))/(d+e));
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(COPY1), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(COPY2), 
			d*(1.0-  exp(-1.0*a*(d+e)))/(d+e));
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_2_BIAS), 
			exp(-1.0*a*(d+e)));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(COPY2), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(COPY1), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(BOTH_2_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(BOTH_PRESENT_FIXED), 1.0);


		
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(COPY1),
			(-2*e + (exp(2*a*(2 + 2*b + g))*b*d*(2 + 2*b + g) - exp(a*(2 + 2*b + g + d + e))*(d + e)*
			(-2 - g + b*(-2 + d) + d + e) + 
           exp(a*(4 + 4*b + 2*g + d + e))*(-((1 + b)*(2 + 2*b + g - d)*d) + (-g + b*(-2 + d) + 2*d)*e + (e*e)))/
         exp(a*(4 + 4*b + 2*g + d + e)))/((2 + 2*b + g)*(d + e)*(-2 - 2*b - g + d + e)));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(COPY2),
			(-2*e + (exp(2*a*(2 + 2*b + g))*b*(2 + 2*b + g)*d - exp(a*(2 + 2*b + g + d + e))*(d + e)*
			(-2 - g + b*(-2 + d) + d + e) + 
           exp(a*(4 + 4*b + 2*g + d + e))*(-((1 + b)*(2 + 2*b + g - d)*d) + (-g + b*(-2 + d) + 2*d)*e + (e*e)))/
         exp(a*(4 + 4*b + 2*g + d + e)))/((2 + 2*b + g)*(d + e)*(-2 - 2*b - g + d + e)));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_1_BIAS),
			b*(exp(-1.0*a*(d + e))-exp(-1.0*a*(2 + 2*b + g)))/(2 + 2*b + g - d - e));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_2_BIAS),
			b*(exp(-1.0*a*(d + e))-exp(-1.0*a*(2 + 2*b + g)))/(2 + 2*b + g - d - e));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_PRESENT_FIXED),
			(2*exp(2*a*(2 + 2*b + g))*b*(2 + 2*b + g)*e + 
        exp(a*(2 + 2*b + g + d + e))*(d + e)*(g*(2 + 2*b + g - d) - (2*b + g)*e) - 
        exp(a*(4 + 4*b + 2*g + d + e))*(2 + 2*b + g - d - e)*(2*b*e + g*(d + e)))/
      (exp(a*(4 + 4*b + 2*g + d + e))*(2 + 2*b + g)*(d + e)*(-2 - 2*b - g + d + e)));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_PRESENT), 
			exp(-1.0*a*(2+2*b+g)));



}




void Dupl_SubF_3_Rate_AllStates_Base_model::calc_transprobs(Branch *taxa, int rate_num)
{
		int i, j;
		double a, b, g, d, e;

		//Set transprobs so that unused states in this model are never entered
		for (i=0;i<curr_exchange->get_condlike_size(); i++){
			if (dupl_to_loss_state(i) != LOST) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(LOST), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(LOST), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(LOST),loss_state_to_dupl(LOST), 1.0);

					if (dupl_to_loss_state(i) != GENERIC_SINGLE_COPY) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(GENERIC_SINGLE_COPY), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(GENERIC_SINGLE_COPY), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(GENERIC_SINGLE_COPY),
					loss_state_to_dupl(GENERIC_SINGLE_COPY), 1.0);

			if (dupl_to_loss_state(i) != MISSING) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(MISSING), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(MISSING), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(MISSING),
					loss_state_to_dupl(MISSING), 1.0);

			if (dupl_to_loss_state(i) != COPY1_OR_BOTH) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(COPY1_OR_BOTH), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1_OR_BOTH), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1_OR_BOTH),
					loss_state_to_dupl(COPY1_OR_BOTH), 1.0);

			if (dupl_to_loss_state(i) != COPY2_OR_BOTH) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(COPY2_OR_BOTH), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2_OR_BOTH), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2_OR_BOTH),
					loss_state_to_dupl(COPY2_OR_BOTH), 1.0);

				if (dupl_to_loss_state(i) != BOTH_FIXED_SUBF) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(BOTH_FIXED_SUBF), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_FIXED_SUBF), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_FIXED_SUBF),loss_state_to_dupl(BOTH_FIXED_SUBF), 1.0);

			if (dupl_to_loss_state(i) != COPY1_BIAS) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(COPY1_BIAS), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1_BIAS), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1_BIAS),loss_state_to_dupl(COPY1_BIAS), 1.0);

			if (dupl_to_loss_state(i) != COPY2_BIAS) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(COPY2_BIAS), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2_BIAS), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2_BIAS),loss_state_to_dupl(COPY2_BIAS), 1.0);
			

		}


		g=curr_exchange->get_dupl_fix_rate();
		b=curr_exchange->get_dupl_parallel_rate();
		a=taxa->get_brnlen();
		d=curr_exchange->get_loss_rate_scale();
		e=curr_exchange->get_fix_rate_scale();

		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(COPY1), 1.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(COPY2), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(COPY1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(COPY2_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_PRESENT_FIXED), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_2_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_FIXED_SUBF), 0.0);

		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(COPY2), 1.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(COPY1), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(COPY1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(COPY2_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_PRESENT_FIXED), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_2_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_FIXED_SUBF), 0.0);


		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_PRESENT_FIXED), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_FIXED_SUBF), 
			e*(1.0 - exp(-1.0*a*(d+e)))/(d+e));
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_2_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(COPY1_BIAS), 
			d*(1.0 - exp(-1.0*a*(d+e)))/(d+e));
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(COPY1), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(COPY2_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(COPY2), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_1_BIAS), 
			exp(-1.0*a*(d+e)));



		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_PRESENT_FIXED), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_FIXED_SUBF),
			e*(1.0 - exp(-1.0*a*(d+e)))/(d+e));
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(COPY1), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(COPY2_BIAS), 
			d*(1.0-  exp(-1.0*a*(d+e)))/(d+e));
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(COPY2), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(COPY1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_2_BIAS), 
			exp(-1.0*a*(d+e)));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(COPY2), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(COPY1), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(COPY2_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(COPY1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(BOTH_2_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(BOTH_FIXED_SUBF), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(BOTH_PRESENT_FIXED), 1.0);


		
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(COPY1),
			(1 - exp(-(a*(2 + 2*b + g))))/(2 + 2*b + g));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(COPY2),
			(1 - exp(-(a*(2 + 2*b + g))))/(2 + 2*b + g));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(COPY1_BIAS), 
			-((b*d*(-2 - 2*b - g + (2 + 2*b + g)/exp(a*(d + e)) + d + e - (d + e)/exp(a*(2 + 2*b + g))))/
        ((2 + 2*b + g)*(2 + 2*b + g - d - e)*(d + e))));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(COPY2_BIAS), 
			-((b*d*(-2 - 2*b - g + (2 + 2*b + g)/exp(a*(d + e)) + d + e - (d + e)/exp(a*(2 + 2*b + g))))/
        ((2 + 2*b + g)*(2 + 2*b + g - d - e)*(d + e))));


		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_1_BIAS),
			((-exp(-(a*(2 + 2*b + g))) + exp(-(a*(d + e))))*b)/(2 + 2*b + g - d - e)	);

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_2_BIAS),
			((-exp(-(a*(2 + 2*b + g))) + exp(-(a*(d + e))))*b)/(2 + 2*b + g - d - e) );

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_PRESENT_FIXED), 
			(g - g/exp(a*(2 + 2*b + g)))/(2 + 2*b + g));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_FIXED_SUBF),
			(-2*b*e*(-2 - 2*b - g + (2 + 2*b + g)/exp(a*(d + e)) + d + e - (d + e)/exp(a*(2 + 2*b + g))))/
      ((2 + 2*b + g)*(2 + 2*b + g - d - e)*(d + e)));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_PRESENT), 
			exp(-(a*(2 + 2*b + g))));



}






void Dupl_SubF_Only_Base_model::calc_transprobs(Branch *taxa, int rate_num)
{
		int i, j;
		double a, b, g, d;

		set_null_transprobs(taxa, rate_num);


		g=curr_exchange->get_dupl_fix_rate();
		b=curr_exchange->get_dupl_parallel_rate();
		a=taxa->get_brnlen();
		d=curr_exchange->get_loss_rate_scale();

		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(COPY1), 1.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(COPY2), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_PRESENT_FIXED), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_2_BIAS), 0.0);

		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(COPY2), 1.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(COPY1), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_PRESENT_FIXED), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_2_BIAS), 0.0);

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_PRESENT_FIXED),
			g*(1.0 - exp(-1.0*a*(d+g)))/(d+g));
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_2_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(COPY1),
			d*(1.0 - exp(-1.0*a*(d+g)))/(d+g));
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(COPY2), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_1_BIAS), 
				exp(-1.0*a*(d+g)));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_PRESENT_FIXED),
			g*(1.0 - exp(-1.0*a*(d+g)))/(d+g));
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(COPY1), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(COPY2), 
			d*(1.0 - exp(-1.0*a*(d+g)))/(d+g));
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_2_BIAS), 
				exp(-1.0*a*(d+g)));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(COPY2), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(COPY1), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(BOTH_2_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(BOTH_PRESENT_FIXED), 1.0);


		
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(COPY1),
			 -(-((-1.0 + exp(2.0*a*(1 + b)))*(g*g)) + (1.0 + b)*d*
          (-2.0 - 2.0*exp(a*(2.0 + 2.0*b - g - d))*b + exp(2.0*a*(1 + b))*(2.0 + 2*b - d) + d) - 
         (-1 + exp(2*a*(1 + b)))*g*(-2*(1 + b) + (2 + b)*d))/(2.0*exp(2*a*(1.0 + b))*(1.0 + b)*(g + d)*(-2 - 2*b + g + d)));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(COPY2), 
			-(-((-1 + exp(2*a*(1 + b)))*(g*g)) + (1 + b)*d*
          (-2 - 2*exp(a*(2 + 2*b - g - d))*b + exp(2*a*(1 + b))*(2 + 2*b - d) + d) - 
         (-1 + exp(2*a*(1 + b)))*g*(-2*(1 + b) + (2 + b)*d))/(2.0*exp(2*a*(1 + b))*(1 + b)*(g + d)*(-2 - 2*b + g + d)));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_1_BIAS),
			-(((-1 + exp(a*(2 + 2*b - g - d)))*b)/(exp(2*a*(1 + b))*(-2 - 2*b + g + d))));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_2_BIAS),
			-(((-1 + exp(a*(2 + 2*b - g - d)))*b)/(exp(2*a*(1 + b))*(-2 - 2*b + g + d))));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_PRESENT_FIXED),
			(b*g*(2*exp(a*(2 + 2*b - g - d))*(1 + b) - g - d + exp(2*a*(1 + b))*(-2 - 2*b + g + d)))/
			(exp(2*a*(1 + b))*(1 + b)*(g + d)*(-2 - 2*b + g + d)));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_PRESENT), 
			exp(-2.0*a*(1 + b)));



}





void Dupl_Subf_All_States_Base_model::calc_transprobs(Branch *taxa, int rate_num)
{
		int i, j;
		double two_plus_two_beta_gamma, one_plus_two_beta_gamma,
			one_plus_beta, one_plus_beta_gamma, one_plus_gamma, one_plus_two_beta;

		//Set transprobs so that unused states in this model are never entered
		for (i=0;i<curr_exchange->get_condlike_size(); i++){
			if (dupl_to_loss_state(i) != LOST) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(LOST), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(LOST), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(LOST),loss_state_to_dupl(LOST), 1.0);

					if (dupl_to_loss_state(i) != GENERIC_SINGLE_COPY) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(GENERIC_SINGLE_COPY), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(GENERIC_SINGLE_COPY), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(GENERIC_SINGLE_COPY),
					loss_state_to_dupl(GENERIC_SINGLE_COPY), 1.0);

			if (dupl_to_loss_state(i) != MISSING) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(MISSING), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(MISSING), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(MISSING),
					loss_state_to_dupl(MISSING), 1.0);

			if (dupl_to_loss_state(i) != COPY1_OR_BOTH) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(COPY1_OR_BOTH), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1_OR_BOTH), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1_OR_BOTH),
					loss_state_to_dupl(COPY1_OR_BOTH), 1.0);

			if (dupl_to_loss_state(i) != COPY2_OR_BOTH) {
				taxa->set_trpb(rate_num, i,loss_state_to_dupl(COPY2_OR_BOTH), 0.0);
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2_OR_BOTH), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2_OR_BOTH),
					loss_state_to_dupl(COPY2_OR_BOTH), 1.0);
			

		}

		two_plus_two_beta_gamma=2.0+curr_exchange->get_dupl_fix_rate()+2.0*curr_exchange->get_dupl_parallel_rate();
		
		one_plus_two_beta_gamma=1.0+curr_exchange->get_dupl_fix_rate()+2.0*curr_exchange->get_dupl_parallel_rate();
	
		one_plus_beta=1.0+curr_exchange->get_dupl_parallel_rate();
	
		one_plus_beta_gamma=1.0+curr_exchange->get_dupl_fix_rate()+curr_exchange->get_dupl_parallel_rate();
	
		one_plus_gamma=1.0+curr_exchange->get_dupl_fix_rate();

		one_plus_two_beta=1.0+2.0*curr_exchange->get_dupl_parallel_rate();
	

		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(COPY1), 1.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(COPY2), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(COPY1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(COPY2_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_PRESENT_FIXED), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_2_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_FIXED_SUBF), 0.0);


		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(COPY2), 1.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(COPY1), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(COPY1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(COPY2_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_PRESENT_FIXED), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_2_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_FIXED_SUBF), 0.0);


		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_PRESENT_FIXED), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_2_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(COPY1), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(COPY1_BIAS), 
			(1.0-exp(-1.0*taxa->get_brnlen()*one_plus_gamma))/one_plus_gamma);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(COPY2_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(COPY2), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_1_BIAS), 
			exp(-1.0*taxa->get_brnlen()*one_plus_gamma));
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_FIXED_SUBF), 
			curr_exchange->get_dupl_fix_rate()*(1.0-exp(-1.0*taxa->get_brnlen()*one_plus_gamma))/one_plus_gamma);

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_PRESENT_FIXED), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(COPY1), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(COPY2), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(COPY1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(COPY2_BIAS), 
			(1.0-exp(-1.0*taxa->get_brnlen()*one_plus_gamma))/one_plus_gamma);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_2_BIAS), 
			exp(-1.0*taxa->get_brnlen()*one_plus_gamma));		
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_FIXED_SUBF), 
			curr_exchange->get_dupl_fix_rate()*(1.0-exp(-1.0*taxa->get_brnlen()*one_plus_gamma))/one_plus_gamma);

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(COPY2), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(COPY1), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(COPY1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(COPY2_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(BOTH_2_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(BOTH_PRESENT_FIXED), 1.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(BOTH_FIXED_SUBF), 0.0);

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_FIXED_SUBF), loss_state_to_dupl(COPY2), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_FIXED_SUBF), loss_state_to_dupl(COPY1), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_FIXED_SUBF), loss_state_to_dupl(COPY1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_FIXED_SUBF), loss_state_to_dupl(COPY2_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_FIXED_SUBF), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_FIXED_SUBF), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_FIXED_SUBF), loss_state_to_dupl(BOTH_2_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_FIXED_SUBF), loss_state_to_dupl(BOTH_PRESENT_FIXED), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_FIXED_SUBF), loss_state_to_dupl(BOTH_FIXED_SUBF), 1.0);

		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1_BIAS), loss_state_to_dupl(COPY2), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1_BIAS), loss_state_to_dupl(COPY1), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1_BIAS), loss_state_to_dupl(COPY1_BIAS), 1.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1_BIAS), loss_state_to_dupl(COPY2_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1_BIAS), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1_BIAS), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1_BIAS), loss_state_to_dupl(BOTH_2_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1_BIAS), loss_state_to_dupl(BOTH_PRESENT_FIXED), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1_BIAS), loss_state_to_dupl(BOTH_FIXED_SUBF), 0.0);

		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2_BIAS), loss_state_to_dupl(COPY2), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2_BIAS), loss_state_to_dupl(COPY1), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2_BIAS), loss_state_to_dupl(COPY1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2_BIAS), loss_state_to_dupl(COPY2_BIAS), 1.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2_BIAS), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2_BIAS), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2_BIAS), loss_state_to_dupl(BOTH_2_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2_BIAS), loss_state_to_dupl(BOTH_PRESENT_FIXED), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2_BIAS), loss_state_to_dupl(BOTH_FIXED_SUBF), 0.0);

	

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(COPY1),
				 (1.0-exp(-1.0*taxa->get_brnlen()*two_plus_two_beta_gamma))/two_plus_two_beta_gamma);

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(COPY2),
				 (1.0-exp(-1.0*taxa->get_brnlen()*two_plus_two_beta_gamma))/two_plus_two_beta_gamma);

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(COPY1_BIAS),
				 curr_exchange->get_dupl_parallel_rate()*
				 (one_plus_two_beta+one_plus_gamma*exp(-1.0*two_plus_two_beta_gamma*taxa->get_brnlen())-
				 two_plus_two_beta_gamma*exp(-1.0*one_plus_gamma*taxa->get_brnlen()))/
				 (one_plus_two_beta*one_plus_gamma*two_plus_two_beta_gamma));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(COPY2_BIAS),
				curr_exchange->get_dupl_parallel_rate()*
				 (one_plus_two_beta+one_plus_gamma*exp(-1.0*two_plus_two_beta_gamma*taxa->get_brnlen())-
				 two_plus_two_beta_gamma*exp(-1.0*one_plus_gamma*taxa->get_brnlen()))/
				 (one_plus_two_beta*one_plus_gamma*two_plus_two_beta_gamma));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_1_BIAS),
			curr_exchange->get_dupl_parallel_rate()* 
			(exp(-1.0*one_plus_gamma*taxa->get_brnlen())-exp(-1.0*two_plus_two_beta_gamma*taxa->get_brnlen()))/
			one_plus_two_beta);

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_2_BIAS),
			curr_exchange->get_dupl_parallel_rate()* 
			(exp(-1.0*one_plus_gamma*taxa->get_brnlen())-exp(-1.0*two_plus_two_beta_gamma*taxa->get_brnlen()))/
			one_plus_two_beta);

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_PRESENT_FIXED),
				 curr_exchange->get_dupl_fix_rate()*(1.0-exp(-1.0*two_plus_two_beta_gamma*taxa->get_brnlen()))/
				 two_plus_two_beta_gamma);

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_FIXED_SUBF),
				 2.0*curr_exchange->get_dupl_parallel_rate()*curr_exchange->get_dupl_fix_rate()*
				 (one_plus_two_beta+one_plus_gamma*exp(-1.0*two_plus_two_beta_gamma*taxa->get_brnlen())-
				 two_plus_two_beta_gamma*exp(-1.0*one_plus_gamma*taxa->get_brnlen()))/
				 (one_plus_two_beta*one_plus_gamma*two_plus_two_beta_gamma));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_PRESENT), 
			exp(-1.0*two_plus_two_beta_gamma*taxa->get_brnlen()));



}


int Dupl_Subf_All_States_Base_model::state_redundancy_val(DUPL_LOSS_STATES the_state)
{
		//Tells us if we need to sum over internal model states for a given input tip state

	if (observable_dupl_state(the_state)==(BOOL)TRUE) {
		switch (the_state) {
		case COPY1:
		case COPY2:
			return(2);
			break;
		case LOST:
			return(1);
			break;
		case BOTH_PRESENT:
			return(5);
			break;
		case MISSING:
			return(9);
			break;
		case GENERIC_SINGLE_COPY:
			return(4);
			break;
		case COPY1_OR_BOTH:
		case COPY2_OR_BOTH:
			return(7);
			break;
		}
	}
	else {
		cerr<<"ERROR: Tried to get Tip redundancy for non-observable internal model state\n";
		return(1);
	}

}



DUPL_LOSS_STATES Dupl_Subf_All_States_Base_model::get_redund_position_n(DUPL_LOSS_STATES the_state, int value)
{
		//Indicates order of the states to be summed through for ambiguous observed states

	switch(the_state) {
	
	//If it is not an ambiguous state, just return it
	case COPY1:
		if (value == 0)
			return(COPY1);
		else
			return(COPY1_BIAS);
		break;
	case COPY2:
		if (value == 0)
			return(COPY2);
		else
			return(COPY2_BIAS);
		break;
	case LOST:
		return(the_state);
		break;

	//Otherwise, give the possible states in an arbitrary order
	case BOTH_PRESENT:
		if (value == 0)
			return(BOTH_PRESENT);
		else if (value == 1)
			return(BOTH_1_BIAS);
		else if (value == 2)
			return(BOTH_2_BIAS);
		else if (value == 3)
			return(BOTH_PRESENT_FIXED);
		else
			return(BOTH_FIXED_SUBF);
		break;
	case MISSING:
		if (value == 0)
			return(BOTH_PRESENT);
		else if (value == 1)
			return(COPY1);
		else if (value == 2)
			return(COPY2);
		else if (value == 3)
			return(BOTH_1_BIAS);
		else if (value == 4)
			return(BOTH_2_BIAS);
		else if (value == 5)
			return(BOTH_PRESENT_FIXED);
		else if (value == 6)
			return(COPY1_BIAS);
		else if (value == 7)
			return(COPY2_BIAS);
		else
			return(BOTH_FIXED_SUBF);
		break;
	case GENERIC_SINGLE_COPY:
		if (value == 0)
			return(COPY1);
		else if (value == 1)
			return(COPY2);
		else if (value == 2)
			return(COPY1_BIAS);
		else
			return(COPY2_BIAS);
		break;
	case COPY1_OR_BOTH:
		if (value == 0)
			return(COPY1);
		else if (value == 1)
			return(BOTH_PRESENT);
		else if (value == 2)
			return(BOTH_PRESENT_FIXED);
		else if (value == 3)
			return(BOTH_1_BIAS);
		else if (value == 4)
			return(BOTH_2_BIAS);
		else if (value == 5)
			return(BOTH_FIXED_SUBF);
		else 
			return(COPY1_BIAS);
		break;
	case COPY2_OR_BOTH:
		if (value == 0)
			return(COPY2);
		else if (value == 1)
			return(BOTH_PRESENT);
		else if (value == 2)
			return(BOTH_PRESENT_FIXED);
		else if (value == 3)
			return(BOTH_1_BIAS);
		else if (value == 4)
			return(BOTH_2_BIAS);
		else if (value == 5)
			return(BOTH_FIXED_SUBF);
		else 
			return(COPY2_BIAS);
		
		break;


	}


}


BOOL Dupl_Subf_All_States_Base_model::allowed_state(DUPL_LOSS_STATES the_state)
{
	BOOL retval=(BOOL)FALSE;

	switch (the_state) {
	case BOTH_PRESENT:
	case BOTH_1_BIAS:
	case BOTH_2_BIAS:
	case BOTH_PRESENT_FIXED:
	case COPY1:
	case COPY2:
	case COPY1_BIAS:
	case COPY2_BIAS:
	case BOTH_FIXED_SUBF:
	case COPY1_OR_BOTH:
	case COPY2_OR_BOTH:
	case GENERIC_SINGLE_COPY:
	case MISSING:
	case LOST:
		retval=(BOOL)TRUE;
		break;

	}
	return(retval);
}


double Dupl_Subf_All_States_Base_model::find_ut (Branch *taxa)
{
	return(taxa->expect_subs_site()/(2.0+2.0*curr_exchange->get_dupl_parallel_rate()+curr_exchange->get_dupl_fix_rate()));
}


void Dupl_Slow_Loss_Con_Fix_Base_model::calc_transprobs(Branch *taxa, int rate_num)
{
		int i, j;
		double a, b, g, d, e;

		set_null_transprobs(taxa, rate_num);


		g=curr_exchange->get_dupl_fix_rate();
		b=curr_exchange->get_dupl_parallel_rate();
		a=taxa->get_brnlen();
		d=curr_exchange->get_loss_rate_scale();
		e=curr_exchange->get_fix_loss_rate();
	
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(COPY1), 1.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(COPY2), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_PRESENT_FIXED), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY1), loss_state_to_dupl(BOTH_2_BIAS), 0.0);

		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(COPY2), 1.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(COPY1), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_PRESENT_FIXED), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(COPY2), loss_state_to_dupl(BOTH_2_BIAS), 0.0);

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_PRESENT_FIXED), 
			0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_2_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(COPY1), 
			1.0 - exp(-1.0*a*d));
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(COPY2), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_1_BIAS), loss_state_to_dupl(BOTH_1_BIAS), 
			exp(-1.0*a*d));



		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_PRESENT_FIXED),
			0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(COPY1), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(COPY2), 
			1.0 - exp(-1.0*a*d));
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_2_BIAS), loss_state_to_dupl(BOTH_2_BIAS), 
			exp(-1.0*a*d));

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(COPY1), 
			0.5- 0.5*(exp(-2.0*a*e)));
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(COPY2), 
			0.5- 0.5*(exp(-2.0*a*e)));
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(BOTH_PRESENT), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(BOTH_1_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(BOTH_2_BIAS), 0.0);
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT_FIXED), loss_state_to_dupl(BOTH_PRESENT_FIXED), 
			exp(-2.0*a*e));


		
		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(COPY1),			
			((-exp(a*d)*g*(2 + 2*b + g - d)) + 2*exp(a*(-2 - 2*b - g + d + 2*e))*(b*(-2 + d) + (2 + g - d)*(-1 + e)) - 
			2*exp(2*a*e)*b*(2 + 2*b + g - 2*e) - exp(a*(d + 2*e))*(-2 - 2*b - g + d)*(2 + 2*b + g - 2*e))/
			(2.*exp(a*(d + 2*e))*(2 + 2*b + g - d)*(2 + 2*b + g - 2*e)) );

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(COPY2),
			((-exp(a*d)*g*(2 + 2*b + g - d)) + 2*exp(a*(-2 - 2*b - g + d + 2*e))*(b*(-2 + d) + (2 + g - d)*(-1 + e)) - 
			2*exp(2*a*e)*b*(2 + 2*b + g - 2*e) - exp(a*(d + 2*e))*(-2 - 2*b - g + d)*(2 + 2*b + g - 2*e))/
			(2.*exp(a*(d + 2*e))*(2 + 2*b + g - d)*(2 + 2*b + g - 2*e)) );

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_1_BIAS),
			(b - b/exp(a*(2 + 2*b + g - d)))/(exp(a*d)*(2 + 2*b + g - d)) );

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_2_BIAS),
			(b - b/exp(a*(2 + 2*b + g - d)))/(exp(a*d)*(2 + 2*b + g - d)) );

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_PRESENT_FIXED),
			-(((-1 + exp(-(a*(2 + 2*b + g - 2*e))))*g)/(exp(2*a*e)*(2 + 2*b + g - 2*e))) );

		taxa->set_trpb(rate_num, loss_state_to_dupl(BOTH_PRESENT), loss_state_to_dupl(BOTH_PRESENT), 
			exp(-(a*(2 + 2*b + g))) );



}




Dupl_Fix_Parallel_model::Dupl_Fix_Parallel_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree)
{
	if(cexchange->get_dupl_parallel_rate() ==0)
		cexchange->set_dupl_parallel_rate(0.1); 
	if(cexchange->get_dupl_fix_rate() == 0)
		cexchange->set_dupl_fix_rate(0.1);

	assemble (cexchange, cdata, ctree);
}

void Dupl_Fix_Parallel_model::describe_results()
{
	cout<<"Instantaneous relative rate of duplicate transition to parallel loss state: "
		<<curr_exchange->get_dupl_parallel_rate()<<"\n";
	cout<<"Instantaneous relative rate of duplicate fixation: "
		<<curr_exchange->get_dupl_fix_rate()<<"\n";
	
}



void Dupl_Fix_Parallel_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
	int i, brn_cnt=0;

	brn_start=0;
	for(i=0; i<curr_exchange->get_num_branches(); i++) {
		if ((curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE) || ((*curr_tree)[i]->get_brnlen() != 0)) {
				par[brn_cnt]=(*curr_tree)[i]->get_brnlen();
				types[brn_cnt]=BRANCH;
				if (curr_exchange->zero_len_brns_fixed() ==(BOOL)TRUE)
					brn_index[brn_cnt]=i;
				brn_cnt++;
		}
	}


	par[curr_exchange->get_num_params()-2]=curr_exchange->get_dupl_parallel_rate();
	types[curr_exchange->get_num_params()-2]=DUPL_PARALLEL_RATE;
	par[curr_exchange->get_num_params()-1]=curr_exchange->get_dupl_fix_rate();
	types[curr_exchange->get_num_params()-1]=DUPL_FIX_RATE;
}

    

void Dupl_Fix_Parallel_model::num_params_model() 
{
	int i, cnt_zero=0;

	if (curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE)
		curr_exchange->set_num_params(curr_exchange->get_num_branches()+2);
	else {
		for(i=0; i<curr_exchange->get_num_branches(); i++)
		{
			if ((*curr_tree)[i]->expect_subs_site() ==0)
				cnt_zero++;
		}
		brn_index=new int[curr_exchange->get_num_branches()-cnt_zero];
		curr_exchange->set_num_params(curr_exchange->get_num_branches()-cnt_zero+2);
	}
}



Dupl_Fix_Parallel_SubF_model::Dupl_Fix_Parallel_SubF_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree)
{
	if(cexchange->get_dupl_parallel_rate() ==0)
		cexchange->set_dupl_parallel_rate(0.1); 
	if(cexchange->get_dupl_fix_rate() == 0)
		cexchange->set_dupl_fix_rate(0.1);

	assemble (cexchange, cdata, ctree);
}
	

void Dupl_Fix_Parallel_SubF_model::describe_results()
{
	cout<<"Instantaneous relative rate of duplicate transition to parallel loss state: "
		<<curr_exchange->get_dupl_parallel_rate()<<"\n";
	cout<<"Instantaneous relative rate of duplicate fixation: "
		<<curr_exchange->get_dupl_fix_rate()<<"\n";
	
}
	

void Dupl_Fix_Parallel_SubF_model::num_params_model()
{
	int i, cnt_zero=0;

	if (curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE)
		curr_exchange->set_num_params(curr_exchange->get_num_branches()+2);
	else {
		for(i=0; i<curr_exchange->get_num_branches(); i++)
		{
			if ((*curr_tree)[i]->expect_subs_site() ==0)
				cnt_zero++;
		}
		brn_index=new int[curr_exchange->get_num_branches()-cnt_zero];
		curr_exchange->set_num_params(curr_exchange->get_num_branches()-cnt_zero+2);
	}
}


	
void Dupl_Fix_Parallel_SubF_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
	int i, brn_cnt=0;

	brn_start=0;
	for(i=0; i<curr_exchange->get_num_branches(); i++) {
		if ((curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE) || ((*curr_tree)[i]->get_brnlen() != 0)) {
				par[brn_cnt]=(*curr_tree)[i]->get_brnlen();
				types[brn_cnt]=BRANCH;
				if (curr_exchange->zero_len_brns_fixed() ==(BOOL)TRUE)
					brn_index[brn_cnt]=i;
				brn_cnt++;
		}
	}


	par[curr_exchange->get_num_params()-2]=curr_exchange->get_dupl_parallel_rate();
	types[curr_exchange->get_num_params()-2]=DUPL_PARALLEL_RATE;
	par[curr_exchange->get_num_params()-1]=curr_exchange->get_dupl_fix_rate();
	types[curr_exchange->get_num_params()-1]=DUPL_FIX_RATE;
}


Dupl_SubF_3_Rate_model::Dupl_SubF_3_Rate_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree)
{
	if(cexchange->get_dupl_parallel_rate() ==0)
		cexchange->set_dupl_parallel_rate(0.1); 
	if(cexchange->get_dupl_fix_rate() == 0)
		cexchange->set_dupl_fix_rate(0.1);
	if(cexchange->get_loss_rate_scale() ==0)
		cexchange->set_loss_rate_scale(1.0);
		if(cexchange->get_fix_rate_scale() ==0)
		cexchange->set_fix_rate_scale(1.0);

	assemble (cexchange, cdata, ctree);
}
	

void Dupl_SubF_3_Rate_model::describe_results()
{
	cout<<"Relative rate of fixation after converg.: "<<curr_exchange->get_fix_rate_scale()<<endl;
	cout<<"Relative rate of loss after converg.: "<<curr_exchange->get_loss_rate_scale()<<endl;
	cout<<"Instantaneous relative rate of duplicate transition to parallel loss state: "
		<<curr_exchange->get_dupl_parallel_rate()<<"\n";
	cout<<"Instantaneous relative rate of duplicate fixation: "
		<<curr_exchange->get_dupl_fix_rate()<<"\n";
	
}
	

void Dupl_SubF_3_Rate_model::num_params_model()
{
	int i, cnt_zero=0;

	if (curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE)
		curr_exchange->set_num_params(curr_exchange->get_num_branches()+4);
	else {
		for(i=0; i<curr_exchange->get_num_branches(); i++)
		{
			if ((*curr_tree)[i]->expect_subs_site() ==0)
				cnt_zero++;
		}
		brn_index=new int[curr_exchange->get_num_branches()-cnt_zero];
		curr_exchange->set_num_params(curr_exchange->get_num_branches()-cnt_zero+4);
	}
}


	
void Dupl_SubF_3_Rate_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
	int i, brn_cnt=0;

	brn_start=0;
	for(i=0; i<curr_exchange->get_num_branches(); i++) {
		if ((curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE) || ((*curr_tree)[i]->get_brnlen() != 0)) {
				par[brn_cnt]=(*curr_tree)[i]->get_brnlen();
				types[brn_cnt]=BRANCH;
				if (curr_exchange->zero_len_brns_fixed() ==(BOOL)TRUE)
					brn_index[brn_cnt]=i;
				brn_cnt++;
		}
	}


	par[curr_exchange->get_num_params()-4]=curr_exchange->get_dupl_parallel_rate();
	types[curr_exchange->get_num_params()-4]=DUPL_PARALLEL_RATE;
	par[curr_exchange->get_num_params()-3]=curr_exchange->get_dupl_fix_rate();
	types[curr_exchange->get_num_params()-3]=DUPL_FIX_RATE;

	par[curr_exchange->get_num_params()-2]=curr_exchange->get_fix_rate_scale();
	types[curr_exchange->get_num_params()-2]=CON_REL_FIX_RATE;
	par[curr_exchange->get_num_params()-1]=curr_exchange->get_loss_rate_scale();
	types[curr_exchange->get_num_params()-1]=CON_REL_LOSS_RATE;
}






Dupl_SubF_Only_model::Dupl_SubF_Only_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree)
{
	if(cexchange->get_dupl_parallel_rate() ==0)
		cexchange->set_dupl_parallel_rate(0.1); 
	if(cexchange->get_dupl_fix_rate() == 0)
		cexchange->set_dupl_fix_rate(0.1);
	if(cexchange->get_loss_rate_scale() ==0)
		cexchange->set_loss_rate_scale(1.0);
	
	assemble (cexchange, cdata, ctree);
}
	

void Dupl_SubF_Only_model::describe_results()
{
	cout<<"Relative rate of loss after converg.: "<<curr_exchange->get_loss_rate_scale()<<endl;
	cout<<"Instantaneous relative rate of duplicate transition to parallel loss state: "
		<<curr_exchange->get_dupl_parallel_rate()<<"\n";
	cout<<"Instantaneous relative rate of duplicate fixation: "
		<<curr_exchange->get_dupl_fix_rate()<<"\n";
	
}
	

void Dupl_SubF_Only_model::num_params_model()
{
	int i, cnt_zero=0;

	if (curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE)
		curr_exchange->set_num_params(curr_exchange->get_num_branches()+3);
	else {
		for(i=0; i<curr_exchange->get_num_branches(); i++)
		{
			if ((*curr_tree)[i]->expect_subs_site() ==0)
				cnt_zero++;
		}
		brn_index=new int[curr_exchange->get_num_branches()-cnt_zero];
		curr_exchange->set_num_params(curr_exchange->get_num_branches()-cnt_zero+3);
	}
}


	
void Dupl_SubF_Only_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
	int i, brn_cnt=0;

	brn_start=0;
	for(i=0; i<curr_exchange->get_num_branches(); i++) {
		if ((curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE) || ((*curr_tree)[i]->get_brnlen() != 0)) {
				par[brn_cnt]=(*curr_tree)[i]->get_brnlen();
				types[brn_cnt]=BRANCH;
				if (curr_exchange->zero_len_brns_fixed() ==(BOOL)TRUE)
					brn_index[brn_cnt]=i;
				brn_cnt++;
		}
	}


	par[curr_exchange->get_num_params()-3]=curr_exchange->get_dupl_parallel_rate();
	types[curr_exchange->get_num_params()-3]=DUPL_PARALLEL_RATE;
	par[curr_exchange->get_num_params()-2]=curr_exchange->get_dupl_fix_rate();
	types[curr_exchange->get_num_params()-2]=DUPL_FIX_RATE;

	par[curr_exchange->get_num_params()-1]=curr_exchange->get_loss_rate_scale();
	types[curr_exchange->get_num_params()-1]=CON_REL_LOSS_RATE;
}




Dupl_NoSubF_2_Rate_model::Dupl_NoSubF_2_Rate_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree)
{
	if(cexchange->get_dupl_parallel_rate() ==0)
		cexchange->set_dupl_parallel_rate(0.1); 
	if(cexchange->get_dupl_fix_rate() == 0)
		cexchange->set_dupl_fix_rate(0.1);
	if(cexchange->get_loss_rate_scale() ==0)
		cexchange->set_loss_rate_scale(1.0);
	
	assemble (cexchange, cdata, ctree);
}
	

void Dupl_NoSubF_2_Rate_model::describe_results()
{
	cout<<"The current model is a 6-state model without subfunctionalizing fixations, but with convergence and differing loss rates after convergence\n";
	cout<<"Relative rate of loss after converg.: "<<curr_exchange->get_loss_rate_scale()<<endl;
	cout<<"Instantaneous relative rate of duplicate transition to parallel loss state: "
		<<curr_exchange->get_dupl_parallel_rate()<<"\n";
	cout<<"Instantaneous relative rate of duplicate fixation: "
		<<curr_exchange->get_dupl_fix_rate()<<"\n";
	
}
	

void Dupl_NoSubF_2_Rate_model::num_params_model()
{
	int i, cnt_zero=0;

	if (curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE)
		curr_exchange->set_num_params(curr_exchange->get_num_branches()+3);
	else {
		for(i=0; i<curr_exchange->get_num_branches(); i++)
		{
			if ((*curr_tree)[i]->expect_subs_site() ==0)
				cnt_zero++;
		}
		brn_index=new int[curr_exchange->get_num_branches()-cnt_zero];
		curr_exchange->set_num_params(curr_exchange->get_num_branches()-cnt_zero+3);
	}
}


	
void Dupl_NoSubF_2_Rate_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
	int i, brn_cnt=0;

	brn_start=0;
	for(i=0; i<curr_exchange->get_num_branches(); i++) {
		if ((curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE) || ((*curr_tree)[i]->get_brnlen() != 0)) {
				par[brn_cnt]=(*curr_tree)[i]->get_brnlen();
				types[brn_cnt]=BRANCH;
				if (curr_exchange->zero_len_brns_fixed() ==(BOOL)TRUE)
					brn_index[brn_cnt]=i;
				brn_cnt++;
		}
	}


	par[curr_exchange->get_num_params()-3]=curr_exchange->get_dupl_parallel_rate();
	types[curr_exchange->get_num_params()-3]=DUPL_PARALLEL_RATE;
	par[curr_exchange->get_num_params()-2]=curr_exchange->get_dupl_fix_rate();
	types[curr_exchange->get_num_params()-2]=DUPL_FIX_RATE;

	par[curr_exchange->get_num_params()-1]=curr_exchange->get_loss_rate_scale();
	types[curr_exchange->get_num_params()-1]=CON_REL_LOSS_RATE;
}






Dupl_SubF_All_States_model::Dupl_SubF_All_States_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree)
{
	if(cexchange->get_dupl_parallel_rate() ==0)
		cexchange->set_dupl_parallel_rate(0.1); 
	if(cexchange->get_dupl_fix_rate() == 0)
		cexchange->set_dupl_fix_rate(0.1);

	assemble (cexchange, cdata, ctree);
}
	

void Dupl_SubF_All_States_model::describe_results()
{
	cout<<"Instantaneous relative rate of duplicate transition to parallel loss state: "
		<<curr_exchange->get_dupl_parallel_rate()<<"\n";
	cout<<"Instantaneous relative rate of duplicate fixation: "
		<<curr_exchange->get_dupl_fix_rate()<<"\n";
	
}

void Dupl_SubF_All_States_model::get_site_state_probs(double **&prob_array, int taxa_id)
{
	int i, j;
	double sum;


	prob_array=new double * [curr_exchange->get_num_localities()];
	
	for(i=0; i<curr_exchange->get_num_localities(); i++)
		prob_array[i] = new double [curr_exchange->get_condlike_size()];

	for(i=0; i<curr_exchange->get_num_localities(); i++)
		for(j=0; j<curr_exchange->get_condlike_size(); j++)
			prob_array[i][j]=0.0;


	for(i=0; i<curr_exchange->get_num_localities(); i++) {
		if (dupl_to_loss_state((*curr_data)[taxa_id][i]) == BOTH_PRESENT) {
			sum=0.0;

			(*curr_data)[taxa_id].Assign_site(i, loss_state_to_dupl(BOTH_PRESENT));

			partial_prob_w_rate_nonhidden(i, curr_tree->get_leftmost_tip(), 0, 
				  curr_tree->find_root(), taxa_id); 

			for (j=0; j<curr_exchange->get_condlike_size();j++)
				prob_array[i][loss_state_to_dupl(BOTH_PRESENT)]+=root_freq(loss_state_to_dupl(BOTH_PRESENT))*
					curr_tree->find_root()->get_trpb(0, loss_state_to_dupl(BOTH_PRESENT), j)*
					curr_tree->find_root()->get_cond_prob(j);
			sum+=prob_array[i][loss_state_to_dupl(BOTH_PRESENT)];

			(*curr_data)[taxa_id].Assign_site(i, loss_state_to_dupl(BOTH_1_BIAS));

			partial_prob_w_rate_nonhidden(i, curr_tree->get_leftmost_tip(), 0, 
				  curr_tree->find_root(), taxa_id); 

			for (j=0; j<curr_exchange->get_condlike_size();j++)
				prob_array[i][loss_state_to_dupl(BOTH_1_BIAS)]+=root_freq(loss_state_to_dupl(BOTH_PRESENT))*
					curr_tree->find_root()->get_trpb(0, loss_state_to_dupl(BOTH_PRESENT), j)*
					curr_tree->find_root()->get_cond_prob(j);
			sum+=prob_array[i][loss_state_to_dupl(BOTH_1_BIAS)];

			(*curr_data)[taxa_id].Assign_site(i, loss_state_to_dupl(BOTH_2_BIAS));

			partial_prob_w_rate_nonhidden(i, curr_tree->get_leftmost_tip(), 0, 
				  curr_tree->find_root(), taxa_id); 

			for (j=0; j<curr_exchange->get_condlike_size();j++)
				prob_array[i][loss_state_to_dupl(BOTH_2_BIAS)]+=root_freq(loss_state_to_dupl(BOTH_PRESENT))*
					curr_tree->find_root()->get_trpb(0, loss_state_to_dupl(BOTH_PRESENT), j)*
					curr_tree->find_root()->get_cond_prob(j);
			sum+=prob_array[i][loss_state_to_dupl(BOTH_2_BIAS)];

			(*curr_data)[taxa_id].Assign_site(i, loss_state_to_dupl(BOTH_PRESENT_FIXED));

			partial_prob_w_rate_nonhidden(i, curr_tree->get_leftmost_tip(), 0, 
				  curr_tree->find_root(), taxa_id); 

			for (j=0; j<curr_exchange->get_condlike_size();j++)
				prob_array[i][loss_state_to_dupl(BOTH_PRESENT_FIXED)]+=root_freq(loss_state_to_dupl(BOTH_PRESENT))*
					curr_tree->find_root()->get_trpb(0, loss_state_to_dupl(BOTH_PRESENT), j)*
					curr_tree->find_root()->get_cond_prob(j);
			sum+=prob_array[i][loss_state_to_dupl(BOTH_PRESENT_FIXED)];

			(*curr_data)[taxa_id].Assign_site(i, loss_state_to_dupl(BOTH_FIXED_SUBF));

			partial_prob_w_rate_nonhidden(i, curr_tree->get_leftmost_tip(), 0, 
				  curr_tree->find_root(), taxa_id); 

			for (j=0; j<curr_exchange->get_condlike_size();j++)
				prob_array[i][loss_state_to_dupl(BOTH_FIXED_SUBF)]+=root_freq(loss_state_to_dupl(BOTH_PRESENT))*
					curr_tree->find_root()->get_trpb(0, loss_state_to_dupl(BOTH_PRESENT), j)*
					curr_tree->find_root()->get_cond_prob(j);
			sum+=prob_array[i][loss_state_to_dupl(BOTH_FIXED_SUBF)];

			(*curr_data)[taxa_id].Assign_site(i, loss_state_to_dupl(BOTH_PRESENT));
		}
		else if (dupl_to_loss_state((*curr_data)[taxa_id][i]) == COPY1) {
			sum=0;
			
		}
		else if (dupl_to_loss_state((*curr_data)[taxa_id][i]) == COPY2) {
			sum=0;

			(*curr_data)[taxa_id].Assign_site(i, loss_state_to_dupl(COPY2));

			partial_prob_w_rate_nonhidden(i, curr_tree->get_leftmost_tip(), 0, 
				  curr_tree->find_root(), taxa_id); 

			for (j=0; j<curr_exchange->get_condlike_size();j++)
				prob_array[i][loss_state_to_dupl(COPY2)]+=root_freq(loss_state_to_dupl(BOTH_PRESENT))*
					curr_tree->find_root()->get_trpb(0, loss_state_to_dupl(BOTH_PRESENT), j)*
					curr_tree->find_root()->get_cond_prob(j);
			sum+=prob_array[i][loss_state_to_dupl(COPY2)];

			(*curr_data)[taxa_id].Assign_site(i, loss_state_to_dupl(COPY2_BIAS));

			partial_prob_w_rate_nonhidden(i, curr_tree->get_leftmost_tip(), 0, 
				  curr_tree->find_root(), taxa_id); 

			for (j=0; j<curr_exchange->get_condlike_size();j++)
				prob_array[i][loss_state_to_dupl(COPY2_BIAS)]+=root_freq(loss_state_to_dupl(BOTH_PRESENT))*
					curr_tree->find_root()->get_trpb(0, loss_state_to_dupl(BOTH_PRESENT), j)*
					curr_tree->find_root()->get_cond_prob(j);
			sum+=prob_array[i][loss_state_to_dupl(COPY2_BIAS)];
			(*curr_data)[taxa_id].Assign_site(i, loss_state_to_dupl(COPY2));
		}

	}



}
	

void Dupl_SubF_All_States_model::num_params_model()
{
	int i, cnt_zero=0;

	if (curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE)
		curr_exchange->set_num_params(curr_exchange->get_num_branches()+2);
	else {
		for(i=0; i<curr_exchange->get_num_branches(); i++)
		{
			if ((*curr_tree)[i]->expect_subs_site() ==0)
				cnt_zero++;
		}
		brn_index=new int[curr_exchange->get_num_branches()-cnt_zero];
		curr_exchange->set_num_params(curr_exchange->get_num_branches()-cnt_zero+2);
	}
}


	
void Dupl_SubF_All_States_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
	int i, brn_cnt=0;

	brn_start=0;
	for(i=0; i<curr_exchange->get_num_branches(); i++) {
		if ((curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE) || ((*curr_tree)[i]->get_brnlen() != 0)) {
				par[brn_cnt]=(*curr_tree)[i]->get_brnlen();
				types[brn_cnt]=BRANCH;
				if (curr_exchange->zero_len_brns_fixed() ==(BOOL)TRUE)
					brn_index[brn_cnt]=i;
				brn_cnt++;
		}
	}


	par[curr_exchange->get_num_params()-2]=curr_exchange->get_dupl_parallel_rate();
	types[curr_exchange->get_num_params()-2]=DUPL_PARALLEL_RATE;
	par[curr_exchange->get_num_params()-1]=curr_exchange->get_dupl_fix_rate();
	types[curr_exchange->get_num_params()-1]=DUPL_FIX_RATE;
}






Dupl_Fix_Parallel_NoState_model::Dupl_Fix_Parallel_NoState_model (Exchange *cexchange, Tree *ctree, Clade *cgenomes, WGD_Data *chomologs) 
{
	allocate_state_model(cexchange, ctree, cgenomes, chomologs); 
	
	if(cexchange->get_dupl_parallel_rate() ==0)
		cexchange->set_dupl_parallel_rate(0.1); 
	if(cexchange->get_dupl_fix_rate() == 0)
		cexchange->set_dupl_fix_rate(0.1);

	left_cond_probs=right_cond_probs=0;

	assemble (cexchange, curr_data, ctree);
}


void Dupl_Fix_Parallel_NoState_model::describe_results()
{
	cout<<"Instantaneous relative rate of duplicate transition to parallel loss state: "
		<<curr_exchange->get_dupl_parallel_rate()<<"\n";
		cout<<"Instantaneous relative rate of duplicate fixation: "
		<<curr_exchange->get_dupl_fix_rate()<<"\n";
		cout<<"Probability of strand definition switching between two single copy genes: "
		<<curr_exchange->get_strand_switch_prob()<<endl;
}



void Dupl_Fix_Parallel_NoState_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
	int i, brn_cnt=0;

	brn_start=0;
	for(i=0; i<curr_exchange->get_num_branches(); i++) {
		if ((curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE) || ((*curr_tree)[i]->get_brnlen() != 0)) {
				par[brn_cnt]=(*curr_tree)[i]->get_brnlen();
				types[brn_cnt]=BRANCH;
				if (curr_exchange->zero_len_brns_fixed() ==(BOOL)TRUE)
					brn_index[brn_cnt]=i;
				brn_cnt++;
		}
	}


	par[curr_exchange->get_num_params()-3]=log(curr_exchange->get_strand_switch_prob());
	types[curr_exchange->get_num_params()-3]=STRAND_SWITCH_PROB;
	par[curr_exchange->get_num_params()-2]=curr_exchange->get_dupl_parallel_rate();
	types[curr_exchange->get_num_params()-2]=DUPL_PARALLEL_RATE;
	par[curr_exchange->get_num_params()-1]=curr_exchange->get_dupl_fix_rate();
	types[curr_exchange->get_num_params()-1]=DUPL_FIX_RATE;
}

    

void Dupl_Fix_Parallel_NoState_model::num_params_model() 
{
	int i, cnt_zero=0;

	if (curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE)
		curr_exchange->set_num_params(curr_exchange->get_num_branches()+3);
	else {
		for(i=0; i<curr_exchange->get_num_branches(); i++)
		{
			if ((*curr_tree)[i]->expect_subs_site() ==0)
				cnt_zero++;
		}
		brn_index=new int[curr_exchange->get_num_branches()-cnt_zero];
		curr_exchange->set_num_params(curr_exchange->get_num_branches()-cnt_zero+3);
	}
}



Dupl_Fix_Parallel_SubF_NoState_model::Dupl_Fix_Parallel_SubF_NoState_model (Exchange *cexchange, Tree *ctree, Clade *cgenomes, WGD_Data *chomologs) 
{
	allocate_state_model(cexchange, ctree, cgenomes, chomologs); 
	
	if(cexchange->get_dupl_parallel_rate() ==0)
		cexchange->set_dupl_parallel_rate(0.1); 
	if(cexchange->get_dupl_fix_rate() == 0)
		cexchange->set_dupl_fix_rate(0.1);

	left_cond_probs=right_cond_probs=0;

	assemble (cexchange, curr_data, ctree);
}

void Dupl_Fix_Parallel_SubF_NoState_model::describe_results()
{
	cout<<"Instantaneous relative rate of duplicate transition to parallel loss state: "
		<<curr_exchange->get_dupl_parallel_rate()<<"\n";
		cout<<"Instantaneous relative rate of duplicate fixation: "
		<<curr_exchange->get_dupl_fix_rate()<<"\n";
		cout<<"Probability of strand definition switching between two single copy genes: "
		<<curr_exchange->get_strand_switch_prob()<<endl;
}



void Dupl_Fix_Parallel_SubF_NoState_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
	int i, brn_cnt=0;

	brn_start=0;
	for(i=0; i<curr_exchange->get_num_branches(); i++) {
		if ((curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE) || ((*curr_tree)[i]->get_brnlen() != 0)) {
				par[brn_cnt]=(*curr_tree)[i]->get_brnlen();
				types[brn_cnt]=BRANCH;
				if (curr_exchange->zero_len_brns_fixed() ==(BOOL)TRUE)
					brn_index[brn_cnt]=i;
				brn_cnt++;
		}
	}


	par[curr_exchange->get_num_params()-3]=log(curr_exchange->get_strand_switch_prob());
	types[curr_exchange->get_num_params()-3]=STRAND_SWITCH_PROB;
	par[curr_exchange->get_num_params()-2]=curr_exchange->get_dupl_parallel_rate();
	types[curr_exchange->get_num_params()-2]=DUPL_PARALLEL_RATE;
	par[curr_exchange->get_num_params()-1]=curr_exchange->get_dupl_fix_rate();
	types[curr_exchange->get_num_params()-1]=DUPL_FIX_RATE;
}

    

void Dupl_Fix_Parallel_SubF_NoState_model::num_params_model() 
{
	int i, cnt_zero=0;

	if (curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE)
		curr_exchange->set_num_params(curr_exchange->get_num_branches()+3);
	else {
		for(i=0; i<curr_exchange->get_num_branches(); i++)
		{
			if ((*curr_tree)[i]->expect_subs_site() ==0)
				cnt_zero++;
		}
		brn_index=new int[curr_exchange->get_num_branches()-cnt_zero];
		curr_exchange->set_num_params(curr_exchange->get_num_branches()-cnt_zero+3);
	}
}



Dupl_SubF_All_States_NoState_model::Dupl_SubF_All_States_NoState_model (Exchange *cexchange, Tree *ctree, Clade *cgenomes, WGD_Data *chomologs) 
{
	allocate_state_model(cexchange, ctree, cgenomes, chomologs); 
	
	if(cexchange->get_dupl_parallel_rate() ==0)
		cexchange->set_dupl_parallel_rate(0.1); 
	if(cexchange->get_dupl_fix_rate() == 0)
		cexchange->set_dupl_fix_rate(0.1);

	left_cond_probs=right_cond_probs=0;

	assemble (cexchange, curr_data, ctree);
}


void Dupl_SubF_All_States_NoState_model::describe_results()
{
	cout<<"Instantaneous relative rate of duplicate transition to parallel loss state: "
		<<curr_exchange->get_dupl_parallel_rate()<<"\n";
		cout<<"Instantaneous relative rate of duplicate fixation: "
		<<curr_exchange->get_dupl_fix_rate()<<"\n";
		cout<<"Probability of strand definition switching between two single copy genes: "
		<<curr_exchange->get_strand_switch_prob()<<endl;
}



void Dupl_SubF_All_States_NoState_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
	int i, brn_cnt=0;

	brn_start=0;
	for(i=0; i<curr_exchange->get_num_branches(); i++) {
		if ((curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE) || ((*curr_tree)[i]->get_brnlen() != 0)) {
				par[brn_cnt]=(*curr_tree)[i]->get_brnlen();
				types[brn_cnt]=BRANCH;
				if (curr_exchange->zero_len_brns_fixed() ==(BOOL)TRUE)
					brn_index[brn_cnt]=i;
				brn_cnt++;
		}
	}


	par[curr_exchange->get_num_params()-3]=log(curr_exchange->get_strand_switch_prob());
	types[curr_exchange->get_num_params()-3]=STRAND_SWITCH_PROB;
	par[curr_exchange->get_num_params()-2]=curr_exchange->get_dupl_parallel_rate();
	types[curr_exchange->get_num_params()-2]=DUPL_PARALLEL_RATE;
	par[curr_exchange->get_num_params()-1]=curr_exchange->get_dupl_fix_rate();
	types[curr_exchange->get_num_params()-1]=DUPL_FIX_RATE;
}

    

void Dupl_SubF_All_States_NoState_model::num_params_model() 
{
	int i, cnt_zero=0;

	if (curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE)
		curr_exchange->set_num_params(curr_exchange->get_num_branches()+3);
	else {
		for(i=0; i<curr_exchange->get_num_branches(); i++)
		{
			if ((*curr_tree)[i]->expect_subs_site() ==0)
				cnt_zero++;
		}
		brn_index=new int[curr_exchange->get_num_branches()-cnt_zero];
		curr_exchange->set_num_params(curr_exchange->get_num_branches()-cnt_zero+3);
	}
}




Dupl_SubF_Only_NoState_model::Dupl_SubF_Only_NoState_model (Exchange *cexchange, Tree *ctree, Clade *cgenomes, WGD_Data *chomologs) 
{
	allocate_state_model(cexchange, ctree, cgenomes, chomologs); 
	
	if(cexchange->get_dupl_parallel_rate() ==0)
		cexchange->set_dupl_parallel_rate(0.1); 
	if(cexchange->get_dupl_fix_rate() == 0)
		cexchange->set_dupl_fix_rate(0.1);

	if (cexchange->get_loss_rate_scale() == 0.0)
		cexchange->set_loss_rate_scale(1.0);

	left_cond_probs=right_cond_probs=0;

	assemble (cexchange, curr_data, ctree);
}


void Dupl_SubF_Only_NoState_model::describe_results()
{
		cout<<"The current model is a 6-state model with only subfunctionalizing fixations and with convergence and differing loss rates after convergence\n";
	cout<<"Relative rate of loss after converg.: "<<curr_exchange->get_loss_rate_scale()<<endl;
	cout<<"Instantaneous relative rate of duplicate transition to parallel loss state: "
		<<curr_exchange->get_dupl_parallel_rate()<<"\n";
		cout<<"Instantaneous relative rate of duplicate fixation: "
		<<curr_exchange->get_dupl_fix_rate()<<"\n";
		cout<<"Probability of strand definition switching between two single copy genes: "
		<<curr_exchange->get_strand_switch_prob()<<endl;
}



void Dupl_SubF_Only_NoState_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
	int i, brn_cnt=0;

	brn_start=0;
	for(i=0; i<curr_exchange->get_num_branches(); i++) {
		if ((curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE) || ((*curr_tree)[i]->get_brnlen() != 0)) {
				par[brn_cnt]=(*curr_tree)[i]->get_brnlen();
				types[brn_cnt]=BRANCH;
				if (curr_exchange->zero_len_brns_fixed() ==(BOOL)TRUE)
					brn_index[brn_cnt]=i;
				brn_cnt++;
		}
	}


	par[curr_exchange->get_num_params()-4]=curr_exchange->get_dupl_parallel_rate();
	types[curr_exchange->get_num_params()-4]=DUPL_PARALLEL_RATE;
	par[curr_exchange->get_num_params()-3]=curr_exchange->get_dupl_fix_rate();
	types[curr_exchange->get_num_params()-3]=DUPL_FIX_RATE;

	par[curr_exchange->get_num_params()-2]=curr_exchange->get_loss_rate_scale();
	types[curr_exchange->get_num_params()-2]=CON_REL_LOSS_RATE;

	par[curr_exchange->get_num_params()-1]=log(curr_exchange->get_strand_switch_prob());
	types[curr_exchange->get_num_params()-1]=STRAND_SWITCH_PROB;

}

    

void Dupl_SubF_Only_NoState_model::num_params_model() 
{
	int i, cnt_zero=0;

	if (curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE)
		curr_exchange->set_num_params(curr_exchange->get_num_branches()+4);
	else {
		for(i=0; i<curr_exchange->get_num_branches(); i++)
		{
			if ((*curr_tree)[i]->expect_subs_site() ==0)
				cnt_zero++;
		}
		brn_index=new int[curr_exchange->get_num_branches()-cnt_zero];
		curr_exchange->set_num_params(curr_exchange->get_num_branches()-cnt_zero+4);
	}
}






Dupl_SubF_3_Rate_NoState_model::Dupl_SubF_3_Rate_NoState_model (Exchange *cexchange, Tree *ctree, Clade *cgenomes, WGD_Data *chomologs) 
{
	allocate_state_model(cexchange, ctree, cgenomes, chomologs); 
	
	if(cexchange->get_dupl_parallel_rate() ==0)
		cexchange->set_dupl_parallel_rate(0.1); 
	if(cexchange->get_dupl_fix_rate() == 0)
		cexchange->set_dupl_fix_rate(0.1);

	cexchange->set_loss_rate_scale(1.0);
	cexchange->set_fix_rate_scale(cexchange->get_dupl_fix_rate());

	left_cond_probs=right_cond_probs=0;

	assemble (cexchange, curr_data, ctree);
}


void Dupl_SubF_3_Rate_NoState_model::describe_results()
{
		cout<<"The current model is a 6-state model with subfunctionalizing and non-subf. fixations and with convergence and differing loss rates after convergence\n";
	cout<<"Relative rate of fixation after converg.: "<<curr_exchange->get_fix_rate_scale()<<endl;
	cout<<"Relative rate of loss after converg.: "<<curr_exchange->get_loss_rate_scale()<<endl;
	cout<<"Instantaneous relative rate of duplicate transition to parallel loss state: "
		<<curr_exchange->get_dupl_parallel_rate()<<"\n";
		cout<<"Instantaneous relative rate of duplicate fixation: "
		<<curr_exchange->get_dupl_fix_rate()<<"\n";
		cout<<"Probability of strand definition switching between two single copy genes: "
		<<curr_exchange->get_strand_switch_prob()<<endl;
}



void Dupl_SubF_3_Rate_NoState_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
	int i, brn_cnt=0;

	brn_start=0;
	for(i=0; i<curr_exchange->get_num_branches(); i++) {
		if ((curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE) || ((*curr_tree)[i]->get_brnlen() != 0)) {
				par[brn_cnt]=(*curr_tree)[i]->get_brnlen();
				types[brn_cnt]=BRANCH;
				if (curr_exchange->zero_len_brns_fixed() ==(BOOL)TRUE)
					brn_index[brn_cnt]=i;
				brn_cnt++;
		}
	}


	par[curr_exchange->get_num_params()-5]=curr_exchange->get_dupl_parallel_rate();
	types[curr_exchange->get_num_params()-5]=DUPL_PARALLEL_RATE;
	par[curr_exchange->get_num_params()-4]=curr_exchange->get_dupl_fix_rate();
	types[curr_exchange->get_num_params()-4]=DUPL_FIX_RATE;

	par[curr_exchange->get_num_params()-3]=curr_exchange->get_fix_rate_scale();
	types[curr_exchange->get_num_params()-3]=CON_REL_FIX_RATE;
	par[curr_exchange->get_num_params()-2]=curr_exchange->get_loss_rate_scale();
	types[curr_exchange->get_num_params()-2]=CON_REL_LOSS_RATE;

	par[curr_exchange->get_num_params()-1]=log(curr_exchange->get_strand_switch_prob());
	types[curr_exchange->get_num_params()-1]=STRAND_SWITCH_PROB;

}

    

void Dupl_SubF_3_Rate_NoState_model::num_params_model() 
{
	int i, cnt_zero=0;

	if (curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE)
		curr_exchange->set_num_params(curr_exchange->get_num_branches()+5);
	else {
		for(i=0; i<curr_exchange->get_num_branches(); i++)
		{
			if ((*curr_tree)[i]->expect_subs_site() ==0)
				cnt_zero++;
		}
		brn_index=new int[curr_exchange->get_num_branches()-cnt_zero];
		curr_exchange->set_num_params(curr_exchange->get_num_branches()-cnt_zero+5);
	}
}


Dupl_SubF_3_Rate_AllStates_NoState_model::Dupl_SubF_3_Rate_AllStates_NoState_model (Exchange *cexchange, Tree *ctree, Clade *cgenomes, WGD_Data *chomologs) 
{
	allocate_state_model(cexchange, ctree, cgenomes, chomologs); 
	
	if(cexchange->get_dupl_parallel_rate() ==0)
		cexchange->set_dupl_parallel_rate(0.1); 
	if(cexchange->get_dupl_fix_rate() == 0)
		cexchange->set_dupl_fix_rate(0.1);

	cexchange->set_loss_rate_scale(1.0);
	cexchange->set_fix_rate_scale(cexchange->get_dupl_fix_rate());

	left_cond_probs=right_cond_probs=0;

	assemble (cexchange, curr_data, ctree);
}


void Dupl_SubF_3_Rate_AllStates_NoState_model::describe_results()
{
		cout<<"The current model is a 6-state model without subfunctionalizing and non-subf fixations and with convergence and differing loss rates after convergence\n";
	cout<<"Relative rate of fixation after converg.: "<<curr_exchange->get_fix_rate_scale()<<endl;
	cout<<"Relative rate of loss after converg.: "<<curr_exchange->get_loss_rate_scale()<<endl;
	cout<<"Instantaneous relative rate of duplicate transition to parallel loss state: "
		<<curr_exchange->get_dupl_parallel_rate()<<"\n";
		cout<<"Instantaneous relative rate of duplicate fixation: "
		<<curr_exchange->get_dupl_fix_rate()<<"\n";
		cout<<"Probability of strand definition switching between two single copy genes: "
		<<curr_exchange->get_strand_switch_prob()<<endl;
}



void Dupl_SubF_3_Rate_AllStates_NoState_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
	int i, brn_cnt=0;

	brn_start=0;
	for(i=0; i<curr_exchange->get_num_branches(); i++) {
		if ((curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE) || ((*curr_tree)[i]->get_brnlen() != 0)) {
				par[brn_cnt]=(*curr_tree)[i]->get_brnlen();
				types[brn_cnt]=BRANCH;
				if (curr_exchange->zero_len_brns_fixed() ==(BOOL)TRUE)
					brn_index[brn_cnt]=i;
				brn_cnt++;
		}
	}


	par[curr_exchange->get_num_params()-5]=curr_exchange->get_dupl_parallel_rate();
	types[curr_exchange->get_num_params()-5]=DUPL_PARALLEL_RATE;
	par[curr_exchange->get_num_params()-4]=curr_exchange->get_dupl_fix_rate();
	types[curr_exchange->get_num_params()-4]=DUPL_FIX_RATE;

	par[curr_exchange->get_num_params()-3]=curr_exchange->get_fix_rate_scale();
	types[curr_exchange->get_num_params()-3]=CON_REL_FIX_RATE;
	par[curr_exchange->get_num_params()-2]=curr_exchange->get_loss_rate_scale();
	types[curr_exchange->get_num_params()-2]=CON_REL_LOSS_RATE;

	par[curr_exchange->get_num_params()-1]=log(curr_exchange->get_strand_switch_prob());
	types[curr_exchange->get_num_params()-1]=STRAND_SWITCH_PROB;

}

    

void Dupl_SubF_3_Rate_AllStates_NoState_model::num_params_model() 
{
	int i, cnt_zero=0;

	if (curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE)
		curr_exchange->set_num_params(curr_exchange->get_num_branches()+5);
	else {
		for(i=0; i<curr_exchange->get_num_branches(); i++)
		{
			if ((*curr_tree)[i]->expect_subs_site() ==0)
				cnt_zero++;
		}
		brn_index=new int[curr_exchange->get_num_branches()-cnt_zero];
		curr_exchange->set_num_params(curr_exchange->get_num_branches()-cnt_zero+5);
	}
}



Dupl_2_Rate_NoSubF_NoState_model::Dupl_2_Rate_NoSubF_NoState_model (Exchange *cexchange, Tree *ctree, Clade *cgenomes, WGD_Data *chomologs) 
{
	allocate_state_model(cexchange, ctree, cgenomes, chomologs); 
	
	if(cexchange->get_dupl_parallel_rate() ==0)
		cexchange->set_dupl_parallel_rate(0.1); 
	if(cexchange->get_dupl_fix_rate() == 0)
		cexchange->set_dupl_fix_rate(0.1);

	if (cexchange->get_loss_rate_scale() == 0)
		cexchange->set_loss_rate_scale(1.0);

	left_cond_probs=right_cond_probs=0;

	assemble (cexchange, curr_data, ctree);
}


void Dupl_2_Rate_NoSubF_NoState_model::describe_results()
{
	cout<<"The current model is a 6-state model without subfunctionalizing fixations, but with convergence and differing loss rates after convergence\n";
	cout<<"Relative rate of loss after converg.: "<<curr_exchange->get_loss_rate_scale()<<endl;
	cout<<"Instantaneous relative rate of duplicate transition to parallel loss state: "
		<<curr_exchange->get_dupl_parallel_rate()<<"\n";
		cout<<"Instantaneous relative rate of duplicate fixation: "
		<<curr_exchange->get_dupl_fix_rate()<<"\n";
		cout<<"Probability of strand definition switching between two single copy genes: "
		<<curr_exchange->get_strand_switch_prob()<<endl;
}



void Dupl_2_Rate_NoSubF_NoState_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
	int i, brn_cnt=0;

	par[0]=curr_exchange->get_loss_rate_scale();
	types[0]=CON_REL_LOSS_RATE;
	
	
	par[1]=curr_exchange->get_dupl_parallel_rate();
	types[1]=DUPL_PARALLEL_RATE;
	par[2]=curr_exchange->get_dupl_fix_rate();
	types[2]=DUPL_FIX_RATE;
	
	
	par[3]=log(curr_exchange->get_strand_switch_prob());
	types[3]=STRAND_SWITCH_PROB;
	
	brn_start=4;
	for(i=0; i<curr_exchange->get_num_branches(); i++) {
		if ((curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE) || ((*curr_tree)[i]->get_brnlen() != 0)) {
				par[brn_cnt+brn_start]=(*curr_tree)[i]->get_brnlen()*2.0;
				types[brn_cnt+brn_start]=BRANCH;
				if (curr_exchange->zero_len_brns_fixed() ==(BOOL)TRUE)
					brn_index[brn_cnt]=i;
				brn_cnt++;
		}
	}

	

}

    

void Dupl_2_Rate_NoSubF_NoState_model::num_params_model() 
{
	int i, cnt_zero=0;

	if (curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE)
		curr_exchange->set_num_params(curr_exchange->get_num_branches()+4);
	else {
		for(i=0; i<curr_exchange->get_num_branches(); i++)
		{
			if ((*curr_tree)[i]->expect_subs_site() ==0)
				cnt_zero++;
		}
		brn_index=new int[curr_exchange->get_num_branches()-cnt_zero];
		curr_exchange->set_num_params(curr_exchange->get_num_branches()-cnt_zero+4);
	}
}



Dupl_2_Rate_NoSubF_AllStates_NoState_model::Dupl_2_Rate_NoSubF_AllStates_NoState_model (Exchange *cexchange, Tree *ctree, Clade *cgenomes, WGD_Data *chomologs) 
{
	allocate_state_model(cexchange, ctree, cgenomes, chomologs); 
	
	if(cexchange->get_dupl_parallel_rate() ==0)
		cexchange->set_dupl_parallel_rate(0.1); 
	if(cexchange->get_dupl_fix_rate() == 0)
		cexchange->set_dupl_fix_rate(0.1);

	if (cexchange->get_loss_rate_scale() == 0)
		cexchange->set_loss_rate_scale(1.0);

	left_cond_probs=right_cond_probs=0;


	assemble (cexchange, curr_data, ctree);
}


void Dupl_2_Rate_NoSubF_AllStates_NoState_model::describe_results()
{
	cout<<"Relative rate of loss after converg.: "<<curr_exchange->get_loss_rate_scale()<<endl;
	cout<<"Instantaneous relative rate of duplicate transition to parallel loss state: "
		<<curr_exchange->get_dupl_parallel_rate()<<"\n";
		cout<<"Instantaneous relative rate of duplicate fixation: "
		<<curr_exchange->get_dupl_fix_rate()<<"\n";
		cout<<"Probability of strand definition switching between two single copy genes: "
		<<curr_exchange->get_strand_switch_prob()<<endl;
}



void Dupl_2_Rate_NoSubF_AllStates_NoState_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
	int i, brn_cnt=0;

	brn_start=0;
	for(i=0; i<curr_exchange->get_num_branches(); i++) {
		if ((curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE) || ((*curr_tree)[i]->get_brnlen() != 0)) {
				par[brn_cnt]=(*curr_tree)[i]->get_brnlen();
				types[brn_cnt]=BRANCH;
				if (curr_exchange->zero_len_brns_fixed() ==(BOOL)TRUE)
					brn_index[brn_cnt]=i;
				brn_cnt++;
		}
	}


	par[curr_exchange->get_num_params()-4]=curr_exchange->get_dupl_parallel_rate();
	types[curr_exchange->get_num_params()-4]=DUPL_PARALLEL_RATE;
	par[curr_exchange->get_num_params()-3]=curr_exchange->get_dupl_fix_rate();
	types[curr_exchange->get_num_params()-3]=DUPL_FIX_RATE;


	par[curr_exchange->get_num_params()-2]=curr_exchange->get_loss_rate_scale();
	types[curr_exchange->get_num_params()-2]=CON_REL_LOSS_RATE;

	par[curr_exchange->get_num_params()-1]=log(curr_exchange->get_strand_switch_prob());
	types[curr_exchange->get_num_params()-1]=STRAND_SWITCH_PROB;

}

    

void Dupl_2_Rate_NoSubF_AllStates_NoState_model::num_params_model() 
{
	int i, cnt_zero=0;

	if (curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE)
		curr_exchange->set_num_params(curr_exchange->get_num_branches()+4);
	else {
		for(i=0; i<curr_exchange->get_num_branches(); i++)
		{
			if ((*curr_tree)[i]->expect_subs_site() ==0)
				cnt_zero++;
		}
		brn_index=new int[curr_exchange->get_num_branches()-cnt_zero];
		curr_exchange->set_num_params(curr_exchange->get_num_branches()-cnt_zero+4);
	}
}





Dupl_Slow_Loss_Con_Fix_model::Dupl_Slow_Loss_Con_Fix_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree)
{
	if(cexchange->get_dupl_parallel_rate() ==0)
		cexchange->set_dupl_parallel_rate(0.1); 
	if(cexchange->get_dupl_fix_rate() == 0)
		cexchange->set_dupl_fix_rate(0.1);
	if(cexchange->get_loss_rate_scale() ==0)
		cexchange->set_loss_rate_scale(1.0);
	if(cexchange->get_fix_loss_rate() ==0)
		cexchange->set_fix_loss_rate(1.0);
	
	assemble (cexchange, cdata, ctree);
}
	

void Dupl_Slow_Loss_Con_Fix_model::describe_results()
{
	cout<<"The current model is a 6-state model with slow losses from the convergent and pseudo-fixed duplicates\n";
	cout<<"Relative rate of loss after converg.: "<<curr_exchange->get_loss_rate_scale()<<endl;
	cout<<"Instantaneous relative rate of duplicate transition to parallel loss state: "
		<<curr_exchange->get_dupl_parallel_rate()<<"\n";
	cout<<"Instantaneous relative rate of duplicate pseudo-fixation: "
		<<curr_exchange->get_dupl_fix_rate()<<"\n";
	cout<<"Instantaneous relative rate of duplicate loss from pseudo-fixation: "
		<<curr_exchange->get_fix_loss_rate()<<"\n";
	
}
	

void Dupl_Slow_Loss_Con_Fix_model::num_params_model()
{
	int i, cnt_zero=0;

	if (curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE)
		curr_exchange->set_num_params(curr_exchange->get_num_branches()+4);
	else {
		for(i=0; i<curr_exchange->get_num_branches(); i++)
		{
			if ((*curr_tree)[i]->expect_subs_site() ==0)
				cnt_zero++;
		}
		brn_index=new int[curr_exchange->get_num_branches()-cnt_zero];
		curr_exchange->set_num_params(curr_exchange->get_num_branches()-cnt_zero+4);
	}
}


	
void Dupl_Slow_Loss_Con_Fix_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
	int i, brn_cnt=0;

	brn_start=0;
	for(i=0; i<curr_exchange->get_num_branches(); i++) {
		if ((curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE) || ((*curr_tree)[i]->get_brnlen() != 0)) {
				par[brn_cnt]=(*curr_tree)[i]->get_brnlen();
				types[brn_cnt]=BRANCH;
				if (curr_exchange->zero_len_brns_fixed() ==(BOOL)TRUE)
					brn_index[brn_cnt]=i;
				brn_cnt++;
		}
	}

	par[curr_exchange->get_num_params()-4]=curr_exchange->get_fix_loss_rate();
	types[curr_exchange->get_num_params()-4]=CON_FIX_LOSS_RATE;

	par[curr_exchange->get_num_params()-3]=curr_exchange->get_dupl_parallel_rate();
	types[curr_exchange->get_num_params()-3]=DUPL_PARALLEL_RATE;
	par[curr_exchange->get_num_params()-2]=curr_exchange->get_dupl_fix_rate();
	types[curr_exchange->get_num_params()-2]=DUPL_FIX_RATE;

	par[curr_exchange->get_num_params()-1]=curr_exchange->get_loss_rate_scale();
	types[curr_exchange->get_num_params()-1]=CON_REL_LOSS_RATE;
}






Dupl_Slow_Loss_Con_Fix_NoState_model::Dupl_Slow_Loss_Con_Fix_NoState_model (Exchange *cexchange, Tree *ctree, Clade *cgenomes, WGD_Data *chomologs) 
{
	allocate_state_model(cexchange, ctree, cgenomes, chomologs); 
	
	if(cexchange->get_dupl_parallel_rate() ==0)
		cexchange->set_dupl_parallel_rate(0.1); 
	if(cexchange->get_dupl_fix_rate() == 0)
		cexchange->set_dupl_fix_rate(0.1);

	if (cexchange->get_loss_rate_scale() == 0)
		cexchange->set_loss_rate_scale(1.0);

	if(cexchange->get_fix_loss_rate() == 0)
		cexchange->set_fix_loss_rate(1.0);

	left_cond_probs=right_cond_probs=0;

	assemble (cexchange, curr_data, ctree);
}


void Dupl_Slow_Loss_Con_Fix_NoState_model::describe_results()
{
	cout<<"The current model is a 6-state model with slow losses from the convergent and pseudo-fixed duplicates and inferred tracking\n";
	cout<<"Relative rate of loss after converg.: "<<curr_exchange->get_loss_rate_scale()<<endl;
	cout<<"Instantaneous relative rate of duplicate transition to parallel loss state: "
		<<curr_exchange->get_dupl_parallel_rate()<<"\n";
	cout<<"Instantaneous relative rate of duplicate pseudo-fixation: "
		<<curr_exchange->get_dupl_fix_rate()<<"\n";
	cout<<"Instantaneous relative rate of duplicate loss from pseudo-fixation: "
		<<curr_exchange->get_fix_loss_rate()<<"\n";
	cout<<"Probability of strand definition switching between two single copy genes: "
		<<curr_exchange->get_strand_switch_prob()<<endl;
}



void Dupl_Slow_Loss_Con_Fix_NoState_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
	int i, brn_cnt=0;

	brn_start=0;
	for(i=0; i<curr_exchange->get_num_branches(); i++) {
		if ((curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE) || ((*curr_tree)[i]->get_brnlen() != 0)) {
				par[brn_cnt]=(*curr_tree)[i]->get_brnlen();
				types[brn_cnt]=BRANCH;
				if (curr_exchange->zero_len_brns_fixed() ==(BOOL)TRUE)
					brn_index[brn_cnt]=i;
				brn_cnt++;
		}
	}

	par[curr_exchange->get_num_params()-5]=curr_exchange->get_fix_loss_rate();
	types[curr_exchange->get_num_params()-5]=CON_FIX_LOSS_RATE;

	par[curr_exchange->get_num_params()-4]=curr_exchange->get_dupl_parallel_rate();
	types[curr_exchange->get_num_params()-4]=DUPL_PARALLEL_RATE;
	par[curr_exchange->get_num_params()-3]=curr_exchange->get_dupl_fix_rate();
	types[curr_exchange->get_num_params()-3]=DUPL_FIX_RATE;

		par[curr_exchange->get_num_params()-2]=curr_exchange->get_loss_rate_scale();
	types[curr_exchange->get_num_params()-2]=CON_REL_LOSS_RATE;

	par[curr_exchange->get_num_params()-1]=log(curr_exchange->get_strand_switch_prob());
	types[curr_exchange->get_num_params()-1]=STRAND_SWITCH_PROB;

}

    

void Dupl_Slow_Loss_Con_Fix_NoState_model::num_params_model() 
{
	int i, cnt_zero=0;

	if (curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE)
		curr_exchange->set_num_params(curr_exchange->get_num_branches()+5);
	else {
		for(i=0; i<curr_exchange->get_num_branches(); i++)
		{
			if ((*curr_tree)[i]->expect_subs_site() ==0)
				cnt_zero++;
		}
		brn_index=new int[curr_exchange->get_num_branches()-cnt_zero];
		curr_exchange->set_num_params(curr_exchange->get_num_branches()-cnt_zero+5);
	}
}



//NEW CODE!!!
SNP_model::SNP_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree)
{
	assemble (cexchange, cdata, ctree);
}


double SNP_model::root_freq(int site)
{
	//Note that we make all states equally probable at the root here--this
	//is not terribly realistic--should do a parameter
	switch (snp_to_snpstate(site)) {
	case DATA_MISSING:
		return(0.0);
		break;
	case TYPE_A:
	case TYPE_B:
	case HETEROZYGOUS:
	case SNP_ABSENT:
		return(0.25);
		break;
	default:
		return(0.0);
		break;
	}

}

void SNP_model::initialize_arrays()
{
}
 


void SNP_model::partial_prob_w_rate(int locale, Branch *lsib, int rate_num, Branch *stop_id)
 //The heart of the likelihood calculation: uses the tranisition probablity matrices stored in
  //the tree to calculate the likelihood of part of that tree.  Uses post-order tree transversal.
  //This is modified version for the Dupl_models that allows the tip states to be ambiguous--hence
  //probability is summed across those ambiguous states
{
  int i, j, k;
  long double temp_p1, temp_p2, remove_conprob, scale_factor;
  Branch *rsib, *parent;
  BOOL all_zero;
  
  rsib=lsib->get_sibling();
  parent=lsib->get_parent();
  

  
  if (lsib->is_tip()==(BOOL)FALSE)
    {
      if (rsib->is_tip()==(BOOL)FALSE) 
		{
			partial_prob_w_rate(locale, curr_tree->find_left_tip(rsib), rate_num, rsib);
			for(i=0; i<curr_exchange->get_condlike_size(); i++)
			{
				temp_p1=temp_p2=0;
				for(j=0; j<curr_exchange->get_condlike_size(); j++)
				{
					temp_p1+=lsib->get_cond_prob(j)*(long double)lsib->get_trpb(rate_num, i, j);
					temp_p2+=rsib->get_cond_prob(j)*(long double)rsib->get_trpb(rate_num, i, j);
				} 
				parent->set_cond_prob(i, (temp_p1/(long double)SQRT_LIKE_SCALE)*(temp_p2/(long double)SQRT_LIKE_SCALE));
			}	   
		}
    
      else
		{
			for(i=0; i<curr_exchange->get_condlike_size(); i++)
			{
				temp_p1=0;
				for(j=0; j<curr_exchange->get_condlike_size(); j++) 
					temp_p1+=lsib->get_cond_prob(j)*(long double)lsib->get_trpb(rate_num, i, j);
		  
				parent->set_cond_prob(i, temp_p1*(long double)(get_tip_prob(locale, rsib, rate_num, i))); 
			}
		}
    }

  else
    {
      if (rsib->is_tip()==(BOOL)FALSE) 
		{
			partial_prob_w_rate(locale, curr_tree->find_left_tip(rsib), rate_num, rsib);
	
			for(i=0; i<curr_exchange->get_condlike_size(); i++)
			{
				temp_p2=0;
				
				for(j=0; j<curr_exchange->get_condlike_size(); j++)
					temp_p2+=rsib->get_cond_prob(j)*(long double)rsib->get_trpb(rate_num, i, j); 
		  
				parent->set_cond_prob(i, (long double)(get_tip_prob(locale, lsib, rate_num, i))*temp_p2);  
			}   
		}

      else
		  for(i=0; i<curr_exchange->get_condlike_size(); i++) 
				parent->set_cond_prob(i, (long double)LIKE_SCALE*(long double)(get_tip_prob(locale, lsib, rate_num, i))*
				(long double)(get_tip_prob(locale, rsib, rate_num, i)));  
	  
    }

  //all_zero=(BOOL)TRUE;

  //for(i=0; i<curr_exchange->get_condlike_size(); i++)
	 // if (parent->get_cond_prob(i) != 0.0) all_zero=(BOOL)FALSE;
  

 // if(all_zero == (BOOL)TRUE) 
	//cerr<<"Failure point\n";


	if (curr_tree->get_constrain_brn() != 0) {
		if (curr_tree->get_constrain_brn() == parent) {
			for(i=0; i<curr_exchange->num_unallowed_states(); i++)
				parent->set_cond_prob(snpstate_to_snp(curr_exchange->get_unallowed_state(i)), 0.0);
		}
	}
	
  if (parent!=stop_id)
      partial_prob_w_rate(locale, parent, rate_num, stop_id);
 
}  //End Like_model::partial_prob_w_rate









double SNP_model::get_tip_prob(int locale, Branch *taxa, int rate_num, int cond_state)
{	
	double ret_val;

	if (snp_to_snpstate((*curr_data)[taxa->get_taxa_id()][locale]) != DATA_MISSING) 
		ret_val = taxa->get_trpb(rate_num, cond_state, (*curr_data)[taxa->get_taxa_id()][locale]);
	
	else 
		ret_val=
			taxa->get_trpb(rate_num, cond_state, snpstate_to_snp(TYPE_A))+taxa->get_trpb(rate_num, cond_state, snpstate_to_snp(TYPE_B))+
			taxa->get_trpb(rate_num, cond_state, snpstate_to_snp(HETEROZYGOUS))+taxa->get_trpb(rate_num, cond_state, snpstate_to_snp(SNP_ABSENT));

	return(ret_val);	
}



void SNP_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
	int i;


	par[0]=log(curr_exchange->get_snp_gains_to_losses());
	types[0]=SNP_GAIN_LOSS_RATIO;

	par[1]=log(curr_exchange->get_snp_loss_to_switch());
	types[1]=SNP_LOSS_SWITCH_RATIO;
}

void SNP_model::describe_results()
{
	cout<<"Ratio of gains to losses: "<<curr_exchange->get_snp_gains_to_losses()<<endl;
	cout<<"Ratio of losses to SNP state changes: "<<curr_exchange->get_snp_loss_to_switch()<<endl;
}


void SNP_model::num_params_model() 
{
	curr_exchange->set_num_params(2);	
}


void SNP_model::calc_transprobs(Branch *taxa, int rate_num)
{
	int i, j;
	double beta, alpha;

	//Force the root branch length to 0 in all cases
	if (taxa == curr_tree->find_root())  taxa->set_brnlen(0.0);

	if (taxa->get_brnlen() <= MIN_BRLEN) {
		taxa->set_brnlen(0.0);
		for (i=0;i<curr_exchange->get_condlike_size(); i++){
			for (j=0;j<curr_exchange->get_condlike_size(); j++){
				if (i==j) taxa->set_trpb(rate_num,i, j, 1.0);
				else taxa->set_trpb(rate_num, i,j , 0.0);

			}
		}
	}
	else {
		for (i=0;i<curr_exchange->get_condlike_size(); i++){
			if (snp_to_snpstate(i) != DATA_MISSING) {
				taxa->set_trpb(rate_num, i,snpstate_to_snp(DATA_MISSING), 0.0);
				taxa->set_trpb(rate_num, snpstate_to_snp(DATA_MISSING), i, 0.0);
			}
			else 
				taxa->set_trpb(rate_num, snpstate_to_snp(DATA_MISSING),snpstate_to_snp(DATA_MISSING), 1.0);
		}

		beta=curr_exchange->get_snp_gains_to_losses();
		alpha=curr_exchange->get_snp_loss_to_switch();

		taxa->set_trpb(rate_num, snpstate_to_snp(TYPE_A), snpstate_to_snp(TYPE_A), 
			beta/(3.0*beta+1.0)+ 1.0/6.0*exp(-1.0*(3.0*alpha+1.0)*taxa->get_brnlen())+
			0.5*exp(-1.0*(alpha+1.0)*taxa->get_brnlen())+ exp(-1.0*(3.0*beta+1.0)*taxa->get_brnlen())/(3.0*(3.0*beta+1)));
		taxa->set_trpb(rate_num, snpstate_to_snp(TYPE_B), snpstate_to_snp(TYPE_B), 
			beta/(3.0*beta+1.0)+ 1.0/6.0*exp(-1.0*(3.0*alpha+1.0)*taxa->get_brnlen())+
			0.5*exp(-1.0*(alpha+1.0)*taxa->get_brnlen())+ exp(-1.0*(3.0*beta+1.0)*taxa->get_brnlen())/(3.0*(3.0*beta+1)));

		taxa->set_trpb(rate_num, snpstate_to_snp(TYPE_A), snpstate_to_snp(TYPE_B), 
			beta/(3.0*beta+1.0)+1.0/6.0*exp(-1.0*(3.0*alpha+1.0)*taxa->get_brnlen())-
			0.5*exp(-1.0*(alpha+1.0)*taxa->get_brnlen())+exp(-1.0*(3.0*beta+1.0)*taxa->get_brnlen())/(3.0*(3.0*beta+1)));
		taxa->set_trpb(rate_num, snpstate_to_snp(TYPE_B), snpstate_to_snp(TYPE_A), 
			beta/(3.0*beta+1.0)+1.0/6.0*exp(-1.0*(3.0*alpha+1.0)*taxa->get_brnlen())-
			0.5*exp(-1.0*(alpha+1.0)*taxa->get_brnlen())+exp(-1.0*(3.0*beta+1.0)*taxa->get_brnlen())/(3.0*(3.0*beta+1)));


		taxa->set_trpb(rate_num, snpstate_to_snp(TYPE_A), snpstate_to_snp(HETEROZYGOUS), 
			beta/(3.0*beta+1.0)-1.0/3.0*exp(-1.0*(3.0*alpha+1.0)*taxa->get_brnlen())+ 
			exp(-1.0*(3.0*beta+1.0)*taxa->get_brnlen())/(3.0*(3.0*beta+1)));
		taxa->set_trpb(rate_num, snpstate_to_snp(TYPE_B), snpstate_to_snp(HETEROZYGOUS), 
			beta/(3.0*beta+1.0)-1.0/3.0*exp(-1.0*(3.0*alpha+1.0)*taxa->get_brnlen())+ 
			exp(-1.0*(3.0*beta+1.0)*taxa->get_brnlen())/(3.0*(3.0*beta+1)));

		taxa->set_trpb(rate_num, snpstate_to_snp(TYPE_A), snpstate_to_snp(SNP_ABSENT), 
			1.0/(3.0*beta+1.0)*(1.0-exp(-1.0*(3.0*beta+1.0)*taxa->get_brnlen())));
		taxa->set_trpb(rate_num, snpstate_to_snp(TYPE_B), snpstate_to_snp(SNP_ABSENT), 
			1.0/(3.0*beta+1.0)*(1.0-exp(-1.0*(3.0*beta+1.0)*taxa->get_brnlen())));
		taxa->set_trpb(rate_num, snpstate_to_snp(HETEROZYGOUS), snpstate_to_snp(SNP_ABSENT), 
			1.0/(3.0*beta+1.0)*(1.0-exp(-1.0*(3.0*beta+1.0)*taxa->get_brnlen())));

		taxa->set_trpb(rate_num, snpstate_to_snp(HETEROZYGOUS), snpstate_to_snp(HETEROZYGOUS),
			beta/(3.0*beta+1.0)+2.0/3.0*exp(-1.0*(3.0*alpha+1.0)*taxa->get_brnlen())+
			exp(-1.0*(3.0*beta+1.0)*taxa->get_brnlen())/(3.0*(3.0*beta+1)));

		taxa->set_trpb(rate_num, snpstate_to_snp(HETEROZYGOUS), snpstate_to_snp(TYPE_A), 
				beta/(3.0*beta+1.0)-1.0/3.0*exp(-1.0*(3.0*alpha+1.0)*taxa->get_brnlen())+
				exp(-1.0*(3.0*beta+1.0)*taxa->get_brnlen())/(3.0*(3.0*beta+1)));
			
		taxa->set_trpb(rate_num, snpstate_to_snp(HETEROZYGOUS), snpstate_to_snp(TYPE_B), 
				beta/(3.0*beta+1.0)-1.0/3.0*exp(-1.0*(3.0*alpha+1.0)*taxa->get_brnlen())+
				exp(-1.0*(3.0*beta+1.0)*taxa->get_brnlen())/(3.0*(3.0*beta+1)));

		taxa->set_trpb(rate_num, snpstate_to_snp(SNP_ABSENT), snpstate_to_snp(SNP_ABSENT),
				1.0/(3.0*beta+1)*(1.0+3.0*beta*exp(-1.0*(3.0*beta+1)*taxa->get_brnlen())));

		taxa->set_trpb(rate_num, snpstate_to_snp(SNP_ABSENT), snpstate_to_snp(TYPE_A),
			beta/(3.0*beta+1)*(1.0-exp(-1.0*(3.0*beta+1.0)*taxa->get_brnlen())));
		taxa->set_trpb(rate_num, snpstate_to_snp(SNP_ABSENT), snpstate_to_snp(TYPE_B),
			beta/(3.0*beta+1)*(1.0-exp(-1.0*(3.0*beta+1.0)*taxa->get_brnlen())));
		taxa->set_trpb(rate_num, snpstate_to_snp(SNP_ABSENT), snpstate_to_snp(HETEROZYGOUS),
			beta/(3.0*beta+1)*(1.0-exp(-1.0*(3.0*beta+1.0)*taxa->get_brnlen())));
	}
		
}
	

double SNP_model::find_ut(Branch *taxa)
{
	return(taxa->expect_subs_site());
}



void SNP_model::store_conprobs()
{
	int i, j, k, l;
	double prob, lnL=0;
	
	
#ifdef MPI_CODE_VERSION
	for(i=mpi_info.start_site; i<mpi_info.end_site; i++) 
#else
	for(i=0; i<curr_exchange->get_num_localities(); i++)
#endif
	{
		prob=prob_w_rate(i,0);
		for(k=0; k<curr_exchange->get_num_branches(); k++) {
			if ((*curr_tree)[k]->get_save_conprob_index() != -1) {
				for(j=0; j<curr_exchange->get_condlike_size(); j++)
					save_conprobs[(*curr_tree)[k]->get_save_conprob_index()][i][j]=(*curr_tree)[k]->get_cond_prob(j);
			}
		}
		lnL+=(log(prob)-log(LIKE_SCALE));
	}
	
	cout<<"PreBranch lnL: "<<lnL<<endl;
}



double SNP_model::newton_branch_opt(double tol)
{
	//The name is wrong since we aren't using Newton's method here,
	//but it hooks in properly with the Powell code
	int i,j, conprobs_index;
	long pretime, posttime;
	double low_t, guess_t, actual_t, low_lnL, guess_lnL, actual_lnL, final_min, last_lnL, old_brn, old_lnL;
	double (*thefunc)(double);
	Branch *decend_brn;

	thefunc=&snp_brn_min;
	
	for(i=0; i<curr_exchange->get_num_branches(); i++)
	{
		if ((*curr_tree)[i] != curr_tree->find_root()) {
			if ((swap_brn_opt==(BOOL)FALSE) || (i%10 ==last_brn_remain)) {
				brn_op_cbrn=(*curr_tree)[i];

				for(j=0; j<curr_exchange->get_num_branches(); j++)
					(*curr_tree)[j]->set_save_conprob_index(-1);

				conprobs_index=2;
				//Start after first two scratch arrays
				brn_op_cbrn->set_save_conprob_index(conprobs_index++);
				brn_op_cbrn->get_sibling()->set_save_conprob_index(conprobs_index++);
				decend_brn=brn_op_cbrn->get_parent();

				while (decend_brn != curr_tree->find_root()) {
					decend_brn->get_sibling()->set_save_conprob_index(conprobs_index++);
					decend_brn=decend_brn->get_parent();
				}
				//decend_brn->get_sibling()->set_save_conprob_index(conprobs_index);
#ifdef DO_TIME
				//pretime=RunTime();
#endif
				
				store_conprobs();
#ifdef DO_TIME
				//posttime=RunTime();
				cout<<"Store conprobs time: "<<posttime-pretime<<endl;
#endif
				old_brn=brn_op_cbrn->get_brnlen();
				old_lnL=optimize_a_branch(old_brn);
				low_t=1e-10;
				guess_t=0.5;
#ifdef DO_TIME
				//pretime=RunTime();
#endif
			
				mnbrak(&low_t, &guess_t, &actual_t, &low_lnL, &guess_lnL, &actual_lnL, thefunc);
				last_lnL =actual_lnL+10;
				while (last_lnL-actual_lnL >FLOAT_TOL) {
					last_lnL=actual_lnL;
					actual_lnL=brent(low_t, actual_t, guess_t, thefunc, tol, &final_min); 
	   
				}
#ifdef DO_TIME
				//posttime=RunTime();
				cout<<"Brn opt time: "<<posttime-pretime<<endl;
#endif
				if (actual_lnL < old_lnL)
					brn_op_cbrn->set_brnlen(fabs(final_min));
				else
					brn_op_cbrn->set_brnlen(old_brn);
				calc_transprobs(brn_op_cbrn, 0);
				cout<<"Opt val of brn: "<<i<<" is "<<brn_op_cbrn->get_brnlen()<<" lnL: "<<actual_lnL<<endl;
			}
		}
	}
	last_brn_remain++;
	if (last_brn_remain == 10) last_brn_remain=0;
	return(actual_lnL);
}



double SNP_model::optimize_a_branch(double brnlen)
{
	//When optimizing a branch, we only need to recalculate on a path down from
	//the branch in question to the root, using the saved conditional probs
	//for the rest of the tree
	int i, j, k;
	double lnL=0, lnL_s=0, temp_p1, temp_p2;
	Branch *curr_branch, *sibling;

	curr_branch=brn_op_cbrn;
	curr_branch->set_brnlen(fabs(brnlen));
	calc_transprobs(curr_branch, 0);
	
	//Move in the initial condition probabilities
	if (curr_branch->is_tip() == (BOOL)FALSE) {
#ifdef MPI_CODE_VERSION
		for(i=mpi_info.start_site; i<mpi_info.end_site; i++)
#else
		for(i=0; i<curr_exchange->get_num_localities(); i++)
#endif
			for(j=0; j<curr_exchange->get_condlike_size(); j++)
				save_conprobs[0][i][j]=save_conprobs[curr_branch->get_save_conprob_index()][i][j];
	}

	while (curr_branch != curr_tree->find_root()) {
		sibling=curr_branch->get_sibling();
		
		if (curr_branch->is_tip() == (BOOL)TRUE) {
			if (sibling->is_tip() == (BOOL)TRUE) {
#ifdef MPI_CODE_VERSION
				for(i=mpi_info.start_site; i<mpi_info.end_site; i++) {
#else
				for(i=0; i<curr_exchange->get_num_localities(); i++) {
#endif
					for(j=0; j<curr_exchange->get_condlike_size(); j++)
							save_conprobs[1][i][j]= (long double)LIKE_SCALE*
							get_tip_prob(i, curr_branch, 0, j)*get_tip_prob(i, sibling, 0, j);
					}
			}
			else {
#ifdef MPI_CODE_VERSION
				for(i=mpi_info.start_site; i<mpi_info.end_site; i++) {
#else
				for(i=0; i<curr_exchange->get_num_localities(); i++) {
#endif
					for(j=0; j<curr_exchange->get_condlike_size(); j++) {
						temp_p1=0;
				
						for(k=0; k<curr_exchange->get_condlike_size(); k++)
							temp_p1+=save_conprobs[sibling->get_save_conprob_index()][i][k]*sibling->get_trpb(0, j, k); 
			  
						save_conprobs[1][i][j]=
							get_tip_prob(i, curr_branch, 0, j)*temp_p1;  
					}   
				}
			}
  		}
		
		else {
			//Internal branch
			if (sibling->is_tip()==(BOOL)FALSE) 	{
#ifdef MPI_CODE_VERSION
				for(i=mpi_info.start_site; i<mpi_info.end_site; i++) {
#else
				for(i=0; i<curr_exchange->get_num_localities(); i++){
#endif
					for(j=0; j<curr_exchange->get_condlike_size(); j++) {			
						temp_p1=temp_p2=0;
						for(k=0; k<curr_exchange->get_condlike_size(); k++)
						{
							temp_p1+=save_conprobs[0][i][k]*curr_branch->get_trpb(0, j, k);
							temp_p2+=save_conprobs[sibling->get_save_conprob_index()][i][k]*sibling->get_trpb(0, j, k);
						} 
						save_conprobs[1][i][j]=(temp_p1/SQRT_LIKE_SCALE)*(temp_p2/SQRT_LIKE_SCALE);
					}	  
				}
			}
    
			else
			{
#ifdef MPI_CODE_VERSION
				for(i=mpi_info.start_site; i<mpi_info.end_site; i++) {
#else
					for(i=0; i<curr_exchange->get_num_localities(); i++) {
#endif
					for(j=0; j<curr_exchange->get_condlike_size(); j++) {			
						temp_p1=0;
						for(k=0; k<curr_exchange->get_condlike_size(); k++) 
							temp_p1+=save_conprobs[0][i][k]*curr_branch->get_trpb(0, j, k);
			  
						save_conprobs[1][i][j]= temp_p1*get_tip_prob(i, sibling, 0, j); 
					}
				}
			}

		}

		curr_branch=curr_branch->get_parent();

		//Rest cond probs for next loop
#ifdef MPI_CODE_VERSION
		for(i=mpi_info.start_site; i<mpi_info.end_site; i++)
#else
		for(i=0; i<curr_exchange->get_num_localities(); i++)
#endif
			for(j=0; j<curr_exchange->get_condlike_size(); j++)
				save_conprobs[0][i][j]=save_conprobs[1][i][j];
	} 

#ifdef MPI_CODE_VERSION
		for(i=mpi_info.start_site; i<mpi_info.end_site; i++) {
#else
		for(i=0; i<curr_exchange->get_num_localities(); i++){
#endif
		temp_p1=0;
		for (j=0; j<curr_exchange->get_condlike_size(); j++) {
			for (k=0; k<curr_exchange->get_condlike_size();k++) {
				temp_p1+=root_freq(k)*curr_branch->get_trpb(0, k, j)*
				save_conprobs[1][i][j];
			}
		}

		if (curr_exchange->likelihood_is_scaled() == (BOOL)FALSE)
			lnL+=log(temp_p1);
		else
			lnL+=(log(temp_p1)-log(LIKE_SCALE));
	}
#ifdef MPI_CODE_VERSION
	MPI_Allreduce(&lnL, &lnL_s, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			lnL=lnL_s;
#endif
	return(-lnL);
}
