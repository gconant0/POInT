//Copyright 2005-2014 Gavin Conant

#include <iostream>
#include <math.h>
#include <iomanip>
#include <fstream>
#include "genome_ploidy_like.h"
#include "powell.h"

#ifdef _OPEN_MP_VERSION_
#include "omp.h"
#endif



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

int cmpWGX(const void *x, const void *y)
{
  double xx = *(double*)x, yy = *(double*)y;
  if (xx < yy) return -1;
  if (xx > yy) return  1;
  return 0;
}


Ploidy_Like_model::Ploidy_Like_model() : Like_model()
{
    switch_permutes=0;
    tracking_state_lookups=0;
    std::cerr<<"Error: call to default constructor of Ploidy_Like_model\n";
    has_degen=0;
    gene_pos=0;
}


Ploidy_Like_model::Ploidy_Like_model(Exchange *cexchange, Tree *ctree, Clade *cgenomes, WGX_Data *chomologs, Phylo_Matrix *cmatrix)
{
    switch_permutes=0;
    tracking_state_lookups=0;
    the_matrix=cmatrix;
    allocate_state_model(cexchange,ctree, cgenomes, chomologs);
}

double Ploidy_Like_model::root_freq(int site)
{
    if (the_matrix->get_nth_state(site)->is_root_state() == TRUE)
        return(1.0);
    else
        return(0.0);
}



void Ploidy_Like_model::partial_prob_w_rate(int locale, Branch *lsib, int rate_num, Branch *stop_id)
 //The heart of the likelihood calculation: uses the tranisition probablity matrices stored in
  //the tree to calculate the likelihood of part of that tree.  Uses post-order tree transversal.
  //This is modified version for the Dupl_models that allows the tip states to be ambiguous--hence
  //probability is summed across those ambiguous states
{
  int i, j;
  double temp_p1, temp_p2;
  Branch *rsib, *parent;
  
  rsib=lsib->get_sibling();
  parent=lsib->get_parent();
  

  
  if (lsib->is_tip()==(BOOL)FALSE)
    {
      if (rsib->is_tip()==(BOOL)FALSE) 
		{
			partial_prob_w_rate(locale, curr_tree->find_left_tip(rsib), rate_num, rsib);
			for(i=0; i<the_matrix->get_num_states(); i++)
			{
				temp_p1=temp_p2=0;
				for(j=0; j<the_matrix->get_num_states(); j++)
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
			for(i=0; i<the_matrix->get_num_states(); i++)
			{
				temp_p1=0;
				for(j=0; j<the_matrix->get_num_states(); j++)
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
	
			for(i=0; i<the_matrix->get_num_states(); i++)
			{
				temp_p2=0;
				
				for(j=0; j<the_matrix->get_num_states(); j++)
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
          for(i=0; i<the_matrix->get_num_states(); i++) {
#ifdef _OPEN_MP_VERSION_
			parent->set_cond_prob_locale(omp_get_thread_num(), i, get_tip_prob(locale, lsib, rate_num, i)*get_tip_prob(locale, rsib, rate_num, i));
#else
             // cout<<the_matrix->get_nth_state(i)->get_state_name()<<": "<<get_tip_prob(locale, lsib, rate_num, i)<<"* "<<get_tip_prob(locale, rsib, rate_num, i)<<endl;
            parent->set_cond_prob(i, get_tip_prob(locale, lsib, rate_num, i)*get_tip_prob(locale, rsib, rate_num, i));
#endif
          }
	  
    }

  if (parent!=stop_id)
      partial_prob_w_rate(locale, parent, rate_num, stop_id);
 
}  //End Like_model::partial_prob_w_rate





double Ploidy_Like_model::get_tip_prob(int locale, Branch *taxa, int rate_num, int cond_state)
{
    int i;
    double retval=0;
    //cout<<"At tip, base state is "<<the_matrix->get_nth_state((*curr_data)[taxa->get_taxa_id()][locale])->get_state_name()<<endl;
    //cout<<" Num redund for this state "<<the_matrix->get_nth_state((*curr_data)[taxa->get_taxa_id()][locale])->num_state_redunds()<<endl;
    for(i=0; i<the_matrix->get_nth_state((*curr_data)[taxa->get_taxa_id()][locale])->num_state_redunds(); i++) {
        retval += taxa->get_trpb(rate_num, cond_state, the_matrix->get_nth_state((*curr_data)[taxa->get_taxa_id()][locale])->get_ith_redund_state(i));
    //return(taxa->get_trpb(rate_num, cond_state, the_matrix->get_nth_state((*curr_data)[taxa->get_taxa_id()][locale])->get_state_id()));
    //cout<<"Site num"<<locale<<" and "<<rate_num<<" and branch "<<taxa->get_taxa_id()<<" Tib prob from ith redund "<<i<<" e.g. "<<the_matrix->get_nth_state(cond_state)->get_state_name()<<" to tip state "<<the_matrix->get_nth_state(the_matrix->get_nth_state((*curr_data)[taxa->get_taxa_id()][locale])->get_ith_redund_state(i))->get_state_name()<<" is "<<taxa->get_trpb(rate_num, cond_state, the_matrix->get_nth_state(the_matrix->get_nth_state((*curr_data)[taxa->get_taxa_id()][locale])->get_ith_redund_state(i))->get_state_id())<<endl;
    }
    return(retval);
}



double Ploidy_Like_model::find_appropriate_ln_like()
{
	int i, j, k, l, num_prob_indices;
#ifdef MAC_TIME
	uint64_t  pretime, pptime, walktime, posttime, preptime, prepaftertime, siteaftertime, makelnltime, makelnlaftertime, sitetotal=0;
#endif
#ifdef DO_TIME
	unsigned long pretime, pptime, walktime, posttime, preptime, prepaftertime, siteaftertime, makelnltime, makelnlaftertime, sitetotal=0;
#endif
	double lnL, sum_prob, site_prob, trans_prob;

	
#ifdef MAC_TIME
	pretime= mach_absolute_time();
#endif
	
#ifdef DO_TIME
	pretime= RunTime();
#endif
	
	//if(the_tracks->tracking_current() ==(BOOL)FALSE)
	//	the_tracks->update_tracking();
	
    for(i=0; i<curr_exchange->get_num_branches(); i++) 
    
	calc_pattern_probs();

#ifdef MAC_TIME
	pptime= mach_absolute_time();
    std::cout<<"Pattern prob calc time: "<<pptime-pretime<<endl;
#endif	
	
#ifdef DO_TIME
	pptime= RunTime();
	std::cout<<"Pattern prob calc time: "<<pptime-pretime<<endl;
#endif	
//	cout<<"0\t";
    
    
#if 0
    for(k=0; k<curr_exchange->get_num_branches(); k++) {
        cout<<"Branch: "<<(*curr_tree)[k]->get_name()<<endl;
        for(i=0; i<the_matrix->get_num_states(); i++) {
            for(j=0; j<the_matrix->get_num_states(); j++)
                cout<<(*curr_tree)[k]->get_trpb(0, i, j)<<"\t";
            cout<<endl;
        }
        
    }

#endif
   
    set_site_base_states(0);

#ifdef _OPEN_MP_VERSION_
#pragma omp parallel for private (i, j, site_prob)
#endif
    for(i=0; i<the_tracks->get_num_possible_track_orders(); i++)
        #ifdef _OPEN_MP_VERSION_
        last_cumulative_probs[i]=pattern_probs[site_dupl_indices[0]][get_tracked_pattern_subindex_from_state(site_dupl_indices[0], i, taxa_tracks_thread[omp_get_thread_num()])];
#else
        last_cumulative_probs[i]=pattern_probs[site_dupl_indices[0]][get_tracked_pattern_subindex_from_state(site_dupl_indices[0], i, taxa_tracks)];
#endif

    
	//std::cout<<"Site 0 Check prob "<<last_cumulative_probs[0]<<" and 32: "<<last_cumulative_probs[31]<<endl;

	for(i=1; i<the_homologs->get_num_homologs(); i++) {
		//cout<<i<<"\t";
		set_site_base_states(i);
        build_transition_vector(i);
#ifdef MAC_TIME
		prepaftertime= mach_absolute_time();
		std::cout<<"Time for prep: "<<prepaftertime-preptime<<endl;
#endif
		
#ifdef DO_TIME
		prepaftertime= RunTime();
		std::cout<<"Time for prep: "<<prepaftertime-preptime<<endl;
#endif

#ifdef _OPEN_MP_VERSION_
#pragma omp parallel for private (j, k, l, sum_prob, site_prob, trans_prob)
#endif
		for(j=0; j<the_tracks->get_num_possible_track_orders(); j++) {
#ifdef _OPEN_MP_VERSION_
            site_prob=pattern_probs[site_dupl_indices[i]][get_tracked_pattern_subindex_from_state(site_dupl_indices[i], j, taxa_tracks_thread[omp_get_thread_num()])];
#else
            site_prob=pattern_probs[site_dupl_indices[i]][get_tracked_pattern_subindex_from_state(site_dupl_indices[i], j, taxa_tracks)];
#endif
           // cout<<"Site "<<i<<" Track "<<j<<" Dupl index: "<<site_dupl_indices[i]<<" is subpattern: "<<get_tracked_pattern_subindex_from_state(site_dupl_indices[i], j, taxa_tracks);
            //for(l=0; l<curr_exchange->get_num_taxa(); l++) cout<<" Taxa "<<l<<" track# "<<taxa_tracks[l]<<"\t";
            //cout<<" Prob: "<<site_prob<<endl;
            
            for(k=0; k<the_tracks->get_num_possible_track_orders(); k++) {
#ifdef _OPEN_MP_VERSION_
                set_track_states(k, taxa_tracks_thread2[omp_get_thread_num()]);
#else
                set_track_states(k, taxa_tracks2);
#endif
                trans_prob=1.0;
                //for(l=0; l<curr_exchange->get_num_taxa(); l++) cout<<" Taxa "<<l<<" track2# "<<taxa_tracks2[l]<<"\t";
#ifdef _OPEN_MP_VERSION_
                for(l=0; l<curr_exchange->get_num_taxa(); l++) trans_prob *= track_transprobs[l][taxa_tracks_thread2[omp_get_thread_num()][l]][taxa_tracks_thread[omp_get_thread_num()][l]];
#else
                for(l=0; l<curr_exchange->get_num_taxa(); l++) trans_prob *= track_transprobs[l][taxa_tracks2[l]][taxa_tracks[l]];
#endif
               //cout<<"S"<<i<<": "<<j<<"->"<<k<<": "<<trans_prob<<endl;
#ifdef _OPEN_MP_VERSION_
				sum_probs_thread[omp_get_thread_num()][k]=log(trans_prob)+last_cumulative_probs[k];
            }
			sum_prob=get_ln_sum_prob_sort(sum_probs_thread[omp_get_thread_num()]);
#else
				sum_probs[k]=log(trans_prob)+last_cumulative_probs[k];
            }
			sum_prob=get_ln_sum_prob_sort(sum_probs);
#endif		
        //cout<<i<<" T:"<<j<<": "<<site_prob<<"->"<<sum_prob<<endl;
			//CHECK THIS!!!!!
			cumulative_probs[j] = site_prob + sum_prob;
        }
		//std::cout<<"Site "<<i<<" Check prob "<<cumulative_probs[0]<<endl;
#ifdef MAC_TIME
		siteaftertime= mach_absolute_time();
		sitetotal+=siteaftertime-prepaftertime;
		//std::cout<<"Time for one site calc: "<<i<<": "<<siteaftertime-prepaftertime<<endl;
#endif
		
#ifdef DO_TIME
		siteaftertime= RunTime();
		sitetotal+=siteaftertime-prepaftertime;
		std::cout<<"Time for one site calc: "<<i<<": "<<siteaftertime-prepaftertime<<endl;
#endif
		
#ifdef _OPEN_MP_VERSION_
#pragma omp parallel for private (j)
#endif
		for(j=0; j<the_tracks->get_num_possible_track_orders(); j++)
			last_cumulative_probs[j]=cumulative_probs[j];

	}
	

		lnL=get_ln_sum_prob_sort(cumulative_probs);
    //cout<<"ln likelihood: "<<lnL<<endl;

#ifdef MAC_TIME
	posttime=mach_absolute_time();
	std::cout<<"Runtime for likelihood calc: "<<posttime-pretime<<endl;
	std::cout<<"Time in main loop: "<<sitetotal<<endl;
#endif	
	
#ifdef DO_TIME
	posttime=RunTime();
	std::cout<<"Runtime for likelihood calc: "<<posttime-pretime<<endl;
	std::cout<<"Time in main loop: "<<sitetotal<<endl;
#endif	
	
	return(lnL);
}


void Ploidy_Like_model::describe_results()
{
    int i;
    
    for(i=0; i<the_matrix->get_num_params(); i++)
        cout<<"Parameter "<<i<<" ("<<the_matrix->get_param_name(i)<<")= "<<the_matrix->get_param(i)<<endl;
}



double Ploidy_Like_model::get_post_prob(int site, int pattern)
{
	if ((pattern >=0) &&(pattern < the_tracks->get_num_possible_track_orders()))
		return(post_probs[site][pattern]);
	else
		return(post_probs[site][0]);

}




void Ploidy_Like_model::print_tracking_probs(string outfile)
{
	int i, j, k, l, *taxa_vals;
	ofstream fout;

	get_gene_conditional_probs();
    taxa_vals=new int[curr_exchange->get_num_taxa()];
    
	fout.open(outfile.c_str());
	if (!fout.fail()) {

        for(i=0; i<curr_exchange->get_num_taxa(); i++) {
            for(j=0; j<the_homologs->get_dupl_level(); j++)
                fout<<(*the_tracks)[i].get_genome()->get_name()<<"\t";
        }
		fout<<"#";
        for(j=0; j<the_tracks->get_num_possible_track_orders(); j++) {
            set_track_states(j,taxa_vals);
            for(i=0; i<curr_exchange->get_num_taxa(); i++) {
                fout<<(*the_tracks)[i].get_genome()->get_name();
                  for(k=0; k<the_homologs->get_dupl_level(); k++) {
                        fout<<tracking_permutes[taxa_vals[i]][k];
                   
                  }
                if (i != curr_exchange->get_num_taxa()-1)
                    fout<<"_";
            }
            if (j != the_tracks->get_num_possible_track_orders())
                fout<<"\t#";
        }

        fout<<"\n";

            
		for (i=0; i<the_homologs->get_num_homologs(); i++) {
			for(j=0; j<curr_exchange->get_num_taxa(); j++) {
                for(k=0; k<the_homologs->get_dupl_level(); k++) {
                    if ((*the_tracks)[j].get_gene_track(i, k)->my_locus != 0)
						fout<<(*the_genomes)[j][(*the_tracks)[j].get_gene_track(i, k)->
						my_locus->get_contig((*the_tracks)[j].get_gene_track(i, k)->index_num)]
					[(*the_tracks)[j].get_gene_track(i, k)->
						my_locus->get_gene((*the_tracks)[j].get_gene_track(i, k)->index_num)].get_name()<<"\t";
				else
					fout<<"NONE\t";
                }
            }
            
            for(j=0; j<the_tracks->get_num_possible_track_orders(); j++) {
				if (j != the_tracks->get_num_possible_track_orders()-1)
					fout<<post_probs[i][j]<<"\t";
				else
					fout<<post_probs[i][j]<<"\n";
            }
        }
        fout.close();
    }
    delete[] taxa_vals;
}



void Ploidy_Like_model::print_transpoint_branch_probs(string filename)
{
    int i, m,start_state, end_state;
    Branch *non_zero_brn;
    ofstream fout;
    
    i=0;
    while((*curr_tree)[i]->get_brnlen() ==0) i++;
    non_zero_brn=(*curr_tree)[i];
    
    calc_transpoint_branch_probs();
    
    
    fout.open(filename.c_str());
    fout<<"Site\t";
    for(m=0; m<curr_exchange->get_num_branches(); m++) {
        fout<<(*curr_tree)[m]->get_name()<<"\t";
        fout<<(*curr_tree)[m]->expect_subs_site()<<"\t";
    }
    fout<<endl;
    fout<<"Site";
    
    for(m=0; m<curr_exchange->get_num_branches(); m++) {
        for(start_state=0; start_state<the_matrix->get_num_states(); start_state++) {
            for(end_state=0; end_state<the_matrix->get_num_states(); end_state++) {
                if (non_zero_brn->get_trpb(0, start_state, end_state) !=0.0)
                    fout<<"\t"<<the_matrix->get_nth_state(start_state)->get_state_name()<<"-->"<<the_matrix->get_nth_state(end_state)->get_state_name();
            }
        }
    }
    fout<<endl;
    
    for(i=0; i<the_homologs->get_num_homologs(); i++) {
        fout<<i;
        for(m=0; m<curr_exchange->get_num_branches(); m++) {
            for(start_state=0; start_state<the_matrix->get_num_states(); start_state++) {
                for(end_state=0; end_state<the_matrix->get_num_states(); end_state++) {
                    if (non_zero_brn->get_trpb(0, start_state, end_state) !=0.0) {
                        fout<<"\t"<<branch_state_probs[i][m][start_state][end_state];
                    }
                }
                
            }
        }
        fout<<endl;
    }
    
    
    fout.close();	
}






//NOT FINISHED!!!
void Ploidy_Like_model::allocate_state_model (Exchange *cexchange, Tree *ctree, Clade *cgenomes, WGX_Data *chomologs)
{
	int i, j, k, curr, pattern_pos, *elements, divisor_so_far;

	if (cexchange->get_strand_switch_prob() == 0.0)
		cexchange->set_strand_switch_prob(0.01);
	curr_exchange=cexchange;
    

#ifdef _OPEN_MP_VERSION_
	curr_data = new Sequence_dataset(cexchange->get_num_taxa(), cexchange->get_num_open_mp_threads(), ARBITRARY);
#else
	curr_data = new Sequence_dataset(cexchange->get_num_taxa(), 1, ARBITRARY);
#endif
	
	//std::cout<<"Curr data has size "<<(*curr_data)[0].Sequence_size()<<" sites"<<" : Threads by exchange: "<<cexchange->get_num_open_mp_threads()<<endl;
	
    switch_prob_param_num=the_matrix->get_num_params()-1;
    
	the_genomes=cgenomes;
	the_homologs=chomologs;
    
   
    has_degen=new BOOL [the_genomes->get_num_genomes()];
    gene_pos=new int [the_genomes->get_num_genomes()];
        //for (j=0; j<cexchange->get_num_taxa(); j++)
          //  cout<<"Site "<<9073<<" taxa "<<j<<" num dupls: "<<(*the_homologs)[9073][j].get_dupl_count()<<endl;
   
	the_tracks = new WGX_Tracks(chomologs, cgenomes);
    
    
    num_zero_brns=0;
    if (curr_exchange->zero_len_brns_fixed() == TRUE) {
        for(i=0; i<curr_exchange->get_num_branches(); i++) {
            if (curr_exchange->zero_len_brns_fixed() ==(BOOL)TRUE) {
                if ((*ctree)[i]->expect_subs_site()  == 0) num_zero_brns++;
            }
        }
        cout<<"Fixing "<<num_zero_brns<<" branches to have 0 length\n";
       // brn_index=new int[curr_exchange->get_num_branches()-num_zero_brns];
    }
    
	site_probs=0;
	do_transpoint_probs=(BOOL)FALSE;

	powN[0]=1;
    pow2[0]=1;
    for(i=1; i<MAX_TAXA; i++) {
		powN[i] =powN[i-1]*the_homologs->get_dupl_level();
        pow2[i]=pow2[i-1]*2;
    }
    
    num_dupl_patterns=powN[cexchange->get_num_taxa()];
    total_pattern_states=new int [num_dupl_patterns];
    site_states=new Model_State * [cexchange->get_num_taxa()];

    num_level_patterns=new int [the_homologs->get_dupl_level()+1];
    
    for(i=0; i<=the_homologs->get_dupl_level(); i++)
        num_level_patterns[i]=n_choose_k(the_homologs->get_dupl_level(), i);
    
    //Level 0 is all lost.  We dont use it, but here for completeness
        
    taxa_dupl_pattern_states = new int * [num_dupl_patterns];
    for(i=0; i<num_dupl_patterns; i++)
        taxa_dupl_pattern_states[i]=new int [cexchange->get_num_taxa()];
	
    
    for(i=0; i<num_dupl_patterns; i++) {
        curr=i;
        total_pattern_states[i]=1;
        for(j=0; j<cexchange->get_num_taxa(); j++) {
            //Num of copies present for the jth taxa is the modulus of the jth order power of i
            taxa_dupl_pattern_states[i][j] = (curr%the_homologs->get_dupl_level())  +1;
            
            //cout<<"For dupl_state "<<i<<" taxa "<<j<<" has "<<taxa_dupl_pattern_states[i][j]<<" copies."<<" This implies "<<num_level_patterns[taxa_dupl_pattern_states[i][j]]<<"states for the taxa\n";
            
            total_pattern_states[i]*=num_level_patterns[taxa_dupl_pattern_states[i][j]];
            
            //Prepare for the next taxa by removing the lowest order factor that we just stored
            curr = curr/the_homologs->get_dupl_level();
        }
       // cout<<"Total substates for dupllevel "<<i<<" is "<<total_pattern_states[i]<<endl;
    }

    //pattern_counters=new int [cexchange->get_num_taxa()];
    
    pattern_substates=new int ** [num_dupl_patterns];
    
    for(i=0; i<num_dupl_patterns; i++) {
        pattern_substates[i]=new int * [cexchange->get_num_taxa()];
    
        for(j=0; j<cexchange->get_num_taxa(); j++)
            pattern_substates[i][j] = new int [total_pattern_states[i]];
        
        divisor_so_far=1;
        
        for (j=0; j<cexchange->get_num_taxa(); j++) {
            for(k=0; k<total_pattern_states[i]; k++) {
                pattern_substates[i][j][k]= (k/divisor_so_far)%num_level_patterns[taxa_dupl_pattern_states[i][j]];
                //cout<<"Dupl level "<<i<<" taxa "<<j<<" sub index "<<k<<" yields "<<pattern_substates[i][j][k]<<endl;
            }
            divisor_so_far*=num_level_patterns[taxa_dupl_pattern_states[i][j]];
        }
        
        
    }
    
    elements=new int [the_homologs->get_dupl_level()];
    for(i=0; i<the_homologs->get_dupl_level(); i++) elements[i]=i;
    
    enum_permutations(the_homologs->get_dupl_level(), num_track_patterns, elements, tracking_permutes);
  
    cout<<"Set up permutation matrix with "<<num_track_patterns<<" permutations each of size "<<the_homologs->get_dupl_level()<<endl;
    
    delete[] elements;
    
    
    switch_permutes=new double [the_homologs->get_dupl_level()+1];
    for(i=0; i<the_homologs->get_dupl_level()+1; i++)
        switch_permutes[i] = 1.0/((double)recurse_factorial(i));
    
    taxa_tracks=new int [cexchange->get_num_taxa()];
    taxa_tracks2=new int [cexchange->get_num_taxa()];

    
    track_transprobs=new double **[cexchange->get_num_taxa()];
    
    for(i=0; i<cexchange->get_num_taxa(); i++) {
        track_transprobs[i]=new double *[num_track_patterns];
        for(j=0; j<num_track_patterns; j++) track_transprobs[i][j]=new double[num_track_patterns];
    }

    
#ifdef _OPEN_MP_VERSION_
    taxa_tracks_thread = new int * [cexchange->get_num_open_mp_threads()];
    taxa_tracks_thread2 = new int * [cexchange->get_num_open_mp_threads()];
    for (i=0; i<cexchange->get_num_open_mp_threads(); i++) {
        taxa_tracks_thread[i]=new int [cexchange->get_num_taxa()];
        taxa_tracks_thread2[i]=new int [cexchange->get_num_taxa()];
    }
#endif
        
        //i=9073;
    //for (j=0; j<cexchange->get_num_taxa(); j++)
      //  cout<<"Prior to setup Site "<<i<<" taxa "<<j<<" num dupls: "<<(*the_homologs)[i][j].get_dupl_count()<<endl;
    
    setup_tracking_state_lookups(cexchange);
    
    
	state_trans_probs=new double [the_tracks->get_num_possible_track_orders()];
	cumulative_probs=new double [the_tracks->get_num_possible_track_orders()];
	last_cumulative_probs=new double [the_tracks->get_num_possible_track_orders()];
	sum_probs=new double [the_tracks->get_num_possible_track_orders()];
	track_switch_probs=new double [the_tracks->get_num_possible_track_orders()];
#ifdef _OPEN_MP_VERSION_
	sum_probs_thread=new double * [cexchange->get_num_open_mp_threads()];
	for (i=0; i<cexchange->get_num_open_mp_threads(); i++) 
		sum_probs_thread[i] = new double [the_tracks->get_num_possible_track_orders()];
#endif
	
	sum_probs_used=new BOOL [the_tracks->get_num_possible_track_orders()];

#ifdef _OPEN_MP_VERSION_
    sum_probs_used_thread=new BOOL * [cexchange->get_num_open_mp_threads()];
	for (i=0; i<cexchange->get_num_open_mp_threads(); i++) 
		sum_probs_used_thread[i]=new BOOL [the_tracks->get_num_possible_track_orders()];
#endif
	
	
	post_probs=0;
	branch_state_probs=0;
	
	
	//First index is assignment of dupl states per taxa
	pattern_probs=new double *[num_dupl_patterns];

	for(i=0; i<num_dupl_patterns; i++) {
	//Next index is number of arrangements of i duplicates in n taxa
		pattern_probs[i] = new double [total_pattern_states[i]];
	}
    

	assemble(cexchange, curr_data, ctree);
}


void Ploidy_Like_model::num_params_model()
{
    curr_exchange->set_num_params(the_matrix->get_num_params()+curr_exchange->get_num_branches()-num_zero_brns);
    std::cout<<"Setting num parameters to "<<curr_exchange->get_num_params()<<std::endl;
}


void Ploidy_Like_model::setup_tracking_state_lookups(Exchange *cexchange)
{
    int i, j, k, l, *lookup, *new_lookup, new_mask, curr;
    Model_State *my_state, *link_state;
    
    lookup=new int [the_homologs->get_dupl_level()];
    new_lookup=new int [the_homologs->get_dupl_level()];
    
    for(k=0; k<=the_homologs->get_dupl_level(); k++)
        the_matrix->initialize_all_cross_refs(k, num_track_patterns);
    
    for(i=0; i<the_matrix->get_num_states(); i++) {
    
        my_state=the_matrix->get_nth_state(i);
        for(k=0; k<the_homologs->get_dupl_level(); k++) {
            if ((my_state->get_binary_rep() & pow2[k]) != 0) lookup[k]=1;
            else lookup[k]=0;
        }
        //std::cout<<"My state is "<<my_state->get_state_name()<<" and my binary state is "<<my_state->get_binary_rep()<<std::endl;
        for(l=0; l<num_track_patterns; l++) {
            new_mask=0;
            for(k=0; k<the_homologs->get_dupl_level(); k++) {
                new_lookup[k]=lookup[tracking_permutes[l][k]];
            }
            for(k=0; k<the_homologs->get_dupl_level(); k++) {
                if (new_lookup[k] ==1) new_mask += pow2[k];
            }
           // std::cout<<"l: "<<l<<" ORIG lk: ";
            //for(k=0; k<the_homologs->get_dupl_level(); k++) cout<<lookup[k]<<"\t";
          //  cout<<" New lookup: ";
        //    for(k=0; k<the_homologs->get_dupl_level(); k++) cout<<new_lookup[k]<<"\t";
           // cout<<" Cross referencing state "<<my_state->get_state_name()<<" using mask "<<new_mask<<" to state "<<the_matrix->get_masked_state(new_mask)->get_state_name()<<std::endl;
            
            my_state->set_cross_ref(l, the_matrix->get_masked_state(new_mask));
        }
        
    }
    
    
    tracking_state_lookups=new Model_State ***[the_homologs->get_dupl_level()];
   
    
    for(i=0; i<the_homologs->get_dupl_level(); i++) {
        tracking_state_lookups[i]=new Model_State** [num_level_patterns[i]];
        for (j=0; j<num_level_patterns[i]; j++) {
            tracking_state_lookups[i][j]=new Model_State*[num_track_patterns];
        }
    }
    
    for(i=1; i<the_homologs->get_dupl_level(); i++) {
        for (j=0; j<num_level_patterns[i]; j++) {
            my_state=the_matrix->get_nth_level_ith_state(i, j);
            
            for(l=0; l<num_track_patterns; l++)
                tracking_state_lookups[i][j][l]=my_state->get_cross_ref(l);
                
            
        }
    }
    delete[] lookup;
    delete[] new_lookup;
    
    site_dupl_indices=new int [the_homologs->get_num_homologs()];
    //cout<<"Site\t";
    //for(j=0; j<curr_exchange->get_num_taxa(); j++)
     //   cout<<(*the_genomes)[j].get_name()<<"\t";
    //cout<<endl;
    
    for(i=0; i<the_homologs->get_num_homologs(); i++) {
        new_mask=0;
        //cout<<i<<"\t";
            for(j=0; j<cexchange->get_num_taxa(); j++) {
         //       cout<<(*the_homologs)[i][j].get_dupl_count()<<"\t";
                //Num of copies present for the jth taxa is the modulus of the jth order power of i
                curr = powN[j]*((*the_homologs)[i][j].get_dupl_count()-1);
                new_mask+=curr;
            }
        site_dupl_indices[i]=new_mask;
        
        
       // cout<<"\t"<<site_dupl_indices[i]<<endl;
    }

   
}





void Ploidy_Like_model::calc_pattern_probs()
{
	int i, j, l, start_state, end_state;
	
    
    //New code for new indexing order
    for(i=0; i<num_dupl_patterns; i++){
#ifdef _OPEN_MP_VERSION_
#pragma omp parallel for private(j, l, start_state, end_state)
#endif
        for(j=0; j<total_pattern_states[i]; j++) {

                setup_data(i, j);
#ifdef _OPEN_MP_VERSION_
                pattern_probs[i][j]=log(prob_w_rate(omp_get_thread_num(),0));
#else
                pattern_probs[i][j]=log(prob_w_rate(0,0));
#endif
            
            if (do_transpoint_probs == (BOOL)TRUE) {
                for(l=0; l<curr_exchange->get_num_branches(); l++) {
                    for (start_state=0; start_state<the_matrix->get_num_states(); start_state++) {
                        for(end_state=0; end_state<the_matrix->get_num_states(); end_state++) {
//#ifdef _OPEN_MP_VERSION_
//                        branch_transpoint_probs[i][j][l][start_state][end_state] =
//                        calculate_transpoint_prob((*curr_tree)[l], start_state, end_state, omp_get_thread_num());
//#else
                        branch_transpoint_probs[i][j][l][start_state][end_state] =
                        calculate_transpoint_prob((*curr_tree)[l], start_state, end_state, 0);
                            //cout<<"Pattern "<<i<<" subpattern "<<j<<" Branch "<<(*curr_tree)[l]->get_name()<<" has state "
                           // <<the_matrix->get_nth_state(start_state)->get_state_name()<<" to "<<the_matrix->get_nth_state(end_state)->get_state_name()<<" at prob "<< branch_transpoint_probs[i][j][l][start_state][end_state]<<endl;
//#endif
                        }

                    }
                }
                
            }

            //cout<<": "<<i<<", "<<j<<"Prob resulting "<<pattern_probs[i][j]<<endl;
        }

    }
}


void Ploidy_Like_model::get_gene_conditional_probs()
//Calculates the probabilities of various trackings for each "pillar" given those to the left and right
{
    int i, j, k, l, last_mask=0, this_mask;
    double sum, site_prob, sum_prob, trans_prob;

    
    post_probs=new double*[the_homologs->get_num_homologs()];
    left_cond_probs=new double*[the_homologs->get_num_homologs()];
    right_cond_probs=new double*[the_homologs->get_num_homologs()];

    
    for(i=0; i<the_homologs->get_num_homologs(); i++) {
        post_probs[i]=new double[the_tracks->get_num_possible_track_orders()];
        left_cond_probs[i]=new double[the_tracks->get_num_possible_track_orders()];
        right_cond_probs[i]=new double[the_tracks->get_num_possible_track_orders()];
    }
    
    
   // the_tracks->update_tracking();
    recalculate_transprobs();
    calc_pattern_probs();
   
    set_site_base_states(0);
    //cout<<"SITE 0 left: ";
    
#ifdef _OPEN_MP_VERSION_
#pragma omp parallel for private (i, j, site_prob)
#endif
    for(i=0; i<the_tracks->get_num_possible_track_orders(); i++) {
#ifdef _OPEN_MP_VERSION_
        left_cond_probs[0][i]=pattern_probs[site_dupl_indices[0]][get_tracked_pattern_subindex_from_state(site_dupl_indices[0], i, taxa_tracks_thread[omp_get_thread_num()])];
#else
        left_cond_probs[0][i]=pattern_probs[site_dupl_indices[0]][get_tracked_pattern_subindex_from_state(site_dupl_indices[0], i, taxa_tracks)];
#endif
        //cout<<left_cond_probs[0][i]<<"\t";
       // if (taxa_tracks[0] == taxa_tracks[1]) cout<<"Site "<<0<<" pattern: "<<i<<" left cond prob: "<<left_cond_probs[0][i]<<endl;
    }
    //cout<<endl;
    
    set_site_base_states(the_homologs->get_num_homologs()-1);
    
#ifdef _OPEN_MP_VERSION_
#pragma omp parallel for private (i, j, site_prob)
#endif
    for(i=0; i<the_tracks->get_num_possible_track_orders(); i++)
#ifdef _OPEN_MP_VERSION_
        right_cond_probs[the_homologs->get_num_homologs()-1][i]=pattern_probs[site_dupl_indices[the_homologs->get_num_homologs()-1]][get_tracked_pattern_subindex_from_state(site_dupl_indices[the_homologs->get_num_homologs()-1], i, taxa_tracks_thread[omp_get_thread_num()])];
#else
        right_cond_probs[the_homologs->get_num_homologs()-1][i]=pattern_probs[site_dupl_indices[the_homologs->get_num_homologs()-1]][get_tracked_pattern_subindex_from_state(site_dupl_indices[the_homologs->get_num_homologs()-1], i, taxa_tracks)];
#endif
    
    //std::cout<<"Site 0 Check prob "<<last_cumulative_probs[0]<<" and 32: "<<last_cumulative_probs[31]<<endl;
    
    
    for(i=1; i<the_homologs->get_num_homologs(); i++) {
        set_site_base_states(i);
        build_transition_vector(i);
        
        
#ifdef _OPEN_MP_VERSION_
#pragma omp parallel for private (j, k, l, sum_prob, site_prob, trans_prob)
#endif
        for(j=0; j<the_tracks->get_num_possible_track_orders(); j++) {
            left_cond_probs[i][j]=0;
#ifdef _OPEN_MP_VERSION_
            site_prob=pattern_probs[site_dupl_indices[i]][get_tracked_pattern_subindex_from_state(site_dupl_indices[i], j, taxa_tracks_thread[omp_get_thread_num()])];
#else
            site_prob=pattern_probs[site_dupl_indices[i]][get_tracked_pattern_subindex_from_state(site_dupl_indices[i], j, taxa_tracks)];
#endif
            //if (taxa_tracks[0] == taxa_tracks[1]) cout<<"Site "<<i<<" pattern: "<<j<<" prob: "<<site_prob<<endl;
            
            for(k=0; k<the_tracks->get_num_possible_track_orders(); k++) {
#ifdef _OPEN_MP_VERSION_
                set_track_states(k, taxa_tracks_thread2[omp_get_thread_num()]);
#else
                set_track_states(k, taxa_tracks2);
#endif
                trans_prob=1.0;
#ifdef _OPEN_MP_VERSION_
                for(l=0; l<curr_exchange->get_num_taxa(); l++) trans_prob *= track_transprobs[l][taxa_tracks_thread2[omp_get_thread_num()][l]][taxa_tracks_thread[omp_get_thread_num()][l]];
#else
                for(l=0; l<curr_exchange->get_num_taxa(); l++) trans_prob *= track_transprobs[l][taxa_tracks2[l]][taxa_tracks[l]];
#endif
#ifdef _OPEN_MP_VERSION_
                sum_probs_thread[omp_get_thread_num()][k]=log(trans_prob)+left_cond_probs[i-1][k];
            }
            sum_prob=get_ln_sum_prob_sort(sum_probs_thread[omp_get_thread_num()]);
#else
                sum_probs[k]=log(trans_prob)+left_cond_probs[i-1][k];
            }
            sum_prob=get_ln_sum_prob_sort(sum_probs);
#endif
            left_cond_probs[i][j] = site_prob+sum_prob;
        
        //if (taxa_tracks[0] == taxa_tracks[1]) cout<<"Site "<<i<<" pattern: "<<j<<" left cond prob: "<<left_cond_probs[i][j]<<" site = "<<site_prob<<" sums: "<<sum_prob<<endl;
        }

        
    }
    
    for(i=the_homologs->get_num_homologs()-2; i>0; i--) {
        set_site_base_states(i);
        build_forward_transition_vector(i);
        
        
#ifdef _OPEN_MP_VERSION_
#pragma omp parallel for private (j, k, l, sum_prob, site_prob, trans_prob)
#endif
        for(j=0; j<the_tracks->get_num_possible_track_orders(); j++) {
            right_cond_probs[i][j]=0;
#ifdef _OPEN_MP_VERSION_
            site_prob=pattern_probs[site_dupl_indices[i]][get_tracked_pattern_subindex_from_state(site_dupl_indices[i], j, taxa_tracks_thread[omp_get_thread_num()])];
#else
            site_prob=pattern_probs[site_dupl_indices[i]][get_tracked_pattern_subindex_from_state(site_dupl_indices[i], j, taxa_tracks)];
#endif
           // set_track_states(j, taxa_tracks);
            for(k=0; k<the_tracks->get_num_possible_track_orders(); k++) {
#ifdef _OPEN_MP_VERSION_
                set_track_states(k, taxa_tracks_thread2[omp_get_thread_num()]);
#else
                set_track_states(k, taxa_tracks2);
#endif
                trans_prob=1.0;
               
                //This should give the transprobs from state j in site i to state k in site i+1
#ifdef _OPEN_MP_VERSION_
                for(l=0; l<curr_exchange->get_num_taxa(); l++) trans_prob *= track_transprobs[l][taxa_tracks_thread2[omp_get_thread_num()][l]][taxa_tracks_thread[omp_get_thread_num()][l]];
#else
               for(l=0; l<curr_exchange->get_num_taxa(); l++) trans_prob *= track_transprobs[l][taxa_tracks2[l]][taxa_tracks[l]];
#endif
                
#ifdef _OPEN_MP_VERSION_
                sum_probs_thread[omp_get_thread_num()][k]=log(trans_prob)+right_cond_probs[i+1][k];
#else
                sum_probs[k]=log(trans_prob)+right_cond_probs[i+1][k];
#endif
            }
#ifdef _OPEN_MP_VERSION_
            sum_prob = get_ln_sum_prob_sort(sum_probs_thread[omp_get_thread_num()]);
#else
            sum_prob = get_ln_sum_prob_sort(sum_probs);
#endif
            
            right_cond_probs[i][j] = site_prob + sum_prob;
            //if ((i>9071) && (i < 9074)) {
//#pragma omp critical
            //if (taxa_tracks[0] == taxa_tracks[1]) cout<<"Site "<<i<<" pattern: "<<j<<" sub pattern: " <<get_tracked_pattern_subindex_from_state(site_dupl_indices[i], j, taxa_tracks)<<" right cond prob: "<<right_cond_probs[i][j]<<" site = "<<site_prob<<" sums: "<<sum_prob<<endl;
            //}
            
        }
        
    }
    

    set_site_base_states(0);
    build_forward_transition_vector(0);


#ifdef _OPEN_MP_VERSION_
#pragma omp parallel for private (j, k, site_prob, sum_prob, l, trans_prob)
#endif
    for(j=0; j<the_tracks->get_num_possible_track_orders(); j++) {
        post_probs[0][j]=0;
#ifdef _OPEN_MP_VERSION_
        site_prob=pattern_probs[site_dupl_indices[0]][get_tracked_pattern_subindex_from_state(site_dupl_indices[0], j, taxa_tracks_thread[omp_get_thread_num()])];
#else
        site_prob=pattern_probs[site_dupl_indices[0]][get_tracked_pattern_subindex_from_state(site_dupl_indices[0], j, taxa_tracks)];
#endif

        for(k=0; k<the_tracks->get_num_possible_track_orders(); k++) {
#ifdef _OPEN_MP_VERSION_
            set_track_states(k, taxa_tracks_thread2[omp_get_thread_num()]);
#else
            set_track_states(k, taxa_tracks2);
#endif
            trans_prob=1.0;
    
            //This should give the transprobs from state j in site i to state k in site i+1
#ifdef _OPEN_MP_VERSION_
            for(l=0; l<curr_exchange->get_num_taxa(); l++) trans_prob *= track_transprobs[l][taxa_tracks_thread2[omp_get_thread_num()][l]][taxa_tracks_thread[omp_get_thread_num()][l]];
#else
            for(l=0; l<curr_exchange->get_num_taxa(); l++) trans_prob *= track_transprobs[l][taxa_tracks2[l]][taxa_tracks[l]];
#endif
    
#ifdef _OPEN_MP_VERSION_
            sum_probs_thread[omp_get_thread_num()][k]=log(trans_prob)+right_cond_probs[1][k];
#else
            sum_probs[k]=log(trans_prob)+right_cond_probs[1][k];
#endif
        }
#ifdef _OPEN_MP_VERSION_
            sum_prob = get_ln_sum_prob_sort(sum_probs_thread[omp_get_thread_num()]);
#else
            sum_prob = get_ln_sum_prob_sort(sum_probs);
#endif


        post_probs[0][j] = site_prob + sum_prob;
        //if (taxa_tracks[0] == taxa_tracks[1])
     //   #pragma omp critical
            //cout<<"Site 0 pattern: "<<j<<"  post prob: "<<post_probs[0][j]<<" site = "<<site_prob<<" sums: "<<sum_prob<<endl;
        
    }


    
    for(i=1; i<the_homologs->get_num_homologs(); i++) {
        set_site_base_states(i);
        build_transition_vector(i);
        
#ifdef _OPEN_MP_VERSION_
#pragma omp parallel for private (j, k, l, sum_prob, site_prob, trans_prob)
#endif
        for(j=0; j<the_tracks->get_num_possible_track_orders(); j++) {
#ifdef _OPEN_MP_VERSION_
            set_track_states(j, taxa_tracks_thread[omp_get_thread_num()]);
#else
            set_track_states(j, taxa_tracks);
#endif
            post_probs[i][j]=0;
             for(k=0; k<the_tracks->get_num_possible_track_orders(); k++) {
#ifdef _OPEN_MP_VERSION_
                 set_track_states(k, taxa_tracks_thread2[omp_get_thread_num()]);
#else
                 set_track_states(k, taxa_tracks2);
#endif
                 trans_prob=1.0;
                 
                 //This should give the transprobs from state j in site i to state k in site i+1
#ifdef _OPEN_MP_VERSION_
                 for(l=0; l<curr_exchange->get_num_taxa(); l++) trans_prob *= track_transprobs[l][taxa_tracks_thread2[omp_get_thread_num()][l]][taxa_tracks_thread[omp_get_thread_num()][l]];
#else
                 for(l=0; l<curr_exchange->get_num_taxa(); l++) trans_prob *= track_transprobs[l][taxa_tracks2[l]][taxa_tracks[l]];
#endif
#ifdef _OPEN_MP_VERSION_
                sum_probs_thread[omp_get_thread_num()][k]=log(trans_prob)+left_cond_probs[i-1][k];
#else
                sum_probs[k]=log(trans_prob)+left_cond_probs[i-1][k];
#endif
            }
#ifdef _OPEN_MP_VERSION_
            sum_prob = get_ln_sum_prob_sort(sum_probs_thread[omp_get_thread_num()]);
#else
            sum_prob = get_ln_sum_prob_sort(sum_probs);
#endif
            
            post_probs[i][j] = right_cond_probs[i][j] + sum_prob;
            //if (taxa_tracks[0] == taxa_tracks[1])
            #pragma omp critical
            if ((i==1000) && (j< 6)) cout<<"Site "<<i<<" pattern: "<<j<<" (Ref "<<taxa_tracks[0]<<")  post prob: "<<post_probs[i][j]<<"  RCP= "<<right_cond_probs[i][j]<<" sums: "<<sum_prob<<endl;
            
        }
        
    }
    

    //Normalize the probs:
    for(i=0; i<the_homologs->get_num_homologs(); i++) {
        for(j=0; j<the_tracks->get_num_possible_track_orders(); j++) sum_probs[j]=post_probs[i][j];
        
        sum=get_ln_sum_prob_sort(sum_probs);
        //cout<<"Sum for site "<<i<<" is "<<sum<<endl;
        for(j=0; j<the_tracks->get_num_possible_track_orders(); j++) {
           // cout<<"\t"<<post_probs[i][j]<<": ";
            post_probs[i][j]=exp(post_probs[i][j]-sum);
           // cout<<"\t"<<post_probs[i][j];
        }
        //cout<<endl;
    }
    
}





void Ploidy_Like_model::calc_transprobs(Branch *taxa, int rate_num)
{
    int i,j;
    
    the_matrix->calc_transprobs(taxa);
#if 0
    if (taxa->get_parent() == 0) {
        cout<<"State\t";
        for(i=0; i<the_matrix->get_num_states(); i++)
            cout<<the_matrix->get_nth_state(i)->get_state_name()<<"\t";
        cout<<endl;
        for(i=0; i<the_matrix->get_num_states(); i++) {
            cout<<the_matrix->get_nth_state(i)->get_state_name()<<"\t";
            for (j=0; j<the_matrix->get_num_states(); j++)
                cout<<taxa->get_trpb(0,i,j)<<"\t";
            cout<<endl;
        }
    }
#endif
    
}


double Ploidy_Like_model::find_ut (Branch *taxa)
{
    return(taxa->expect_subs_site());
}



void Ploidy_Like_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
    int i, brn_cnt=0;
    
    brn_index=new int[curr_exchange->get_num_branches()-num_zero_brns];
    
    for(i=0; i<the_matrix->get_num_params(); i++) {
        par[i]=the_matrix->get_scaled_param(i);
        types[i]=ARBITRARY_PARAM;
    }
   
    brn_start=the_matrix->get_num_params();
    for(i=0; i<curr_exchange->get_num_branches(); i++) {
        if ((curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE) || ((*curr_tree)[i]->expect_subs_site()  != 0)) {
            //cout<<"Initializing branch "<<i<<" to "<<(*curr_tree)[i]->get_brnlen()<<endl;
            par[brn_cnt+brn_start]=(*curr_tree)[i]->get_brnlen();
            types[brn_cnt+brn_start]=BRANCH;
            
            if (curr_exchange->zero_len_brns_fixed() ==(BOOL)TRUE) {
                if ((*curr_tree)[i]->expect_subs_site()  != 0)
                    brn_index[brn_cnt]=i;
            }
            
            brn_cnt++;
        }
        
    }
}


void Ploidy_Like_model::set_arb_param(int param_num, double param_val)
{
    the_matrix->set_param(param_num, param_val);
}


void Ploidy_Like_model::setup_data(int dupl_pattern_index, int subposition_index)
{
	int i, site=0;
#ifdef _OPEN_MP_VERSION_
	site=omp_get_thread_num();
#endif
	
	for(i=0; i<curr_exchange->get_num_taxa(); i++)
        (*curr_data)[i].Assign_site(site, the_matrix->get_nth_level_ith_state(taxa_dupl_pattern_states[dupl_pattern_index][i], pattern_substates[dupl_pattern_index][i][subposition_index])->get_state_id());
#if 0
    cout<<"Dupl l: "<<dupl_pattern_index<<" Subpos: "<<subposition_index<<endl;
    for(i=0; i<curr_exchange->get_num_taxa(); i++) {
        cout<<(*the_genomes)[i].get_name()<<": "<<the_matrix->get_nth_state((*curr_data)[i][site])->get_state_name()<<"\t";
    }
    cout<<endl;
#endif
}


        
void Ploidy_Like_model::set_site_base_states(int site_num)
{
	int i, j, new_mask;


	for(i=0; i<curr_exchange->get_num_taxa(); i++) {
        new_mask=0;
    
        for(j=0; j<the_homologs->get_dupl_level(); j++) {
            if ((*the_tracks)[i].get_gene_track(site_num, j)->my_locus != 0) new_mask += pow2[j];
        }
        site_states[i]=the_matrix->get_masked_state(new_mask);
        //cout<<"Site "<<site_num<<" taxa "<<i<<" setting base state to "<<site_states[i]->get_state_name()<<" from mask "<<new_mask<<endl;
        
    }
	
}

void Ploidy_Like_model::set_track_states(int net_track_index, int *taxa_vals)
{
    int j,  divisor_so_far;
    
    divisor_so_far=1;
    //cout<<"Tracking "<<net_track_index;
    for (j=0; j<curr_exchange->get_num_taxa(); j++) {
        taxa_vals[j] = (net_track_index/divisor_so_far)%the_tracks->get_single_genome_num_trackings();
       // cout<<" Taxa "<<j<<" = "<<taxa_vals[j]<<"\t";
        divisor_so_far*=the_tracks->get_single_genome_num_trackings();
    }
    //cout<<endl;
}

int Ploidy_Like_model::get_tracked_pattern_subindex_from_state(int dupl_index, int net_track_index, int *taxa_vals)
{
    int j, multiplier_so_far, ret_index, my_dupl;
    
    set_track_states(net_track_index,taxa_vals);
    
    multiplier_so_far=1;
    ret_index=0;
    
   // cout<<"DuplMask: "<<dupl_index<<" Track: "<<net_track_index;
    
    for (j=0; j<curr_exchange->get_num_taxa(); j++) {
        my_dupl=((dupl_index/powN[j])%the_homologs->get_dupl_level())+1;
     // cout<<"\t"<<j<<" dupl: "<<my_dupl<<" Crossref: "<<site_states[j]->get_cross_ref(taxa_vals[j])->get_state_name();
        ret_index += (site_states[j]->get_cross_ref(taxa_vals[j])->get_level_state_id())*multiplier_so_far;

        multiplier_so_far*=num_level_patterns[my_dupl];
        //multiplier_so_far*=taxa_dupl_pattern_states[my_dupl][j];
    }
   // cout<<" Index: "<<ret_index<<endl;
    return(ret_index);
}


void Ploidy_Like_model::build_transition_vector(int locus)
{
    int i, j, k, taxon, num_diff, num_break, num_join, num_no_break, num_swap, num_break_swap;
    double  prob , track_prob, lost_prob, probx;
    
    
    if (the_homologs->get_dupl_level()== 3) {
        for(taxon=0; taxon<curr_exchange->get_num_taxa(); taxon++) {
            num_join=0;
            for(k=0; k<the_homologs->get_dupl_level(); k++) {
                if ((*the_tracks)[taxon].has_back_link(locus, k) == TRUE)
                    num_join++;
            }
            
            switch (num_join) {
                case 0:
                    probx = 1.0/((double)the_tracks->get_single_genome_num_trackings());
                    break;
                case 1:
                    probx = 1.0/(2.0 +4.0*the_matrix->get_param(switch_prob_param_num));
                    break;
                case 2:
                    probx = 1.0/(1.0+2.0*the_matrix->get_param(switch_prob_param_num)+3.0*(the_matrix->get_param(switch_prob_param_num)*the_matrix->get_param(switch_prob_param_num)));
                    break;
                case 3:
                    probx = 1.0/(1.0+3.0*(the_matrix->get_param(switch_prob_param_num)*the_matrix->get_param(switch_prob_param_num)) +
                                 2.0*(the_matrix->get_param(switch_prob_param_num)*the_matrix->get_param(switch_prob_param_num)*the_matrix->get_param(switch_prob_param_num)));
                    break;
                
            }
            
             for(i=0; i<the_tracks->get_single_genome_num_trackings(); i++) {
                 for(j=0; j<the_tracks->get_single_genome_num_trackings(); j++) {
                     num_swap=0;
                     num_break_swap=0;
                     for(k=0; k<the_homologs->get_dupl_level(); k++) {
                         if ((tracking_permutes[i][k] != tracking_permutes[j][k]) &&
                             ((*the_tracks)[taxon].has_back_link(locus, tracking_permutes[j][k]) == TRUE)) num_break_swap++;
                             
                     }
                    
                     track_transprobs[taxon][i][j]=probx;
                     for (k=0; k<num_break_swap; k++) track_transprobs[taxon][i][j]=track_transprobs[taxon][i][j]*the_matrix->get_param(switch_prob_param_num);
                
                 }
             }
        }
    }
        
 if (the_homologs->get_dupl_level()== 4) {
     for(taxon=0; taxon<curr_exchange->get_num_taxa(); taxon++) {
         num_join=0;
         for(k=0; k<the_homologs->get_dupl_level(); k++) {
             if ((*the_tracks)[taxon].has_back_link(locus, k) == TRUE)
                 num_join++;
         }
         
         switch (num_join) {
             case 0:
                 probx = 1.0/((double)the_tracks->get_single_genome_num_trackings());
                 break;
             case 1:
                 probx = 1.0/(6.0 +18.0*the_matrix->get_param(switch_prob_param_num));
                 break;
             case 2:
                 probx = 1.0/(2.0 +8.0*the_matrix->get_param(switch_prob_param_num)+
                              14.0*(the_matrix->get_param(switch_prob_param_num)*the_matrix->get_param(switch_prob_param_num)));
                 break;
             case 3:
                 probx = 1.0/(1.0 + 3.0*the_matrix->get_param(switch_prob_param_num)+
                              9.0*(the_matrix->get_param(switch_prob_param_num)*the_matrix->get_param(switch_prob_param_num))+
                              11.0*(the_matrix->get_param(switch_prob_param_num)*the_matrix->get_param(switch_prob_param_num)*the_matrix->get_param(switch_prob_param_num)));
                 break;
             case 4:
                 probx = 1.0/(1.0 + 6.0*(the_matrix->get_param(switch_prob_param_num)*the_matrix->get_param(switch_prob_param_num))+8.0*(the_matrix->get_param(switch_prob_param_num)*the_matrix->get_param(switch_prob_param_num)*the_matrix->get_param(switch_prob_param_num))+
                              9.0*(the_matrix->get_param(switch_prob_param_num)*the_matrix->get_param(switch_prob_param_num)*the_matrix->get_param(switch_prob_param_num)*the_matrix->get_param(switch_prob_param_num)));
                 break;
                 
         }
         for(i=0; i<the_tracks->get_single_genome_num_trackings(); i++) {
             for(j=0; j<the_tracks->get_single_genome_num_trackings(); j++) {
                 num_swap=0;
                 num_break_swap=0;
                 for(k=0; k<the_homologs->get_dupl_level(); k++) {
                     if ((tracking_permutes[i][k] != tracking_permutes[j][k]) &&
                         ((*the_tracks)[taxon].has_back_link(locus, tracking_permutes[j][k]) == TRUE)) num_break_swap++;
                     
                 }
                 
                 track_transprobs[taxon][i][j]=probx;
                 for (k=0; k<num_break_swap; k++) track_transprobs[taxon][i][j]=track_transprobs[taxon][i][j]*the_matrix->get_param(switch_prob_param_num);
                 
             }
         }
     }
     
#if 0
        track_prob = 1.0/((double)(the_homologs->get_dupl_level()-1));
       
        
        for(taxon=0; taxon<curr_exchange->get_num_taxa(); taxon++) {
            for(i=0; i<the_tracks->get_single_genome_num_trackings(); i++) {
                num_no_break=0;
                //Compute how to distribute the probabilty lost to track breaks among the non-breaking trackings
                for(j=0; j<the_tracks->get_single_genome_num_trackings(); j++) {
                    num_break=0;
                    for(k=0; k<the_homologs->get_dupl_level(); k++) {
                        if ((tracking_permutes[i][k] != tracking_permutes[j][k])  && ((*the_tracks)[taxon].has_back_link(locus, tracking_permutes[j][k])))
                            num_break++;
                    }
                    if (num_break == 0) num_no_break++;
                }

                lost_prob=0;
                for(j=0; j<the_tracks->get_single_genome_num_trackings(); j++) {
                    prob=1.0/((double)the_tracks->get_single_genome_num_trackings());
                    num_break=0;
                    for(k=0; k<the_homologs->get_dupl_level(); k++) {
                        if ((tracking_permutes[i][k] != tracking_permutes[j][k])  && ((*the_tracks)[taxon].has_back_link(locus, tracking_permutes[j][k])))

                            num_break++;
                    }
                    if (num_break != 0) {
                        track_transprobs[taxon][i][j]=prob*pow(the_matrix->get_param(switch_prob_param_num), num_break);
                        lost_prob += prob * (1.0-pow(the_matrix->get_param(switch_prob_param_num), num_break));
                    }
                }
                
                for(j=0; j<the_tracks->get_single_genome_num_trackings(); j++) {
                    prob=1.0/((double)the_tracks->get_single_genome_num_trackings());
                    num_break=0;
                    for(k=0; k<the_homologs->get_dupl_level(); k++) {
                        if ((tracking_permutes[i][k] != tracking_permutes[j][k])  && ((*the_tracks)[taxon].has_back_link(locus, tracking_permutes[j][k])))

                            num_break++;
                    }
                    if (num_break == 0) {
                        track_transprobs[taxon][i][j]=prob + lost_prob/((double)num_no_break);
                    }
                    
                  // cout<<"Arb size: "<<(*the_genomes)[taxon].get_name()<<": "<<locus<<", "<<i<<"->"<<j<<" breaks: "<<num_break<<"= "<<track_transprobs[taxon][i][j]<<endl;
                }
                prob=0;
                for(j=0; j<the_tracks->get_single_genome_num_trackings(); j++) prob +=track_transprobs[taxon][i][j];
                if ((prob <0.9999) || (prob >1.0001)) cerr<<"ERROR: for taxa "<<taxon<<" and tracking "<<i<<" transitions probs do not sum to 1.0: "<<prob<<endl;
            }
        }
#endif
    }
   if (the_homologs->get_dupl_level()== 2) {
        for(taxon=0; taxon<curr_exchange->get_num_taxa(); taxon++) {
            for(i=0; i<the_tracks->get_single_genome_num_trackings(); i++) {
                for(j=0; j<the_tracks->get_single_genome_num_trackings(); j++) {
                    num_join=0;
                    num_swap=0;
                    for(k=0; k<the_homologs->get_dupl_level(); k++) {
                        if (tracking_permutes[i][k] != tracking_permutes[j][k]) num_swap++;
                        
                        if ((*the_tracks)[taxon].has_back_link(locus, tracking_permutes[j][k]) == TRUE)
                            num_join++;
                    }
                    
                    if (num_swap == 0) {
                        if (num_join > 0)
                            track_transprobs[taxon][i][j]=1.0-the_matrix->get_param(switch_prob_param_num);
                        else
                            track_transprobs[taxon][i][j]=0.5;

                    }
                    else {
                        if (num_join > 0)
                            track_transprobs[taxon][i][j]=the_matrix->get_param(switch_prob_param_num);
                        else
                            track_transprobs[taxon][i][j]=0.5;
                    }
                   // cout<<"Size 2: "<<(*the_genomes)[taxon].get_name()<<": "<<locus<<", "<<i<<"->"<<j<<" breaks: "<<num_break<<"= "<<track_transprobs[taxon][i][j]<<endl;
                    
                }
            }
            
        }
    }
}

void Ploidy_Like_model::build_forward_transition_vector(int locus)
{
    int i, j, k, taxon, num_diff, num_break, num_join, num_swap, num_no_break,num_break_swap;
    double  prob , track_prob, lost_prob,  probx;
    
    if (the_homologs->get_dupl_level()== 3) {
        for(taxon=0; taxon<curr_exchange->get_num_taxa(); taxon++) {
            num_join=0;
            for(k=0; k<the_homologs->get_dupl_level(); k++) {
                if ((*the_tracks)[taxon].has_back_link(locus+1, k) == TRUE)
                    num_join++;
            }
            
            switch (num_join) {
                case 0:
                    probx = 1.0/((double)the_tracks->get_single_genome_num_trackings());
                    break;
                case 1:
                    probx = 1.0/(2.0 +4.0*the_matrix->get_param(switch_prob_param_num));
                    break;
                case 2:
                    probx = 1.0/(1.0+2.0*the_matrix->get_param(switch_prob_param_num)+3.0*(the_matrix->get_param(switch_prob_param_num)*the_matrix->get_param(switch_prob_param_num)));
                    break;
                case 3:
                    probx = 1.0/(1.0+3.0*(the_matrix->get_param(switch_prob_param_num)*the_matrix->get_param(switch_prob_param_num)) +
                                 2.0*(the_matrix->get_param(switch_prob_param_num)*the_matrix->get_param(switch_prob_param_num)*the_matrix->get_param(switch_prob_param_num)));
                    break;
                    
            }
            
            for(i=0; i<the_tracks->get_single_genome_num_trackings(); i++) {
                for(j=0; j<the_tracks->get_single_genome_num_trackings(); j++) {
                    num_swap=0;
                    num_break_swap=0;
                    for(k=0; k<the_homologs->get_dupl_level(); k++) {
                        if ((tracking_permutes[i][k] != tracking_permutes[j][k]) &&
                            ((*the_tracks)[taxon].has_back_link(locus+1, tracking_permutes[j][k]) == TRUE)) num_break_swap++;
                        
                    }
                    
                    track_transprobs[taxon][i][j]=probx;
                    for (k=0; k<num_break_swap; k++) track_transprobs[taxon][i][j]=track_transprobs[taxon][i][j]*the_matrix->get_param(switch_prob_param_num);
                    
                }
            }
        }
    }
    if (the_homologs->get_dupl_level()== 4) {
        for(taxon=0; taxon<curr_exchange->get_num_taxa(); taxon++) {
            num_join=0;
            for(k=0; k<the_homologs->get_dupl_level(); k++) {
                if ((*the_tracks)[taxon].has_back_link(locus+1, k) == TRUE)
                    num_join++;
            }
            
            switch (num_join) {
                case 0:
                    probx = 1.0/((double)the_tracks->get_single_genome_num_trackings());
                    break;
                case 1:
                    probx = 1.0/(6.0 +18.0*the_matrix->get_param(switch_prob_param_num));
                    break;
                case 2:
                    probx = 1.0/(2.0 +8.0*the_matrix->get_param(switch_prob_param_num)+
                                 14.0*(the_matrix->get_param(switch_prob_param_num)*the_matrix->get_param(switch_prob_param_num)));
                    break;
                case 3:
                    probx = 1.0/(1.0 + 3.0*the_matrix->get_param(switch_prob_param_num)+
                                 9.0*(the_matrix->get_param(switch_prob_param_num)*the_matrix->get_param(switch_prob_param_num))+
                                 11.0*(the_matrix->get_param(switch_prob_param_num)*the_matrix->get_param(switch_prob_param_num)*the_matrix->get_param(switch_prob_param_num)));
                    break;
                case 4:
                    probx = 1.0/(1.0 + 6.0*(the_matrix->get_param(switch_prob_param_num)*the_matrix->get_param(switch_prob_param_num))+8.0*(the_matrix->get_param(switch_prob_param_num)*the_matrix->get_param(switch_prob_param_num)*the_matrix->get_param(switch_prob_param_num))+
                                 9.0*(the_matrix->get_param(switch_prob_param_num)*the_matrix->get_param(switch_prob_param_num)*the_matrix->get_param(switch_prob_param_num)*the_matrix->get_param(switch_prob_param_num)));
                    break;
                    
            }
            for(i=0; i<the_tracks->get_single_genome_num_trackings(); i++) {
                for(j=0; j<the_tracks->get_single_genome_num_trackings(); j++) {
                    num_swap=0;
                    num_break_swap=0;
                    for(k=0; k<the_homologs->get_dupl_level(); k++) {
                        if ((tracking_permutes[i][k] != tracking_permutes[j][k]) &&
                            ((*the_tracks)[taxon].has_back_link(locus+1, tracking_permutes[j][k]) == TRUE)) num_break_swap++;
                        
                    }
                    
                    track_transprobs[taxon][i][j]=probx;
                    for (k=0; k<num_break_swap; k++) track_transprobs[taxon][i][j]=track_transprobs[taxon][i][j]*the_matrix->get_param(switch_prob_param_num);
                    
                }
            }
        }
    }

    
    if (the_homologs->get_dupl_level()== 2) {
        for(taxon=0; taxon<curr_exchange->get_num_taxa(); taxon++) {
            for(i=0; i<the_tracks->get_single_genome_num_trackings(); i++) {
                for(j=0; j<the_tracks->get_single_genome_num_trackings(); j++) {
                    num_join=0;
                    num_swap=0;
                    for(k=0; k<the_homologs->get_dupl_level(); k++) {
                        if (tracking_permutes[i][k] != tracking_permutes[j][k]) num_swap++;
                        
                        if ((*the_tracks)[taxon].has_back_link(locus+1, tracking_permutes[j][k]) == TRUE)
                            num_join++;
                    }
                    
                    if (num_swap == 0) {
                        if (num_join > 0)
                            track_transprobs[taxon][i][j]=1.0-the_matrix->get_param(switch_prob_param_num);
                        else
                            track_transprobs[taxon][i][j]=0.5;
                        
                    }
                    else {
                        if (num_join > 0)
                            track_transprobs[taxon][i][j]=the_matrix->get_param(switch_prob_param_num);
                        else
                            track_transprobs[taxon][i][j]=0.5;
                    }
                }
            }
            
        }
    }
}


double Ploidy_Like_model::get_ln_sum_prob_sort(double *vals)
{
	return(get_ln_sum_prob_sort(vals, the_tracks->get_num_possible_track_orders()));
}



double Ploidy_Like_model::get_ln_sum_prob_sort(double *vals, int num)
{
	int i, n;
	double  cumulative;
	
	qsort(vals, num, sizeof(double), cmpWGX);
    cumulative=vals[0];

	for(i=1; i<num; i++) cumulative=logadd(cumulative, vals[i]);

	
	return(cumulative);
}



double Ploidy_Like_model::find_smallest(double *vals, int num)
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



double Ploidy_Like_model::get_ln_sum_prob_zero_sort(double *vals, int num)
{
	int i;
	double cumulative;
	
	qsort(vals, num, sizeof(double), cmpWGX);
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
	



double Ploidy_Like_model::find_smallest_zero(double *vals)
{
	int i, start, index;
	double retval;
	
	start=0;
	while(sum_probs_used[start] ==(BOOL)TRUE) start++;
	
	index=start;
	retval=vals[start];
	for (i=start+1; i<the_tracks->get_num_possible_track_orders(); i++)
	{
		if ((sum_probs_used[i] ==(BOOL)FALSE) && ((vals[i] < retval) || (vals[i] == 1.0))) {
			index=i;
			retval=vals[i];
		}
	}

	sum_probs_used[index]=(BOOL)TRUE;
	return(retval);
}


double Ploidy_Like_model::recurse_shared_tracking(int best_track, int locus, int curr_depth, int total_depth, BOOL *&omitted_taxa, BOOL *&incoming_omit)
{
    
    //Note that the degeneracy computation assumes a genome duplication
    int taxa_id, taxa_id2, track_id;
    double new_val, ret_val=0;
    BOOL *new_omits, *local_best, match;

    //cout<<"Entering recurse_shared_tracking: "<<curr_depth<<" looking for "<<total_depth<<endl;
    new_omits=new BOOL [curr_exchange->get_num_taxa()];
    omitted_taxa=new BOOL [curr_exchange->get_num_taxa()];
    
    
    
    if (curr_depth == 0) {
        for(taxa_id=0; taxa_id<curr_exchange->get_num_taxa(); taxa_id++)
            new_omits[taxa_id]=FALSE;
    }
    else {
        for(taxa_id=0; taxa_id<curr_exchange->get_num_taxa(); taxa_id++)
            new_omits[taxa_id]=incoming_omit[taxa_id];
    }
    
    if ((curr_depth+1) < total_depth) {
        for(taxa_id=0; taxa_id<curr_exchange->get_num_taxa(); taxa_id++) {
            if (new_omits[taxa_id]==FALSE) {
                new_omits[taxa_id]=TRUE;
                new_val=recurse_shared_tracking(best_track, locus, curr_depth+1, total_depth, local_best, new_omits);
                
                if (new_val > ret_val) {
                    ret_val=new_val;
                    for(taxa_id2=0; taxa_id2<curr_exchange->get_num_taxa(); taxa_id2++)
                        omitted_taxa[taxa_id2]=local_best[taxa_id2];
                }
                delete[] local_best;
                new_omits[taxa_id]=FALSE;
            }
        }
        delete[] new_omits;
        return(ret_val);
    }
    else {
        set_track_states(best_track, taxa_tracks);
        
        
        for(taxa_id=0; taxa_id<curr_exchange->get_num_taxa(); taxa_id++) {
            if (new_omits[taxa_id]==FALSE) {
                new_omits[taxa_id]=TRUE;
                
                new_val=0;
                
                for(track_id=0; track_id<the_tracks->get_num_possible_track_orders(); track_id++) {
                    set_track_states(track_id, taxa_tracks2);
                    
                    match = TRUE;
                    
                    for(taxa_id2=0; taxa_id2<curr_exchange->get_num_taxa(); taxa_id2++) {
                        if ((new_omits[taxa_id2] == FALSE) && (taxa_tracks[taxa_id2] != taxa_tracks2[taxa_id2]))
                            match=FALSE;
                    }
                    
                    if (match == TRUE) new_val+=get_post_prob(locus, track_id);
                    
                    //if (match == TRUE) new_val+=get_post_prob(locus, track_id);
                }
                
                //cout<<locus<<" omits are ";
                //for(taxa_id2=0; taxa_id2<curr_exchange->get_num_taxa(); taxa_id2++) {
                //    if (new_omits[taxa_id2] == FALSE) cout<<"O";
                 //   else cout<<"X";
               // }
                //cout<<" and prob: "<<new_val<<endl;
                
                if (new_val > ret_val) {
                    ret_val=new_val;
                    for(taxa_id2=0; taxa_id2<curr_exchange->get_num_taxa(); taxa_id2++)
                        omitted_taxa[taxa_id2]=new_omits[taxa_id2];
                }
             new_omits[taxa_id]=FALSE;
            }
        }
      //  cout<<"For depth "<<total_depth<<" best omit prob "<<ret_val<<" ";
      //  for(taxa_id2=0; taxa_id2<curr_exchange->get_num_taxa(); taxa_id2++) {
      //      if (omitted_taxa[taxa_id2] == FALSE) cout<<"O";
      //      else cout<<"X";
      //  }
      //  cout<<endl;
        delete[] new_omits;
        return(ret_val);
    }
}

Ploidy_Like_model::~Ploidy_Like_model()
{
	int i, j, k, l;

	for(i=0; i<num_dupl_patterns; i++) {
        delete[] pattern_probs[i];
        delete[] taxa_dupl_pattern_states[i];
        for(j=0; j<curr_exchange->get_num_taxa(); j++) delete[] pattern_substates[i][j];
        delete[] pattern_substates[i];
	}
    delete[] pattern_probs;
    delete[] pattern_substates;
    delete[] taxa_dupl_pattern_states;

    delete[] taxa_tracks;
    delete[] taxa_tracks2;

    
    delete[] state_trans_probs;
	delete[] cumulative_probs;
	delete[] last_cumulative_probs;
	delete[] sum_probs;
    
    if (switch_permutes !=0) delete[] switch_permutes;
    
    if(tracking_state_lookups !=0) {
        for(i=0; i<the_homologs->get_dupl_level(); i++) {
            for (j=0; j<num_level_patterns[i]; j++)
                delete [] tracking_state_lookups[i][j];
            delete[] tracking_state_lookups[i];
        }
        delete[] tracking_state_lookups;
    }
    
    for(j=0; j<curr_exchange->get_num_taxa(); j++) {
        for(i=0; i<num_track_patterns; i++)
            delete[] track_transprobs[j][i];
        delete[] track_transprobs[j];
    }
    
    if (has_degen!=0) delete[] has_degen;
    if (gene_pos !=0) delete[] gene_pos;
    
    delete[] track_transprobs;
    delete[] site_states;
    
#ifdef _OPEN_MP_VERSION_
    for (i=0; i<curr_exchange->get_num_open_mp_threads(); i++) {
        delete[] taxa_tracks_thread[i];
        delete[] taxa_tracks_thread2[i];
        delete[] sum_probs_thread[i];
    }
	delete[] sum_probs_thread;
    delete[] taxa_tracks_thread;
    delete[] taxa_tracks_thread2;
#endif
	
	
	
	delete[] sum_probs_used;
	
#ifdef _OPEN_MP_VERSION_
	for (i=0; i<curr_exchange->get_num_open_mp_threads(); i++) delete[] sum_probs_used_thread[i];
	delete[] sum_probs_used_thread;
#endif
	

    //delete the_tracks;
    delete curr_data;
    delete[] total_pattern_states;
    delete[] num_level_patterns;
    
}

    
Ploidy_Like_RootDelta_model::Ploidy_Like_RootDelta_model() : Ploidy_Like_model()
{
    switch_permutes=0;
    tracking_state_lookups=0;
    intermediate_cond_probs=0;
    intermediate_cond_probs_thread=0;
    
    std::cerr<<"Error: call to default constructor of Ploidy_Like_RootDelta_model\n";
}


Ploidy_Like_RootDelta_model::Ploidy_Like_RootDelta_model(Exchange *cexchange, Tree *ctree, Clade *cgenomes, WGX_Data *chomologs, Phylo_Matrix *cmatrix, Phylo_Matrix *rmatrix)
{
    int i;
    
    cout<<"Constructing RootDelta model\n";
    switch_permutes=0;
    tracking_state_lookups=0;
    the_matrix=cmatrix;
    root_matrix=rmatrix;
   
    intermediate_cond_probs=0;
    intermediate_cond_probs_thread=0;
    
#ifdef _OPEN_MP_VERSION_
    intermediate_cond_probs_thread=new double*[cexchange->get_num_open_mp_threads()];
    for(i=0; i<cexchange->get_num_open_mp_threads(); i++)
        intermediate_cond_probs_thread[i]=new double[cmatrix->get_num_states()];
#else
    intermediate_cond_probs=new double[cmatrix->get_num_states()];
#endif
    
    allocate_state_model(cexchange,ctree, cgenomes, chomologs);
}

    
double Ploidy_Like_RootDelta_model::root_freq(int site)
{
    if (root_matrix->get_nth_state(site)->is_root_state() == TRUE)
        return(1.0);
    else
        return(0.0);
}

    
    
void Ploidy_Like_RootDelta_model::calc_transprobs(Branch *taxa, int rate_num)
{
    if (taxa != curr_tree->find_root())
        the_matrix->calc_transprobs(taxa);
    else {
        the_matrix->calc_transprobs(taxa);
        root_matrix->calc_transprobs(taxa);
    }
}

void Ploidy_Like_RootDelta_model::num_params_model()
{
    curr_exchange->set_num_params(the_matrix->get_num_params()+root_matrix->get_num_params()+curr_exchange->get_num_branches()-num_zero_brns);
    std::cout<<"Setting num parameters to "<<curr_exchange->get_num_params()<<std::endl;
}
    
void Ploidy_Like_RootDelta_model::intialize_parameters (double par[], PARAM_TYPE types[])
{
    int i, brn_cnt=0;
    
    for(i=0; i<the_matrix->get_num_params(); i++) {
        par[i]=the_matrix->get_scaled_param(i);
        types[i]=ARBITRARY_PARAM;
    }
    second_matrix_start=the_matrix->get_num_params();
    
    brn_index=new int[curr_exchange->get_num_branches()-num_zero_brns];
    
    for(i=0; i<root_matrix->get_num_params(); i++) {
        par[i+second_matrix_start]=root_matrix->get_scaled_param(i);
        types[i+second_matrix_start]=ARBITRARY_PARAM;
    }
    
   
    brn_start=the_matrix->get_num_params()+root_matrix->get_num_params();
    for(i=0; i<curr_exchange->get_num_branches(); i++) {
       
        if ((curr_exchange->zero_len_brns_fixed() ==(BOOL)FALSE) || ((*curr_tree)[i]->expect_subs_site()  != 0)) {
            //cout<<"Initializing branch "<<i<<" to "<<(*curr_tree)[i]->get_brnlen()<<endl;
            par[brn_cnt+brn_start]=(*curr_tree)[i]->get_brnlen();
            types[brn_cnt+brn_start]=BRANCH;
            if (curr_exchange->zero_len_brns_fixed() ==(BOOL)TRUE) {
                if ((*curr_tree)[i]->expect_subs_site()  != 0)
                    brn_index[brn_cnt]=i;
            }
            
            brn_cnt++;
        }
        
    }
    
}


void Ploidy_Like_RootDelta_model::set_arb_param(int param_num, double param_val)
{
    if (param_num<second_matrix_start)
        the_matrix->set_param(param_num, param_val);
    else
        root_matrix->set_param(param_num-second_matrix_start, param_val);
}

    
Ploidy_Like_RootDelta_model::~Ploidy_Like_RootDelta_model()
{
    if (intermediate_cond_probs !=0)
        delete[] intermediate_cond_probs;
    if (intermediate_cond_probs_thread !=0)
        delete[] intermediate_cond_probs_thread;
}
    
    
long double Ploidy_Like_RootDelta_model::prob_w_rate (int locale, int rate_num)
//This is the key likelihood calculating function.  It calls partial_prob_w_rate
//to do the post-order tree transversal and then uses the root frequencies to
//calculate the overall likelihood
{
    int i ,j, t;
    double retval=0;
    
    partial_prob_w_rate(locale, curr_tree->get_leftmost_tip(), rate_num,
                        curr_tree->find_root());
    
  
    for (i=0; i<the_matrix->get_num_states(); i++) {
#ifdef _OPEN_MP_VERSION_
        intermediate_cond_probs_thread[omp_get_thread_num()][i]=0.0;
#else
        intermediate_cond_probs[i]=0.0;
#endif
        for (j=0; j<the_matrix->get_num_states();j++)
#ifdef _OPEN_MP_VERSION_
            intermediate_cond_probs_thread[omp_get_thread_num()][i]+=curr_tree->find_root()->get_trpb(0, i, j)*curr_tree->find_root()->get_cond_prob_locale(locale, j);
#else
            intermediate_cond_probs[i]+=curr_tree->find_root()->get_trpb(0, i, j)*curr_tree->find_root()->get_cond_prob(j);
#endif
        
    }
#if 0
    cout<<"Root top Cond probs for "<<the_matrix->get_num_states()<<" states"<<endl<<flush;
    for(i=0; i<the_matrix->get_num_states(); i++) cout<<the_matrix->get_nth_state(i)->get_state_name()<<"\t";
    cout<<endl<<flush;
    for(i=0; i<the_matrix->get_num_states(); i++) cout<<curr_tree->find_root()->get_cond_prob(i)<<"\t";
    cout<<endl<<flush<<endl<<flush;
    
    cout<<"UpperRoot\t";
    for(i=0; i<the_matrix->get_num_states(); i++)
        cout<<the_matrix->get_nth_state(i)->get_state_name()<<"\t";
    cout<<endl;
    for(i=0; i<the_matrix->get_num_states(); i++) {
        cout<<the_matrix->get_nth_state(i)->get_state_name()<<"\t";
        for (j=0; j<the_matrix->get_num_states(); j++)
            cout<<curr_tree->find_root()->get_trpb(0, i, j)<<"\t";
        cout<<endl;
    }
    
    
    cout<<"Internal Cond probs\n";
    for(i=0; i<the_matrix->get_num_states(); i++) cout<<the_matrix->get_nth_state(i)->get_state_name()<<"\t";
    cout<<endl;
    for(i=0; i<the_matrix->get_num_states(); i++) cout<<intermediate_cond_probs[i]<<"\t";
    cout<<endl;
    
    
    cout<<"RootState\t";
    for(i=0; i<root_matrix->get_num_states(); i++)
        cout<<root_matrix->get_nth_state(i)->get_state_name()<<"\t";
    cout<<endl;
    for(i=0; i<root_matrix->get_num_states(); i++) {
        cout<<root_matrix->get_nth_state(i)->get_state_name()<<"\t";
        for (j=0; j<root_matrix->get_num_states(); j++)
            cout<<root_matrix->get_tp_matrix_num(curr_tree->find_root()->get_brn_num())->get_transprob(i, j, 0)<<"\t";
        cout<<endl;
    }

#endif
    for (i=0; i<curr_exchange->get_condlike_size(); i++)
        for (j=0; j<curr_exchange->get_condlike_size();j++)
#ifdef _OPEN_MP_VERSION_
            retval+=root_freq(j)*root_matrix->get_tp_matrix_num(curr_tree->find_root()->get_brn_num())->get_transprob(j, i, 0)*intermediate_cond_probs_thread[omp_get_thread_num()][i];
#else
            retval+=root_freq(j)*root_matrix->get_tp_matrix_num(curr_tree->find_root()->get_brn_num())->get_transprob(j, i, 0)*intermediate_cond_probs[i];
#endif
    return(retval);
    
} 

    

    
void Ploidy_Like_model::will_calculate_branch_transpoint_probs()
{
    int i, j, l, start_state;
    
    do_transpoint_probs=(BOOL)TRUE;
    branch_transpoint_probs=new double ****[num_dupl_patterns];
#ifdef _OPEN_MP_VERSION_
    transpoint_condprobs1_thread = new double* [curr_exchange->get_num_open_mp_threads()];
    transpoint_condprobs2_thread = new double* [curr_exchange->get_num_open_mp_threads()];
    for(i=0; i<curr_exchange->get_num_open_mp_threads(); i++) {
        transpoint_condprobs1_thread[i] = new double [the_matrix->get_num_states()];
        transpoint_condprobs2_thread[i] = new double [the_matrix->get_num_states()];
    }
#endif
    
    transpoint_condprobs1 = new double[the_matrix->get_num_states()];
    transpoint_condprobs2 = new double[the_matrix->get_num_states()];
    
    
    
    for(i=0; i<num_dupl_patterns; i++) {
        branch_transpoint_probs[i]=new double *** [total_pattern_states[i]];
        for(j=0; j<total_pattern_states[i]; j++) {
            branch_transpoint_probs[i][j]=new double ** [curr_exchange->get_num_branches()];
            for(l=0; l<curr_exchange->get_num_branches(); l++) {
                branch_transpoint_probs[i][j][l]=new double * [the_matrix->get_num_states()];
                for (start_state=0; start_state<the_matrix->get_num_states(); start_state++)
                    branch_transpoint_probs[i][j][l][start_state]=new double [the_matrix->get_num_states()];
                
            }
        }
    }
    
}
    
    

    
    
double Ploidy_Like_model::calculate_transpoint_prob(Branch *taxa, int start_state, int end_state, int locale)
{
    //Note that this function currently does not allow for ambiguous tip duplication states
    int i, j, myid=0, redund_state;
    double retval=0, new_prob, sib_prob, start_to_end_prob=0;
    Branch *curr;
    
#ifdef _OPEN_MP_VERSION_
    myid=omp_get_thread_num();
#endif
    
    for(i=0; i<the_matrix->get_num_states(); i++)
#ifdef _OPEN_MP_VERSION_
        transpoint_condprobs1_thread[myid][i]=0;
#else
        transpoint_condprobs1[i]=0;
#endif

    
    if (taxa != curr_tree->find_root()) {
        
        //Make conditional probs for this branch with only the allowed transitions
        if (taxa->is_tip() == (BOOL)TRUE) {
            redund_state=0;
            //cout<<"For "<<taxa->get_name()<<" tip state is "<<the_matrix->get_nth_state((*curr_data)[taxa->get_taxa_id()][locale])->get_state_name()
           //<<" Redudn is "<<the_matrix->get_nth_state((*curr_data)[taxa->get_taxa_id()][locale])->num_state_redunds()<<endl;
           // cout<<"Looking at "<<the_matrix->get_nth_state(start_state)->get_state_name()<<" to "<<the_matrix->get_nth_state(end_state)->get_state_name()<<endl;
            while ((redund_state<the_matrix->get_nth_state((*curr_data)[taxa->get_taxa_id()][locale])->num_state_redunds()) &&
                   (the_matrix->get_nth_state((*curr_data)[taxa->get_taxa_id()][locale])->get_ith_redund_state(redund_state) != end_state)) redund_state++;
            
            if (redund_state<the_matrix->get_nth_state((*curr_data)[taxa->get_taxa_id()][locale])->num_state_redunds()) {
                start_to_end_prob = taxa->get_trpb(0, start_state, end_state);
            }
            else start_to_end_prob=0;
        }
        else {
#ifdef _OPEN_MP_VERSION_
            start_to_end_prob=taxa->get_cond_prob_locale(myid, end_state)*taxa->get_trpb(0, start_state, end_state);
#else
            start_to_end_prob=taxa->get_cond_prob(end_state)*taxa->get_trpb(0, start_state, end_state);
#endif
        }
        
        if (start_to_end_prob!=0) {
            if (taxa->get_sibling()->is_tip() == (BOOL)TRUE)
                start_to_end_prob *=	get_tip_prob(locale, taxa->get_sibling(), 0, start_state);
            else {
                sib_prob=0;
                for(j=0; j<the_matrix->get_num_states(); j++)
#ifdef _OPEN_MP_VERSION_
                    sib_prob+=taxa->get_sibling()->get_trpb(0, start_state, j)*taxa->get_sibling()->get_cond_prob_locale(myid, j);
#else
                    sib_prob+=taxa->get_sibling()->get_trpb(0, start_state, j)*taxa->get_sibling()->get_cond_prob(j);
#endif
                start_to_end_prob*=sib_prob;
            }
        
            //cout<<"Prob of chosen branch trans and sister is : "<<start_to_end_prob<<endl;
            //Now walk down the tree with these partial conditional probs
            curr=taxa->get_parent();
            
            if (curr_tree->is_root(curr) == (BOOL)FALSE) {

                for(i=0; i<the_matrix->get_num_states(); i++) {
                    new_prob=curr->get_trpb(0, i, start_state)*start_to_end_prob;
                    
                    if (curr->get_sibling()->is_tip() == (BOOL)TRUE)
                        sib_prob = get_tip_prob(locale, curr->get_sibling(), 0, i);
                    else {
                        sib_prob=0;
                        for(j=0; j<the_matrix->get_num_states(); j++)
#ifdef _OPEN_MP_VERSION_
                            sib_prob+=curr->get_sibling()->get_trpb(0, i, j)*curr->get_sibling()->get_cond_prob_locale(myid, j);
#else
                            sib_prob+=curr->get_sibling()->get_trpb(0, i, j)*curr->get_sibling()->get_cond_prob(j);
#endif
                    }
                    
#ifdef _OPEN_MP_VERSION_
                    transpoint_condprobs1_thread[myid][i]=new_prob*sib_prob;
#else
                    transpoint_condprobs1[i]=new_prob*sib_prob;
#endif
                
                }
                
            }
            else {
#ifdef _OPEN_MP_VERSION_
                transpoint_condprobs1_thread[myid][start_state]=start_to_end_prob;
#else
                transpoint_condprobs1[start_state]=start_to_end_prob;
#endif
            }
            
           // for(i=0; i<the_matrix->get_num_states(); i++)
           //     cout<<"State "<<the_matrix->get_nth_state(i)->get_state_name()<<" has cond prob at parent start "<<transpoint_condprobs1[i]<<endl;
            
            curr=curr->get_parent();
            
            while((curr!=0) &&(curr_tree->is_root(curr) == (BOOL)FALSE)) {
                for(i=0; i<the_matrix->get_num_states(); i++) {
                    new_prob=0;
                    
                    for(j=0; j<the_matrix->get_num_states(); j++) {
#ifdef _OPEN_MP_VERSION_
                        new_prob += transpoint_condprobs1_thread[myid][j]*curr->get_trpb(0, i, j);
#else
                        new_prob += transpoint_condprobs1[j]*curr->get_trpb(0, i, j);
#endif
                    }
                    
                    if (curr->get_sibling()->is_tip() == (BOOL)TRUE)
                        sib_prob = get_tip_prob(locale, curr->get_sibling(), 0, i);
                    else {
                        sib_prob=0;
                        for(j=0; j<the_matrix->get_num_states(); j++) {
#ifdef _OPEN_MP_VERSION_
                            sib_prob+=curr->get_sibling()->get_trpb(0, i, j)*curr->get_sibling()->get_cond_prob_locale(myid, j);
#else
                            sib_prob+=curr->get_sibling()->get_trpb(0, i, j)*curr->get_sibling()->get_cond_prob(j);
#endif
                        }
                    }
#ifdef _OPEN_MP_VERSION_
                    transpoint_condprobs2_thread[myid][i]=new_prob*sib_prob;
#else
                    transpoint_condprobs2[i]=new_prob*sib_prob;
#endif
                }
                
                for(i=0; i<the_matrix->get_num_states(); i++)
#ifdef _OPEN_MP_VERSION_
                    transpoint_condprobs1_thread[myid][i]=transpoint_condprobs2_thread[myid][i];
#else
                    transpoint_condprobs1[i]=transpoint_condprobs2[i];
#endif
                
                curr=curr->get_parent();
            }
            
            //for(i=0; i<the_matrix->get_num_states(); i++)
            //cout<<"State "<<the_matrix->get_nth_state(i)->get_state_name()<<" has cond prob at base "<<transpoint_condprobs1[i]<<endl;
        }
    }
    else {
        
        if (start_state == the_matrix->get_first_root_state()->get_state_id())
#ifdef _OPEN_MP_VERSION_
            transpoint_condprobs1_thread[myid][end_state]=curr_tree->find_root()->get_cond_prob_locale(myid, end_state);
#else
            transpoint_condprobs1[end_state]=curr_tree->find_root()->get_cond_prob(end_state);
#endif
    }
    
    
    retval=0;
    for(i=0; i<the_matrix->get_num_states(); i++) {
            
#ifdef _OPEN_MP_VERSION_
        retval += curr_tree->find_root()->get_trpb(0, the_matrix->get_first_root_state()->get_state_id(), i) *
        transpoint_condprobs1_thread[myid][i];
#else
        retval += curr_tree->find_root()->get_trpb(0, the_matrix->get_first_root_state()->get_state_id(), i) *
        transpoint_condprobs1[i];
#endif
        
    }
    
    return(retval);
}
    

    
    
void Ploidy_Like_model::calc_transpoint_branch_probs()
{
    int i, j, k, l, m, n, num_prob_indices, start_state, end_state;
#ifdef MAC_TIME
    uint64_t  pretime, pptime, walktime, posttime, preptime, prepaftertime, siteaftertime, makelnltime, makelnlaftertime, sitetotal=0;
#endif
#ifdef DO_TIME
    unsigned long pretime, pptime, walktime, posttime, preptime, prepaftertime, siteaftertime, makelnltime, makelnlaftertime, sitetotal=0;
#endif
    
    double *transpoint_state_probs, sum_prob, site_prob, full_lnL, prop_prob, trans_prob;
    
    
    branch_state_probs=new double ***[the_homologs->get_num_homologs()];
    for(i=0; i<the_homologs->get_num_homologs(); i++) {
        branch_state_probs[i]=new double **[curr_exchange->get_num_branches()];
        for(j=0; j<curr_exchange->get_num_branches(); j++) {
            branch_state_probs[i][j]=new double*[the_matrix->get_num_states()];
            for(start_state=0; start_state<the_matrix->get_num_states(); start_state++)
                branch_state_probs[i][j][start_state]=new double [the_matrix->get_num_states()];
        }
    }
    
    transpoint_state_probs = new double  [the_tracks->get_num_possible_track_orders()];
    
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
        for(start_state=0; start_state<the_matrix->get_num_states(); start_state++) {
            for(end_state=0; end_state<the_matrix->get_num_states(); end_state++) {
                if ((*curr_tree)[m]->get_trpb(0, start_state, end_state)==0.0) site_prob = 1.0;
                else {
                    //num_prob_indices=get_locus_state_index(0);
                    //pair_site_mask(0, 1);
                    set_site_base_states(0);

                    build_forward_transition_vector(0);
#ifdef _OPEN_MP_VERSION_
#pragma omp parallel for private (j, k,n, sum_prob, site_prob, trans_prob)
#endif
                    
                    for(j=0; j<the_tracks->get_num_possible_track_orders(); j++) {
#ifdef _OPEN_MP_VERSION_
                        set_track_states(j, taxa_tracks_thread[omp_get_thread_num()]);
#else
                        set_track_states(j, taxa_tracks);
#endif
                        for(k=0; k<the_tracks->get_num_possible_track_orders(); k++) {
#ifdef _OPEN_MP_VERSION_
                            set_track_states(k, taxa_tracks_thread2[omp_get_thread_num()]);
#else
                            set_track_states(k, taxa_tracks2);
#endif
                            trans_prob=1.0;
#ifdef _OPEN_MP_VERSION_
                            for(n=0; n<curr_exchange->get_num_taxa(); n++) trans_prob *= track_transprobs[n][taxa_tracks_thread2[omp_get_thread_num()][n]][taxa_tracks_thread[omp_get_thread_num()][n]];
#else
                            for(n=0; n<curr_exchange->get_num_taxa(); n++) trans_prob *= track_transprobs[n][taxa_tracks2[n]][taxa_tracks[n]];
#endif
                              //cout<<"S"<<i<<": "<<j<<"->"<<k<<": "<<trans_prob<<endl;
#ifdef _OPEN_MP_VERSION_
                            sum_probs_thread[omp_get_thread_num()][k]=log(trans_prob)+right_cond_probs[1][k];
                        }
                        sum_prob=get_ln_sum_prob_sort(sum_probs_thread[omp_get_thread_num()]);
#else
                            sum_probs[k]=log(trans_prob)+right_cond_probs[1][k];
                        }
                        sum_prob=get_ln_sum_prob_sort(sum_probs);
#endif
                    
                        
#ifdef _OPEN_MP_VERSION_
                        if (branch_transpoint_probs[site_dupl_indices[0]][get_tracked_pattern_subindex_from_state(site_dupl_indices[0], j, taxa_tracks_thread[omp_get_thread_num()])][m][start_state][end_state] > 0)
                            site_prob =
                            log(branch_transpoint_probs[site_dupl_indices[0]][get_tracked_pattern_subindex_from_state(site_dupl_indices[0], j, taxa_tracks_thread[omp_get_thread_num()])][m][start_state][end_state]);
#else 
                        if (branch_transpoint_probs[site_dupl_indices[0]][get_tracked_pattern_subindex_from_state(site_dupl_indices[0], j, taxa_tracks)][m][start_state][end_state] > 0)
                            site_prob =
                                log(branch_transpoint_probs[site_dupl_indices[0]][get_tracked_pattern_subindex_from_state(site_dupl_indices[0], j, taxa_tracks)][m][start_state][end_state]);
#endif
                    
                        else
                            site_prob = 1.0;
                    
                        if (site_prob != 1.0) {
                            transpoint_state_probs[j] = site_prob + sum_prob;
                           // cout<<"For track pattern "<<j<<" on branch "<<(*curr_tree)[m]->get_name()<<" at site 0 and start = "<<the_matrix->get_nth_state(start_state)->get_state_name()<<" and end = "<<the_matrix->get_nth_state(end_state)->get_state_name()<<" Dupl state is "<<site_dupl_indices[0]<<" prob is "<<transpoint_state_probs[j]<<endl;
                        }
                        else
                            transpoint_state_probs[j] = 1.0;
                    } //End of parallel for loop over j
                
                    site_prob=get_ln_sum_prob_zero_sort(transpoint_state_probs, the_tracks->get_num_possible_track_orders());
                }
                //cout<<"On branch "<<(*curr_tree)[m]->get_name()<<" at site 0 and start = "<<the_matrix->get_nth_state(start_state)->get_state_name()<<" and end = "<<the_matrix->get_nth_state(end_state)->get_state_name()<<" prob is "<<site_prob<<endl;
            if (site_prob != 1.0) {
                     branch_state_probs[0][m][start_state][end_state] = exp(site_prob - full_lnL);
                    cout<<"On branch "<<(*curr_tree)[m]->get_name()<<" at site 0 and start = "<<the_matrix->get_nth_state(start_state)->get_state_name()<<" and end = "<<the_matrix->get_nth_state(end_state)->get_state_name()<<" prob is "<<exp(site_prob - full_lnL)<<endl;
            }
                else
                     branch_state_probs[0][m][start_state][end_state]=0.0;
            
            }
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
            for(start_state=0; start_state<the_matrix->get_num_states(); start_state++) {
                for(end_state=0; end_state<the_matrix->get_num_states(); end_state++) {
                    if ((*curr_tree)[m]->get_trpb(0, start_state, end_state) == 0.0) site_prob = 1.0;
                    else {
                        set_site_base_states(i);
                        
                        build_transition_vector(i);
#ifdef _OPEN_MP_VERSION_
#pragma omp parallel for private (j, k, sum_prob, site_prob)
#endif
                        for(j=0; j<the_tracks->get_num_possible_track_orders(); j++) {
                            
#ifdef _OPEN_MP_VERSION_
                            set_track_states(j, taxa_tracks_thread[omp_get_thread_num()]);
#else
                            set_track_states(j, taxa_tracks);
#endif

                            
                            for(k=0; k<the_tracks->get_num_possible_track_orders(); k++) {
#ifdef _OPEN_MP_VERSION_
                                set_track_states(k, taxa_tracks_thread2[omp_get_thread_num()]);
#else
                                set_track_states(k, taxa_tracks2);
#endif
                                trans_prob=1.0;
#ifdef _OPEN_MP_VERSION_
                                for(n=0; n<curr_exchange->get_num_taxa(); n++) trans_prob *= track_transprobs[n][taxa_tracks_thread2[omp_get_thread_num()][n]][taxa_tracks_thread[omp_get_thread_num()][n]];
#else
                                for(n=0; n<curr_exchange->get_num_taxa(); n++) trans_prob *= track_transprobs[n][taxa_tracks2[n]][taxa_tracks[n]];
#endif

#ifdef _OPEN_MP_VERSION_
                                sum_probs_thread[omp_get_thread_num()][k]=log(trans_prob)+left_cond_probs[i-1][k];
#else
                                sum_probs[k]=log(trans_prob)+left_cond_probs[i-1][k];
#endif
                            }
#ifdef _OPEN_MP_VERSION_
                            sum_prob = get_ln_sum_prob_sort(sum_probs_thread[omp_get_thread_num()]);
#else
                            sum_prob = get_ln_sum_prob_sort(sum_probs);
#endif
                            
                            
#ifdef _OPEN_MP_VERSION_
                            if (branch_transpoint_probs[site_dupl_indices[i]][get_tracked_pattern_subindex_from_state(site_dupl_indices[i], j, taxa_tracks_thread[omp_get_thread_num()])][m][start_state][end_state] > 0)
                                site_prob =
                                log(branch_transpoint_probs[site_dupl_indices[i]][get_tracked_pattern_subindex_from_state(site_dupl_indices[i], j, taxa_tracks_thread[omp_get_thread_num()])][m][start_state][end_state]);
#else
                            if (branch_transpoint_probs[site_dupl_indices[i]][get_tracked_pattern_subindex_from_state(site_dupl_indices[i], j, taxa_tracks)][m][start_state][end_state] > 0)
                                site_prob =
                                log(branch_transpoint_probs[site_dupl_indices[i]][get_tracked_pattern_subindex_from_state(site_dupl_indices[i], j, taxa_tracks)][m][start_state][end_state]);
#endif
                            else
                                site_prob=1.0;

                            if (site_prob != 1.0)
                                transpoint_state_probs[j] = site_prob + sum_prob;
                            else
                                transpoint_state_probs[j] = 1.0;
                        }
                        
                        
                        //pair_site_mask(i, i+1);
                        build_forward_transition_vector(i);
                        
#ifdef _OPEN_MP_VERSION_
#pragma omp parallel for private (j, k, site_prob)
#endif
                        for(j=0; j<the_tracks->get_num_possible_track_orders(); j++) {
#ifdef _OPEN_MP_VERSION_
                            set_track_states(j, taxa_tracks_thread[omp_get_thread_num()]);
#else
                            set_track_states(j, taxa_tracks);
#endif
                            for(k=0; k<the_tracks->get_num_possible_track_orders(); k++) {
#ifdef _OPEN_MP_VERSION_
                                set_track_states(k, taxa_tracks_thread2[omp_get_thread_num()]);
#else
                                set_track_states(k, taxa_tracks2);
#endif
                                trans_prob=1.0;
#ifdef _OPEN_MP_VERSION_
                                for(n=0; n<curr_exchange->get_num_taxa(); n++) trans_prob *= track_transprobs[n][taxa_tracks_thread2[omp_get_thread_num()][n]][taxa_tracks_thread[omp_get_thread_num()][n]];
#else
                                for(n=0; n<curr_exchange->get_num_taxa(); n++) trans_prob *= track_transprobs[n][taxa_tracks2[n]][taxa_tracks[n]];
#endif

                                
#ifdef _OPEN_MP_VERSION_
                                sum_probs_thread[omp_get_thread_num()][k]=log(trans_prob)+right_cond_probs[i+1][k];
#else
                                sum_probs[k]=log(trans_prob)+right_cond_probs[i+1][k];
#endif
                            }
                            if (transpoint_state_probs[j] != 1.0)
#ifdef _OPEN_MP_VERSION_
                                transpoint_state_probs[j] += get_ln_sum_prob_sort(sum_probs_thread[omp_get_thread_num()]);
#else
                            transpoint_state_probs[j] += get_ln_sum_prob_sort(sum_probs);
#endif
                            
                        }
                        
                        site_prob=get_ln_sum_prob_zero_sort(transpoint_state_probs, the_tracks->get_num_possible_track_orders());
                    }
                    
                    if (site_prob != 1.0)
                        branch_state_probs[i][m][start_state][end_state] = exp(site_prob - full_lnL);
                    else
                        branch_state_probs[i][m][start_state][end_state]=0.0;
                    
                    //if (start_state==end_state==2) cout<<"Site "<<i<<" Branch "<<(*curr_tree)[m]->get_name()<<" State: "<<the_matrix->get_nth_state(start_state)->get_state_name()
                    //    <<" self_trans has prob "<<branch_state_probs[i][m][start_state][end_state]<<endl;
                    
                }
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
    
    for(m=0; m<curr_exchange->get_num_branches(); m++) {
        for(start_state=0; start_state<the_matrix->get_num_states(); start_state++) {
            for(end_state=0; end_state<the_matrix->get_num_states(); end_state++) {
                if ((*curr_tree)[m]->get_trpb(0, start_state, end_state) == 0) site_prob=1.0;
                else {
                    set_site_base_states(the_homologs->get_num_homologs()-1);
                    
                    build_transition_vector(the_homologs->get_num_homologs()-1);
#ifdef _OPEN_MP_VERSION_
#pragma omp parallel for private (j, k, site_prob, trans_prob)
#endif
                    for(j=0; j<the_tracks->get_num_possible_track_orders(); j++) {
#ifdef _OPEN_MP_VERSION_
                        set_track_states(j, taxa_tracks_thread[omp_get_thread_num()]);
#else
                        set_track_states(j, taxa_tracks);
#endif
                        for(k=0; k<the_tracks->get_num_possible_track_orders(); k++) {
#ifdef _OPEN_MP_VERSION_
                            set_track_states(k, taxa_tracks_thread2[omp_get_thread_num()]);
#else
                            set_track_states(k, taxa_tracks2);
#endif
                            trans_prob=1.0;
#ifdef _OPEN_MP_VERSION_
                            for(n=0; n<curr_exchange->get_num_taxa(); n++) trans_prob *= track_transprobs[n][taxa_tracks_thread2[omp_get_thread_num()][n]][taxa_tracks_thread[omp_get_thread_num()][n]];
#else
                            for(n=0; n<curr_exchange->get_num_taxa(); n++) trans_prob *= track_transprobs[n][taxa_tracks2[n]][taxa_tracks[n]];
#endif
                        
#ifdef _OPEN_MP_VERSION_
                            sum_probs_thread[omp_get_thread_num()][k]=log(trans_prob)+left_cond_probs[the_homologs->get_num_homologs()-2][k];
#else
                            sum_probs[k]=log(trans_prob)+left_cond_probs[the_homologs->get_num_homologs()-2][k];
#endif
                        }
                        
#ifdef _OPEN_MP_VERSION_
                        sum_prob = get_ln_sum_prob_sort(sum_probs_thread[omp_get_thread_num()]);
#else
                        sum_prob = get_ln_sum_prob_sort(sum_probs);
#endif
           
#ifdef _OPEN_MP_VERSION_
                        if (branch_transpoint_probs[site_dupl_indices[the_homologs->get_num_homologs()-1]][get_tracked_pattern_subindex_from_state(site_dupl_indices[the_homologs->get_num_homologs()-1], j, taxa_tracks_thread[omp_get_thread_num()])][m][start_state][end_state] > 0)
                            site_prob =
                            log(branch_transpoint_probs[site_dupl_indices[the_homologs->get_num_homologs()-1]][get_tracked_pattern_subindex_from_state(site_dupl_indices[the_homologs->get_num_homologs()-1], j, taxa_tracks_thread[omp_get_thread_num()])][m][start_state][end_state]);
#else
                        if (branch_transpoint_probs[site_dupl_indices[the_homologs->get_num_homologs()-1]][get_tracked_pattern_subindex_from_state(site_dupl_indices[the_homologs->get_num_homologs()-1], j, taxa_tracks)][m][start_state][end_state] > 0)
                            site_prob =
                            log(branch_transpoint_probs[site_dupl_indices[the_homologs->get_num_homologs()-1]][get_tracked_pattern_subindex_from_state(site_dupl_indices[the_homologs->get_num_homologs()-1], j, taxa_tracks)][m][start_state][end_state]);
#endif
                        else
                            site_prob = 1.0;
                        
                        if (site_prob != 1.0)
                            transpoint_state_probs[j] = site_prob + sum_prob; 
                        else
                            transpoint_state_probs[j] = 1.0;
                    }
                    
                    site_prob=get_ln_sum_prob_zero_sort(transpoint_state_probs, the_tracks->get_num_possible_track_orders());
                }
                if (site_prob != 1.0)
                    branch_state_probs[the_homologs->get_num_homologs()-1][m][start_state][end_state]= exp(site_prob - full_lnL);
                else
                    branch_state_probs[the_homologs->get_num_homologs()-1][m][start_state][end_state]=0.0;
            }
        }
    }
    delete[] transpoint_state_probs;
}  //End calc_transpoint_branch_probs


