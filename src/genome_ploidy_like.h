//Copyright 1999-2014 Gavin Conant

#include <math.h>
#include "lin_alg.h"
#include "gen_dna_funcs.h"
#include "maxlike.h"
#include "genome_tripl_list.h"
#include "phylo_model_matrix.h"

#ifndef ___GENOME_PLOIDY_LIKE_H___
#define ___GENOME_PLOIDY_LIKE_H___



//#define LIKE_SCALE 1e240
//#define SQRT_LIKE_SCALE 1e120

#define MAX_TAXA 20

class Ploidy_Like_model : public Like_model
{
public:
    WGX_Tracks *the_tracks;
    int **tracking_permutes, *gene_pos;
    double ****branch_state_probs;
    BOOL *has_degen;
    
    Ploidy_Like_model();
    Ploidy_Like_model(Exchange *cexchange, Tree *ctree,
                      Clade *cgenomes, WGX_Data *chomologs, Phylo_Matrix *cmatrix);
	double root_freq(int site);
	double get_basefreq(int base, Branch *tree_branch) {return(0.0);};
	double get_basefreq(int base, int codon_pos, Branch *tree_branch) {return(0.0);};
	//void get_site_state_probs(double **&prob_array, int taxa_id);
    void allocate_state_model(Exchange *cexchange, Tree *ctree, Clade *cgenomes, WGX_Data *chomologs);
    double find_appropriate_ln_like();
    void print_tracking_probs(string outfile);
    void print_transpoint_branch_probs(string filename);
    void build_transition_vector(int locus);
    double get_post_prob(int site, int pattern);
    //virtual void get_site_state_probs(int taxa_id);
    double get_site_prob(int site, Model_State *state);
    void will_calculate_branch_transpoint_probs();
    int get_num_params_model() {return (the_matrix->get_num_params()+curr_exchange->get_num_branches()-num_zero_brns);};
    int get_num_track_patterns()  {return(num_track_patterns);};
    void calc_transprobs(Branch *taxa, int rate_num);
    double find_ut (Branch *taxa);
    void num_params_model();
    void set_arb_param(int param_num, double param_val);
    void intialize_parameters (double par[], PARAM_TYPE types[]);
    void describe_results();
    void get_gene_conditional_probs();
    double get_expect_sub_site(double brlen)  {return(brlen);};
    void set_track_states(int net_track_index, int *taxa_vals);
    double recurse_shared_tracking(int best_track, int locus, int curr_depth, int total_depth, BOOL *&omitted_taxa, BOOL *&incoming_omits);
    double get_track_transprob(int genome, int from_state, int to_state) {return(track_transprobs[genome][from_state][to_state]);};
    void set_site_base_states(int site_num);
    int get_track_position_state(int tracking_id, int position_id) {return(tracking_permutes[tracking_id][position_id]);};
    Clade * get_the_genomes()   {return(the_genomes);};
    WGX_Data * get_the_homologs() {return(the_homologs);};
    
    ~Ploidy_Like_model();
    
protected:
    //Old mask setup--lets try again
    //int **masks, *dupl_states, transition_mask, *gene_counts, *site_indices, *dupl_indices;
    
    int **tracking_masks_by_level, powN[MAX_TAXA], num_dupl_patterns, **taxa_dupl_pattern_states, ***pattern_substates,
    *num_level_patterns, pow2[MAX_TAXA], *total_pattern_states, num_track_patterns, *site_dupl_indices, *taxa_tracks, *taxa_tracks2, switch_prob_param_num, num_zero_brns;
#ifdef _OPEN_MP_VERSION_
    int **taxa_tracks_thread, **taxa_tracks_thread2;
#endif
    //taxa_dupl_pattern_states gives the duplication level for each taxa at each of N^numtaxa possiblities.  Indexed over N^numtaxa first and then by #taxa
    
    //total_pattern_states is indexed over N^numtaxa.  Gives the total # of states for each of the states in the N^numtaxa
    
    //pattern_substates gives the state# for each taxa for each position in the N^numtaxa over all total_pattern_states.  Indices are N^numtaxa, then num_taxa, then total_pattern_states for that level
    
    double **pattern_probs, *state_trans_probs, *cumulative_probs, *last_cumulative_probs, ***track_transprobs,
    *sum_probs, *sum_probs_r, **post_probs, **left_cond_probs, **right_cond_probs, **site_probs, *switch_permutes, *****branch_transpoint_probs, *track_switch_probs;


    
    //double ***pattern_probs, *state_trans_probs, *cumulative_probs, *last_cumulative_probs,
    //*sum_probs, *sum_probs_r, **post_probs, **left_cond_probs, **right_cond_probs, **site_probs,
    //*****branch_transpoint_probs;
#ifdef _OPEN_MP_VERSION_
    double **sum_probs_thread;
#endif
    BOOL *sum_probs_used, *is_ambig;
#ifdef _OPEN_MP_VERSION_
    BOOL **sum_probs_used_thread;
#endif
    
    Clade *the_genomes;
    WGX_Data *the_homologs;
    Phylo_Matrix *the_matrix;
    Model_State ****tracking_state_lookups, **site_states;

    
	BOOL do_transpoint_probs;
	double *transpoint_condprobs1, *transpoint_condprobs2;
#ifdef _OPEN_MP_VERSION_
	double **transpoint_condprobs1_thread, **transpoint_condprobs2_thread;
#endif
    
    
    void initialize_arrays()  {};
	void partial_prob_w_rate(int locale, Branch *lsib, int rate_num, Branch *stop_id);
	double get_tip_prob(int locale, Branch *taxa, int rate_num, int cond_state);
    
     void calc_pattern_probs();
    void setup_data(int dupl_pattern_index, int subposition_index);
    //void build_transition_vector(int locus);
    void build_transition_vector(int locus, BOOL look_back);
    void build_forward_transition_vector(int locus);
    double get_ln_sum_prob_sort(double *vals);
    double get_ln_sum_prob_sort(double *vals, int num);
    double find_smallest(double *vals, int num);
    double get_ln_sum_prob_zero_sort(double *vals, int num);
    double find_smallest_zero(double *vals);
    
    void setup_tracking_state_lookups(Exchange *cexchange);
    int get_tracked_pattern_subindex_from_state(int dupl_index, int net_track_index, int *taxa_vals);
    void calc_transpoint_branch_probs();
    double calculate_transpoint_prob(Branch *taxa, int start_state, int end_state, int locale);
    
};


class Ploidy_Like_RootDelta_model : public Ploidy_Like_model
{
public:
    
    
    Ploidy_Like_RootDelta_model();
    Ploidy_Like_RootDelta_model(Exchange *cexchange, Tree *ctree,
                      Clade *cgenomes, WGX_Data *chomologs, Phylo_Matrix *cmatrix, Phylo_Matrix *rmatrix);
    double root_freq(int site);
    void calc_transprobs(Branch *taxa, int rate_num);
    void num_params_model();
    void set_arb_param(int param_num, double param_val);
    void intialize_parameters (double par[], PARAM_TYPE types[]);
    ~Ploidy_Like_RootDelta_model();
    
protected:
    int second_matrix_start;
    double *intermediate_cond_probs, **intermediate_cond_probs_thread;
    Phylo_Matrix *root_matrix;

    long double prob_w_rate (int locale, int rate_num);
    
};
#endif
