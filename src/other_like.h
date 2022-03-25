//Copyright 1999-2005 Gavin Conant

#include <math.h>
#include "gen_code.h"
#include "lin_alg.h"
#include "gen_dna_funcs.h"
#include "maxlike.h"
#include "score_matrix.h"
#include "genome_list.h"

#ifndef ___OTHER_LIKE_H___
#define ___OTHER_LIKE_H___

#define NUM_TRANSPOINT_STATES 3

enum TRANSPOINT_STATES {DUPLICATE, DUPLICATE_TO_SINGLE, SINGLE};

TRANSPOINT_STATES int_to_transpoint_state(int n);
int transpoint_state_to_int(TRANSPOINT_STATES state);

BOOL is_trans_in_transpoint_state(DUPL_LOSS_STATES starting_state, DUPL_LOSS_STATES ending_state, 
								  TRANSPOINT_STATES the_state);


//#define LIKE_SCALE 1e240
//#define SQRT_LIKE_SCALE 1e120

class Dupl_Base_model : virtual public Like_model
{
public:
	double root_freq(int site);	
	double get_basefreq(int base, Branch *tree_branch) {return(0.0);};
	double get_basefreq(int base, int codon_pos, Branch *tree_branch) {return(0.0);};
	virtual double get_site_prob(int site, DUPL_LOSS_STATES state) {return(0.0);};
	virtual void get_site_state_probs(double **&prob_array, int taxa_id);

protected:
	BOOL do_transpoint_probs;
	double *transpoint_condprobs1, *transpoint_condprobs2;
#ifdef _OPEN_MP_VERSION_
	double **transpoint_condprobs1_thread, **transpoint_condprobs2_thread;
#endif
	void initialize_arrays();
	void partial_prob_w_rate(int locale, Branch *lsib, int rate_num, Branch *stop_id);
	void partial_prob_w_rate_nonhidden(int locale, Branch *lsib, int rate_num, Branch *stop_id, int nonhidden_taxa);
	double get_tip_prob(int locale, Branch *taxa, int rate_num, int cond_state);
	virtual int state_redundancy_val(DUPL_LOSS_STATES the_state)=0;
	virtual DUPL_LOSS_STATES get_redund_position_n(DUPL_LOSS_STATES the_state, int value)=0;
	virtual BOOL allowed_state(DUPL_LOSS_STATES the_state);
	double calculate_transpoint_prob(Branch *taxa, TRANSPOINT_STATES state, int locale, int rate_num);
};

class Dupl_Null_Base_model : virtual public Dupl_Base_model
{
public:
	double get_expect_sub_site(double brlen)  {return(2.0*brlen);};
protected:
	void calc_transprobs(Branch *taxa, int rate_num);
	int state_redundancy_val(DUPL_LOSS_STATES the_state);
	DUPL_LOSS_STATES get_redund_position_n(DUPL_LOSS_STATES the_state, int value);
	double find_ut(Branch *taxa);
	

};

class Dupl_State_Base_model : virtual public Dupl_Base_model
{
protected:
	double **branch_transpoint_probs;
};


class Dupl_NoState_Base_model : virtual public Dupl_Base_model
{
public:
	WGD_Tracks *the_tracks;
	double ***branch_state_probs;

	void allocate_state_model(Exchange *cexchange, Tree *ctree, Clade *cgenomes, WGD_Data *chomologs);
	double find_appropriate_ln_like();
	void print_tracking_probs(char *outfile);
	double get_post_prob(int site, int pattern);
	virtual void get_site_state_probs(int taxa_id);
	double get_site_prob(int site, DUPL_LOSS_STATES state);
	void will_calculate_branch_transpoint_probs();
	void print_transpoint_branch_probs(char * filename);
	void calc_transpoint_branch_probs();
	
	~Dupl_NoState_Base_model();
protected:
	int pow2[64], **masks, *dupl_states, transition_mask, *gene_counts, *site_indices, *dupl_indices; 
	double ***pattern_probs, *state_trans_probs, *cumulative_probs, *last_cumulative_probs, 
		*sum_probs, *sum_probs_r, **post_probs, **left_cond_probs, **right_cond_probs, **site_probs,
		*****branch_transpoint_probs;
#ifdef _OPEN_MP_VERSION_
	double **sum_probs_thread;
#endif
	BOOL *sum_probs_used, *is_ambig;
#ifdef _OPEN_MP_VERSION_
	BOOL **sum_probs_used_thread;
#endif
	
	Clade *the_genomes;
	WGD_Data *the_homologs;	
	
	void get_gene_conditional_probs();
	
	int get_locus_num_dupls(int locus);
	int get_locus_state_index(int locus);
	int get_dupl_state_index();
	void recurse_states(int num_dupls, int &so_far, int pos);
	void get_states(int index_num, int num_dupls, int pos);
	void build_masks();
	void site_mask(int site_num);
	void pair_site_mask(int this_site, int last_site);
	void calc_pattern_probs();
	void setup_data(int num_dupls, int *dupl_positions, int track_index);
	void build_transition_vector(int locus); 
	void build_transition_vector(int locus, BOOL look_back); 
	double get_ln_sum_prob(double *vals);
	double get_ln_sum_prob_sort(double *vals);
	double get_ln_sum_prob(double *vals, int num);
	double get_ln_sum_prob_sort(double *vals, int num);
	double find_smallest(double *vals, int num);
	double get_ln_sum_prob_zero(double *vals, int num);
	double get_ln_sum_prob_zero_sort(double *vals, int num);
	double find_smallest_zero(double *vals);
	
};


class Dupl_model : public Dupl_Null_Base_model, public Dupl_State_Base_model
{
public:
  Dupl_model () {};
  Dupl_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree)
    {assemble (cexchange, cdata, ctree);};
  
  void describe_results();  
  void num_params_model();
  
 protected:
  void intialize_parameters (double par[], PARAM_TYPE types[]);
};



class Dupl_NoState_model : public Dupl_NoState_Base_model, public Dupl_Null_Base_model
{
public:
	Dupl_NoState_model()  {cerr<<"Error: Call to default constructor of class Dupl_NoState_model\n";};
	Dupl_NoState_model (Exchange *cexchange, Tree *ctree, Clade *cgenomes, WGD_Data *chomologs)
			{allocate_state_model(cexchange, ctree, cgenomes, chomologs); left_cond_probs=right_cond_probs=0;}; 
	void describe_results() {};
	void num_params_model();
protected:
	void intialize_parameters (double par[], PARAM_TYPE types[]);
};




class Dupl_Fix_Base_model : virtual public Dupl_Base_model
{
public:
	double get_expect_sub_site(double brlen) {return((2.0+curr_exchange->get_dupl_fix_rate())*brlen);};
protected:
	void calc_transprobs(Branch *taxa, int rate_num);
	virtual int state_redundancy_val(DUPL_LOSS_STATES the_state);
	virtual DUPL_LOSS_STATES get_redund_position_n(DUPL_LOSS_STATES the_state, int value);
	virtual double find_ut (Branch *taxa);
	virtual BOOL allowed_state(DUPL_LOSS_STATES the_state);
};



class Dupl_Fix_model : public Dupl_Fix_Base_model, public Dupl_State_Base_model
{
public:
	Dupl_Fix_model () {};
	Dupl_Fix_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree);
	void describe_results();
	void num_params_model();
	
protected:
	void intialize_parameters (double par[], PARAM_TYPE types[]);
};



class Dupl_Fix_NoState_model : public Dupl_Fix_Base_model, public Dupl_NoState_Base_model
{
public:
	Dupl_Fix_NoState_model () {cerr<<"Error: Call to default constructor of class Dupl_Fix_NoState_model\n";};
	Dupl_Fix_NoState_model (Exchange *cexchange, Tree *ctree, Clade *cgenomes, WGD_Data *chomologs);
	
	void describe_results();
	void num_params_model();

protected:
	void intialize_parameters (double par[], PARAM_TYPE types[]);
};



class Dupl_Parallel_Base_model : virtual public Dupl_Base_model 
{
public:
	double get_expect_sub_site(double brlen)  {return((2.0+2.0*curr_exchange->get_dupl_parallel_rate())*brlen);};
protected:
	void calc_transprobs(Branch *taxa, int rate_num);
	virtual int state_redundancy_val(DUPL_LOSS_STATES the_state);
	virtual DUPL_LOSS_STATES get_redund_position_n(DUPL_LOSS_STATES the_state, int value);
	virtual double find_ut (Branch *taxa);
	virtual BOOL allowed_state(DUPL_LOSS_STATES the_state);
};


class Dupl_Parallel_model : public Dupl_Parallel_Base_model, public Dupl_State_Base_model
{
public:
	Dupl_Parallel_model () {};
	Dupl_Parallel_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree);
	void describe_results();
	void num_params_model();

protected:
	void intialize_parameters (double par[], PARAM_TYPE types[]);
};



class Dupl_Parallel_NoState_model : public Dupl_Parallel_Base_model, public Dupl_NoState_Base_model
{
	public:
	Dupl_Parallel_NoState_model () {};
	Dupl_Parallel_NoState_model (Exchange *cexchange, Tree *ctree, Clade *cgenomes, WGD_Data *chomologs);
		
	void describe_results();
	void num_params_model();

protected:
	void intialize_parameters (double par[], PARAM_TYPE types[]);
};


class Dupl_Parallel_2_Rate_Base_model : virtual public Dupl_Parallel_Base_model
{
protected:
	void calc_transprobs(Branch *taxa, int rate_num);
};


class Dupl_Parallel_2_Rate_model : public Dupl_Parallel_2_Rate_Base_model, public Dupl_State_Base_model
{
	public:
	Dupl_Parallel_2_Rate_model () {};
	Dupl_Parallel_2_Rate_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree);
	void describe_results();
	void num_params_model();

protected:
	void intialize_parameters (double par[], PARAM_TYPE types[]);

};


class Dupl_Parallel_2_Rate_NoState_model : public Dupl_Parallel_2_Rate_Base_model, public Dupl_NoState_Base_model
{
		public:
	Dupl_Parallel_2_Rate_NoState_model () {};
	Dupl_Parallel_2_Rate_NoState_model (Exchange *cexchange, Tree *ctree, Clade *cgenomes, WGD_Data *chomologs);
		
	void describe_results();
	void num_params_model();

protected:
	void intialize_parameters (double par[], PARAM_TYPE types[]);
};


class Dupl_Fix_Parallel_Base_model : virtual public Dupl_Base_model
{
public:
	double get_expect_sub_site(double brlen)    {return((2.0+2.0*curr_exchange->get_dupl_parallel_rate()+
		curr_exchange->get_dupl_fix_rate())*brlen);};	
protected:
	void set_null_transprobs(Branch *taxa, int rate_num);
	void calc_transprobs(Branch *taxa, int rate_num);
	virtual int state_redundancy_val(DUPL_LOSS_STATES the_state);
	virtual DUPL_LOSS_STATES get_redund_position_n(DUPL_LOSS_STATES the_state, int value);
	virtual double find_ut (Branch *taxa);
	virtual BOOL allowed_state(DUPL_LOSS_STATES the_state);
};


class Dupl_Fix_Parallel_SubF_Base_model : virtual public Dupl_Fix_Parallel_Base_model
{

protected:
	void calc_transprobs(Branch *taxa, int rate_num);

};


class Dupl_SubF_3_Rate_Base_model : virtual public Dupl_Fix_Parallel_Base_model
{
protected:
	void calc_transprobs(Branch *taxa, int rate_num);
};


class Dupl_2_Rate_NoSubF_Base_model : virtual public Dupl_Fix_Parallel_Base_model
{
protected:
	void calc_transprobs(Branch *taxa, int rate_num);
};



class Dupl_SubF_Only_Base_model : virtual public Dupl_Fix_Parallel_Base_model
{
protected:
	void calc_transprobs(Branch *taxa, int rate_num);
};





class Dupl_Subf_All_States_Base_model: virtual public Dupl_Base_model
{
	public:
	double get_expect_sub_site(double brlen)    {return((2.0+2.0*curr_exchange->get_dupl_parallel_rate()+
		curr_exchange->get_dupl_fix_rate())*brlen);};	
protected:
	void calc_transprobs(Branch *taxa, int rate_num);
	virtual int state_redundancy_val(DUPL_LOSS_STATES the_state);
	virtual DUPL_LOSS_STATES get_redund_position_n(DUPL_LOSS_STATES the_state, int value);
	virtual double find_ut (Branch *taxa);
	virtual BOOL allowed_state(DUPL_LOSS_STATES the_state);

};


class Dupl_SubF_3_Rate_AllStates_Base_model: virtual public Dupl_Subf_All_States_Base_model
{
protected:
	void calc_transprobs(Branch *taxa, int rate_num);
};

class Dupl_2_Rate_NoSubF_AllStates_Base_model : virtual public Dupl_Subf_All_States_Base_model
{
protected:
	void calc_transprobs(Branch *taxa, int rate_num);
	virtual BOOL allowed_state(DUPL_LOSS_STATES the_state);
};


class Dupl_Slow_Loss_Con_Fix_Base_model : virtual public Dupl_Fix_Parallel_Base_model
{
protected:
	void calc_transprobs(Branch *taxa, int rate_num);
};





class Dupl_Fix_Parallel_model : public Dupl_Fix_Parallel_Base_model, public Dupl_State_Base_model
{
public:
	Dupl_Fix_Parallel_model () {};
	Dupl_Fix_Parallel_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree);
	void describe_results();
	void num_params_model();

protected:
	void intialize_parameters (double par[], PARAM_TYPE types[]);
};


class Dupl_Fix_Parallel_SubF_model : public Dupl_Fix_Parallel_SubF_Base_model, public Dupl_State_Base_model
{
public:
	Dupl_Fix_Parallel_SubF_model () {};
	Dupl_Fix_Parallel_SubF_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree);
	void describe_results();
	void num_params_model();

protected:
	void intialize_parameters (double par[], PARAM_TYPE types[]);
};


class Dupl_SubF_3_Rate_model : public Dupl_SubF_3_Rate_Base_model, public Dupl_State_Base_model
{
public:
	Dupl_SubF_3_Rate_model () {};
	Dupl_SubF_3_Rate_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree);
	void describe_results();
	void num_params_model();

protected:
	void intialize_parameters (double par[], PARAM_TYPE types[]);
};


class Dupl_SubF_3_Rate_NoState_model : public Dupl_SubF_3_Rate_Base_model, public Dupl_NoState_Base_model
{
public:
	Dupl_SubF_3_Rate_NoState_model () {cerr<<"Error: Call to default constructor of class Dupl_SubF_3_Rate_NoState_model\n";};
	Dupl_SubF_3_Rate_NoState_model (Exchange *cexchange, Tree *ctree, Clade *cgenomes, WGD_Data *chomologs);
	void describe_results();
	void num_params_model();

protected:
	void intialize_parameters (double par[], PARAM_TYPE types[]);
};


class Dupl_SubF_3_Rate_AllStates_NoState_model : public Dupl_SubF_3_Rate_AllStates_Base_model, public Dupl_NoState_Base_model
{
public:
	Dupl_SubF_3_Rate_AllStates_NoState_model () {cerr<<"Error: Call to default constructor of class Dupl_SubF_3_Rate_NoState_model\n";};
	Dupl_SubF_3_Rate_AllStates_NoState_model (Exchange *cexchange, Tree *ctree, Clade *cgenomes, WGD_Data *chomologs);
	void describe_results();
	void num_params_model();

protected:
	void intialize_parameters (double par[], PARAM_TYPE types[]);
};


class Dupl_2_Rate_NoSubF_NoState_model : public Dupl_2_Rate_NoSubF_Base_model, public Dupl_NoState_Base_model
{
public:
	Dupl_2_Rate_NoSubF_NoState_model () {cerr<<"Error: Call to default constructor of class Dupl_2_Rate_NoSubF_NoState_model\n";};
	Dupl_2_Rate_NoSubF_NoState_model (Exchange *cexchange, Tree *ctree, Clade *cgenomes, WGD_Data *chomologs);
	void describe_results();
	void num_params_model();

protected:
	void intialize_parameters (double par[], PARAM_TYPE types[]);
};


class Dupl_2_Rate_NoSubF_AllStates_NoState_model : public Dupl_2_Rate_NoSubF_AllStates_Base_model, public Dupl_NoState_Base_model
{
public:
	Dupl_2_Rate_NoSubF_AllStates_NoState_model () {cerr<<"Error: Call to default constructor of class Dupl_SubF_3_Rate_NoState_model\n";};
	Dupl_2_Rate_NoSubF_AllStates_NoState_model (Exchange *cexchange, Tree *ctree, Clade *cgenomes, WGD_Data *chomologs);
	void describe_results();
	void num_params_model();

protected:
	void intialize_parameters (double par[], PARAM_TYPE types[]);
};


class Dupl_NoSubF_2_Rate_model : public Dupl_2_Rate_NoSubF_Base_model, public Dupl_State_Base_model
{
	public:
	Dupl_NoSubF_2_Rate_model () {};
	Dupl_NoSubF_2_Rate_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree);
	void describe_results();
	void num_params_model();

protected:
	void intialize_parameters (double par[], PARAM_TYPE types[]);
};

class Dupl_SubF_Only_model : public Dupl_SubF_3_Rate_Base_model, public Dupl_State_Base_model
{
public:
	Dupl_SubF_Only_model () {};
	Dupl_SubF_Only_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree);
	void describe_results();
	void num_params_model();

protected:
	void intialize_parameters (double par[], PARAM_TYPE types[]);
};


class Dupl_SubF_Only_NoState_model : public Dupl_SubF_Only_Base_model, public Dupl_NoState_Base_model
{
public:
	Dupl_SubF_Only_NoState_model () {cerr<<"Error: Call to default constructor of class Dupl_SubF_All_States_NoState_model\n";};
	Dupl_SubF_Only_NoState_model (Exchange *cexchange, Tree *ctree, Clade *cgenomes, WGD_Data *chomologs);
	void describe_results();
	void num_params_model();

protected:
	void intialize_parameters (double par[], PARAM_TYPE types[]);
};






class Dupl_SubF_All_States_model : public Dupl_Subf_All_States_Base_model, public Dupl_State_Base_model
{
public:
	Dupl_SubF_All_States_model () {};
	Dupl_SubF_All_States_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree);
	void describe_results();
	void num_params_model();
	void get_site_state_probs(double **&prob_array, int taxa_id);

protected:
	void intialize_parameters (double par[], PARAM_TYPE types[]);
};


class Dupl_SubF_All_States_NoState_model : public Dupl_Subf_All_States_Base_model, public Dupl_NoState_Base_model
{
public:
	Dupl_SubF_All_States_NoState_model () {cerr<<"Error: Call to default constructor of class Dupl_SubF_All_States_NoState_model\n";};
	Dupl_SubF_All_States_NoState_model (Exchange *cexchange, Tree *ctree, Clade *cgenomes, WGD_Data *chomologs);
	void describe_results();
	void num_params_model();
//	void get_site_state_probs(double **&prob_array, int taxa_id);

protected:
	void intialize_parameters (double par[], PARAM_TYPE types[]);
};


class Dupl_Fix_Parallel_NoState_model : public Dupl_Fix_Parallel_Base_model, public Dupl_NoState_Base_model
{
public:	
	Dupl_Fix_Parallel_NoState_model () {cerr<<"Error: Call to default constructor of class Dupl_Fix_Parallel_NoState_model\n";};
	Dupl_Fix_Parallel_NoState_model (Exchange *cexchange, Tree *ctree, Clade *cgenomes, WGD_Data *chomologs);
			
	void describe_results();
	void num_params_model();

protected:
	void intialize_parameters (double par[], PARAM_TYPE types[]);
};


class Dupl_Fix_Parallel_SubF_NoState_model : public Dupl_Fix_Parallel_SubF_Base_model, public Dupl_NoState_Base_model
{
public:	
	Dupl_Fix_Parallel_SubF_NoState_model () {cerr<<"Error: Call to default constructor of class Dupl_Fix_Parallel_NoState_model\n";};
	Dupl_Fix_Parallel_SubF_NoState_model (Exchange *cexchange, Tree *ctree, Clade *cgenomes, WGD_Data *chomologs);
			
	void describe_results();
	void num_params_model();

protected:
	void intialize_parameters (double par[], PARAM_TYPE types[]);
};


class Dupl_Slow_Loss_Con_Fix_model : public Dupl_Slow_Loss_Con_Fix_Base_model, public Dupl_State_Base_model
{
public:
	Dupl_Slow_Loss_Con_Fix_model () {};
	Dupl_Slow_Loss_Con_Fix_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree);
	void describe_results();
	void num_params_model();

protected:
	void intialize_parameters (double par[], PARAM_TYPE types[]);

};

class Dupl_Slow_Loss_Con_Fix_NoState_model : public Dupl_Slow_Loss_Con_Fix_Base_model, public Dupl_NoState_Base_model
{
public:
	Dupl_Slow_Loss_Con_Fix_NoState_model () {cerr<<"Error: Call to default constructor of class Dupl_Slow_Loss_Con_Fix_NoState_model\n";};
	Dupl_Slow_Loss_Con_Fix_NoState_model (Exchange *cexchange, Tree *ctree, Clade *cgenomes, WGD_Data *chomologs);
			
	void describe_results();
	void num_params_model();

protected:
	void intialize_parameters (double par[], PARAM_TYPE types[]);
};



class SNP_model : public Like_model
{
public:
	SNP_model () {};
	SNP_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree);
  
	void describe_results();  
	void num_params_model();
  
 
  	double root_freq(int site);	
	double get_basefreq(int base, Branch *tree_branch) {return(0.0);};
	double get_basefreq(int base, int codon_pos, Branch *tree_branch) {return(0.0);};
	double get_expect_sub_site(double brlen)  {return(brlen);};
	double  optimize_a_branch(double brnlen);
	double newton_branch_opt(double tol);

protected:
	void initialize_arrays();
	void intialize_parameters (double par[], PARAM_TYPE types[]);
	void partial_prob_w_rate(int locale, Branch *lsib, int rate_num, Branch *stop_id);
	double get_tip_prob(int locale, Branch *taxa, int rate_num, int cond_state);
	void calc_transprobs(Branch *taxa, int rate_num);
	double find_ut(Branch *taxa);
	void store_conprobs();
};


#endif
