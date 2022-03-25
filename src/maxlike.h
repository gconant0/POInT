//Copyright 1999-2002 Gavin Conant


#ifndef ___MAXLIKE_H___
#define ___MAXLIKE_H___

#include "tree.h"
#include "read_seq.h"
#include "exchange.h"
#include "gen_dna_funcs.h"
#include <iostream>

using namespace::std;

#ifdef MPI_CODE_VERSION
#include "mpi_info.h"
extern Message_pass_interface mpi_info;
#endif


#define NEG_INF -1.0e300
#define LN_HALF -.69314718056

#define BRN_TOL 1e-5

#ifndef FLOAT_TOL
#define FLOAT_TOL 1e-9
#endif

#ifndef LIKE_SCALE
#define LIKE_SCALE 1e240
#define SQRT_LIKE_SCALE 1e120
#endif


unsigned long RunTime();  
extern long pre, post;

class Like_model
//This is the base class for all likelihood models, codon and nucleotide
//Constructors become challenging with multiple inheritance, so the void function
//assemble takes the place of the Like_model constructor
//For multiple inheritance to work correctly there cannot be pure virtual functions in this class
//(i.e. even virtual functions must have a function definition).  The virtual function definitions are
//stabs that print error messages, as, if inheritance have occured properly, they should never be called
{
 public:
  //variables 
  double **q_matrix, **t_matrix;           //Shares data with the Linear_algebra class and the lapack linear algebra routines
  BOOL just_found_opt, swap_brn_opt;

  //Functions
  Like_model ();
  void assemble (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree);
  void recalculate_transprobs();
  void list_transprobs(Branch *taxa);
  void list_rates();
  void reinit_params();                         
  void send_params(double[], int offset);
  void set_save_k_common(double save)         {save_k_common=save;};
  void set_save_k_dupl(double save)           {save_k_dupl=save;}; 
  PARAM_TYPE get_param_type(int param_num);
  virtual double find_appropriate_ln_like();
  double find_ln_like_w_rates();
  double find_two_taxa_ln_like();
  double find_ln_like_single_rate();
  double find_ln_like_fixed_site_rates();
  double find_lnL_on_tree(Tree *ctree);
  void return_lnL (int start, int end, int rate_num, double *new_rate, double new_sitelnLs[]);
  double min_lnL (double params[], int offset);
  void min_single_param(PARAM_TYPE parameter, double tol);
  double optimize_single_param(double trs_trv);
  double newton_branch_opt();
  virtual double newton_branch_opt(double tol);
  void setup_branch_opt(Branch *branch);
  double newton_opt_a_branch(Branch* branch, double tol);
  double newton_single_iter(Branch *branch, double step_size, double tol);
  void optimize_branches ();
  virtual double optimize_a_branch (double new_ut);
  void branch_opt_like_and_deriv(double ut, double &likelihood, double &likelihood_prime, double &likelihood_double_prime, BOOL first);
  double branch_opt_func(double new_part_brn);
  void set_rate (double new_rate);
  void set_rate_by_sites (int pos, int rate);
  void get_rate_posts(int site, double *posts);
  virtual void set_expect_subs();
  virtual double get_expect_sub_site(double brlen);

  //get_Basefreq is overloaded for codon-basefrequencies
  virtual double get_basefreq(int base, Branch *tree_branch);
  virtual double get_basefreq(int base, int codon_pos, Branch *tree_branch);

  virtual double root_freq (int base);
  virtual double get_obs_trs_trv_ratio(); 
  virtual void calc_first_derivatives(Branch *taxa, double ut);
  virtual void calc_second_derivatives(Branch *taxa, double ut);
  virtual void calc_first_derivatives(Branch *taxa, double ut, int rate);
  virtual void calc_second_derivatives(Branch *taxa, double ut, int rate);
  virtual void print_model ();
  virtual void describe_results();
  virtual void describe_nuc_results();
  virtual int rate_param_size();
  virtual void set_ratios(Branch *taxa) {cerr<<"Wrong get_ratio\n";};
  virtual void set_arb_param(int param_num, double param_val) {cerr<<"Error: Cannot use arbitrary parameters in base Like_model\n";};
  virtual ~Like_model();
  
 protected:   //Protected members for inheritance
  
  //Variables
  int  bf_end, curr_brn1, curr_brn2, old_parent, brn_start, start_swap, *brn_order, pns_start, 
	  *prop_pns_index, *prop_rate_index, num_aa_prop_coeff, *num_aa_props_per_rate, num_pns_coeff, num_pi_coeff, matcoeff_start, 
	  last_codonmat, *matrix_index, rate_prob_param_index, *brn_index, save_conprobs_size, last_brn_remain;

  double mu,  *params,  p_stop,  
    pur_pyr_split, a_g_split, c_t_split, 
    ***brncondprobs1, ***brncondprobs2, ***save_conprobs, save_k_common, save_k_dupl,
    last_lnL, last_like_prime, last_like_double_prime,  left_to_tip, 
    *site_lnLs, right_to_tip;

  BOOL first_brnopt;
  PARAM_TYPE *param_types, single_param_id;
  AA_PROPERTIES *prop_index;
  Branch *brn_op_cbrn;
  Tree *curr_tree;
  Sequence_dataset *curr_data;
  Exchange *curr_exchange;
  Constrain_Param_Lookup *lookup_p_tree;
#ifdef MPI_CODE_VERSION
  Message_pass_interface *mpi_interface;
#endif
  //Virtual functions implemented in derived classes
  virtual void calc_codon_matrices (int rate_num, Branch  *taxa);
  virtual void calc_transprobs (Branch *taxa, int rate_num);
  virtual void num_params_model();
  virtual void initialize_arrays();
  virtual void intialize_parameters (double par[], PARAM_TYPE types[]);
  virtual double q_matrix_entry (int start_base, int end_base, int codon_pos, Branch *tree_branch);
  virtual double prob_nonsyn(int start[3], int end[3], Branch *taxa, int rate_num);
  virtual double find_ut(Branch *taxa);
  virtual double get_trs_trv();
  virtual BOOL change_rate(int rate_num, double *rate_info);
  //Likelihood calculating functions
  virtual long double prob_w_rate (int locale, int rate_num);
  //This is a virtual function to allow ambigutity states in the Dupl_model class of models
  virtual void partial_prob_w_rate(int locale, Branch *lsib, int rate_num, Branch *stop_id);

  void copy_transprobs(Branch *from, Branch *to, int rate_num);

  //Utility functions
  void array_copy(int array1[], int array2[], int size);
  int num_rate_prob_params();
  int initialize_rate_prob_params(int start_param_index, double *par, PARAM_TYPE *types);
};


//The next four classes inplement the nucleotide substitution models indicated
//They can be inherited either by nucleotide models (nucleotide_like.h), where
//they specify completely the model of evolution used, or by the codon models,
//where they specify only the nucleotide substitution pattern

class JC_model : virtual public Like_model             
//Juke-Cantor 1969 model  
{
 public:
  double get_expect_sub_site(double brlen);
  void print_model();
  void describe_nuc_results();
 protected:
  double get_basefreq(int base, Branch *tree_branch);
  double get_basefreq(int base, int codon_pos, Branch *tree_branch);
  double get_trs_trv();
  double q_matrix_entry (int start_base, int end_base, int codon_pos, Branch *tree_branch);
  double find_ut(Branch *taxa);
};



class K2P_model : virtual public Like_model
//Kimura 2-Paramter model (1980)
{ 
 public: 
  double get_expect_sub_site(double brlen);
  void print_model();
  void describe_nuc_results();
 protected:
  double get_basefreq(int base, Branch *tree_branch);
  double get_basefreq(int base, int codon_pos, Branch *tree_branch);
  double get_trs_trv();
  double q_matrix_entry (int start_base, int end_base, int codon_pos, Branch *tree_branch);
  double find_ut(Branch *taxa);
};



class HKY_model : virtual public Like_model
//Hasegawa, Kishino, Yano model (1985)
{
 public:
  double get_expect_sub_site(double brlen);
  void print_model();
  void describe_nuc_results();
 protected:
  double get_basefreq(int base, Branch *tree_branch);
  double get_basefreq(int base, int codon_pos, Branch *tree_branch);
  double get_trs_trv();
  double q_matrix_entry (int start_base, int end_base, int codon_pos, Branch *tree_branch);
  double find_ut(Branch *taxa);

};



class GG_98_model : virtual public Like_model
//Galtier and Gaut model (1998)--may have other names
{
 public:
  void print_model();
  void describe_nuc_results();
 protected:
  double get_basefreq(int base, Branch *tree_branch); 
  double get_basefreq(int base, int codon_pos, Branch *tree_branch);
  double get_trs_trv();
  double q_matrix_entry (int start_base, int end_base, int codon_pos, Branch *tree_branch);
  double find_ut(Branch *taxa);

};


#endif





