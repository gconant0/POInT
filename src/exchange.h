//exchange.h
//Copyright 1999-2002 Gavin Conant 

//The exchange object hold global parameters that may be needed
//by numerous other classes.  Although models of evolution are
//set in the Like_model classes, the exchange object holds information
//about site-wise rates, base frequencies, trs/trv ratio etc.


#ifndef ___EXCHANGE_H___
#define ___EXCHANGE_H___
#include "gen_dna_funcs.h"

#ifdef MPI_CODE_VERSION
#include "mpi_info.h"
#endif

#define NUM_AA_PROPS 5


enum SITE_RATE_TYPE {SINGLE_RATE, ARBITRARY_SITE_RATES, CODON_RATES, SITE_SPECFIC_RATE_PROB, GENERIC_RATE_PROB};

class Exchange {
 public:
    Amino_acid_group *curr_groups;
    Exchange();
    Exchange& operator=(Exchange &assign_from);
    int get_num_taxa()                              {return(num_taxa);};
    int get_num_sites()                             {return(num_sites);};
    int get_num_codons()                            {return(num_codons);};
    int get_num_localities()                        {return(this->*num_localities);};
    int get_num_branches()                          {return(num_branches);};
    int get_num_rates()                             {return(num_rates);};
    int get_num_params()                            {return(num_params);};
    int get_condlike_size()                         {return(condlike_size);};
    int get_prop_index_num(AA_PROPERTIES prop);
    int get_num_aa_props()						  {return(NUM_AA_PROPS);};
    int get_num_live_aa_props()					  {return(num_live_properties);};
    int get_site_rate_num(int site)				  {return((this->*site_rate_num_func)(site));};
    int get_branch_kaks1(int brn_num);
    int get_num_p_non_syn()                         {return (num_p_nonsyn);};
    int get_num_matrices()						  {return (num_matrices);};
    int get_num_p_non_syn_used();
    int get_num_p_inter_group_used();
    int get_num_aa_props_used(int prop);
    int get_num_nonsyn_params()					  {return(this->*num_nonsyn_params);};
    int get_num_nonsyn_patterns()					  {return(this->*num_patterns);};
    int get_rate_aa_prop_rate(int rate_num, AA_PROPERTIES prop);
    int get_rate_aa_prop_rate(int rate_num, int prop_num);
    int num_unallowed_states()						  {return(num_unallowed_SNP_states);};
    int get_WGX_depth()                              {return(WGX_depth);};
    int get_max_breaks()                             {return(max_breaks);};
	

  AA_PROPERTIES get_prop_num_n(int n);
  
  BOOL non_standard_code()                        {return(nonstandard_gencode);};
  BOOL modeling_codons()                          {return(model_codons);};
  BOOL pre_optimizing()                           {return(pre_optimize);};
  BOOL aa_classes_modeled()                       {return(aa_class_model);};
  BOOL branch_basefreqs()                         {return(branch_freqs);};
  BOOL have_data()                                {return(have_dataset);};
  BOOL p_non_syn_fixed(int pns_num)               {return(p_non_syn_fix[pns_num][0]);};
  BOOL p_non_syn_fixed(int pns_num, int rate_num) {return(p_non_syn_fix[pns_num][rate_num]);};
  BOOL p_inter_group_fixed(int pns_num)           {return(p_inter_group_fix[pns_num][0]);};
  BOOL p_inter_group_fixed(int pns_num, int rate_num)           {return(p_inter_group_fix[pns_num][rate_num]);};
  BOOL using_codon_position_rates();               
  BOOL using_arbitrary_site_rates(); 
  BOOL using_generic_site_rate_probs();
  BOOL using_codon_basefreqs()                    {return(codon_freqs);};
  BOOL fixed_basefreq()                           {return(basefreqs_fixed);}; 
  BOOL is_mol_clock_3()                           {return(three_tree_mol_clock);}; 
  BOOL have_branch_kaks1()                        {return(is_kaks1);};      
  BOOL have_site_rate_probs();
  BOOL aa_prop_allowed(AA_PROPERTIES prop);		  
  BOOL aa_prop_allowed(int prop);				  
  BOOL is_p_non_syn_used(int pns_num);			 
  BOOL is_p_inter_group_used(int pns_num);		  
  BOOL is_aa_prop_used(int pns_num, AA_PROPERTIES prop);
  BOOL is_aa_prop_used(int pns_num, int prop);
  BOOL is_matrix_coeff_used(int pns_num, int matrix);
  BOOL are_optimizing_rate_props()				  {return(opt_rate_probs);};
  BOOL allow_pos_LCAP_weights()					  {return(allow_LCAP_pos_weights);};
  BOOL is_rooted_tree()							  {return(rooted_tree);};
  BOOL zero_len_brns_fixed()					  {return(fix_zero_brns);};
  BOOL use_track_gaps_as_missing()				      {return(track_gaps_as_missing);};
  BOOL full_conprobs()							  {return(use_full_conprobs);};
  BOOL likelihood_is_scaled()					  {return(scale_likelihood);};
    BOOL use_guessed_order()                        {return(guess_order);};
    BOOL count_break_pos()                          {return(use_break_pos);};
    BOOL use_all_taxa()                             {return(all_taxa_used);};
    BOOL is_taxa_used(int taxa_num);

  CLOCK_TYPE get_clock_type()                          {return(clock_type);};
  SITE_RATE_TYPE get_site_rate_type()				   {return(rate_type);};
  
  double return_basefreq(int freqnum);
  double return_codon_basefreq(int codon_pos, int base);
  double get_trs_trv()                            {return(kappa);};
  double get_obs_trs_trv()                        {return(obs_trs_trv);};
  double get_rate(int rate_num)                   {return(rates[rate_num]);};
  double get_site_rate(int site)                  {return((this->*site_rate_func)(site));};
  double get_site_rate_prob(int site, int rate);   
  double get_p_non_syn();
  double get_p_non_syn(int pns_num);
  double get_p_non_syn(int pns_num, int rate_num);
  double get_p_inter_group(int pns_num, int rate_num);
  double get_p_inter_group(int pns_num);
  double get_aa_prop_fac(int pns_num, 
			 AA_PROPERTIES prop);
  double get_aa_prop_fac(int pns_num, 
			 int prop);
  double get_aa_prop_fac(int pns_num, int rate_num, AA_PROPERTIES);
  double get_aa_prop_fac(int pns_num, int rate_num, int property);
  double get_matrix_coeff(int coeff_num, int pns_num);
  double get_saved_lnL()                          {return(saved_lnL);};
  double get_rate_prob(int rate)				  {return(rate_probs[rate]);};
  double get_dupl_parallel_rate()				  {return(dupl_parallel_rate);};
  double get_dupl_fix_rate()					  {return(dupl_fix_rate);};
  double get_strand_switch_prob()				  {return(strand_switch_prob);};
  double get_fix_rate_scale()					  {return(fix_rate_scale);};
  double get_loss_rate_scale()					  {return(loss_rate_scale);};
  double get_fix_loss_rate()					  {return(fix_loss_rate);};
  double get_snp_gains_to_losses()				  {return(snp_gains_to_losses);};
  double get_snp_loss_to_switch()				  {return(snp_loss_to_switch);};
  char* get_treefile()                            {return(treefile);};
  char* get_datafile()                            {return(Exchange::datafile);};
  char* get_gencode_file()                        {return(Exchange::gencode_file);};
  
  LKMODEL get_model()                             {return(Exchange::current_model);};
  DATAFORMAT get_dataformat()                     {return(current_dataformat);};
  GENETIC_CODE get_genetic_code()                 {return(curr_genetic_code);};  
  Sequence_dataset * get_dataset()                {return(current_dataset);};
	SNP_STATE get_unallowed_state(int pos)        {return (unallowed_states[pos]);};
  
  void set_num_taxa(int taxa);
  void set_num_sites(int sites);
  void set_num_rates(int rates);
  void set_num_params(int ps)                     {num_params=ps;};
  void set_basefreqs (double freq, int freq_num)  {basefreqs[freq_num]=freq;};
  void set_trs_trv (double tvratio)               {kappa=tvratio;};
  void set_obs_trs_trv(double obs_tr)             {obs_trs_trv=obs_tr;};
  void set_rate(int rate_num, double rate)        {rates[rate_num]=rate;};
  void set_dupl_parallel_rate(double rate)		  {dupl_parallel_rate=rate;};
  void set_dupl_fix_rate(double rate)			  {dupl_fix_rate=rate;};
  void set_strand_switch_prob(double p)			  {strand_switch_prob=p;};
  void set_fix_rate_scale(double p)				  {fix_rate_scale=p;};
  void set_loss_rate_scale(double p)		      {loss_rate_scale=p;};
  void set_fix_loss_rate(double p)				  {fix_loss_rate=p;};
  void set_snp_gains_to_losses(double p)		  {snp_gains_to_losses=p;};
  void set_snp_loss_to_switch(double p)		      {snp_loss_to_switch=p;};
  void set_branch_kaks1(int branch);
  void set_p_non_syn(int pns_num, double pns);
  void set_p_non_syn(int pns_num, int rate, double pns);
  void set_matrix_coeff(int coeff_num, int 
	  pns_num, double value)					  {matrix_coeff[pns_num][coeff_num]=value;};
  void set_p_inter_group(int pns_num, double pgroup);
  void set_p_inter_group(int pns_num, int rate_num, double pgroup);
  void set_aa_prop_fact(int pns_num, AA_PROPERTIES property, 
			double value);
  void set_aa_prop_fact(int pns_num, int property, 
			double value);
  void set_aa_prop_fact(int pns_num, int rate_num, AA_PROPERTIES property, 
			double value);
  void set_aa_prop_fact(int pns_num, int rate_num, int property, 
			double value);
  void set_scaling_initial(int pns_num, double p_ns);
  void set_saved_lnL(double lnL)                  {saved_lnL=lnL;};
  void set_pre_optimize(BOOL opt)                 {pre_optimize=opt;};
  void set_aa_class_modeling(BOOL modeling)       {aa_class_model=modeling;};
  void set_branch_basefreqs(BOOL use)             {branch_freqs=use;};
  void set_num_p_nonsyn(int num, double **pns); 
  void set_num_p_nonsyn(int num);
  void set_num_p_nonsyn(int n_pns, int n_patterns);
  void set_num_matrices(int num);
  void set_max_breaks(int mbreaks)              {max_breaks=mbreaks;};
  void set_have_data(BOOL set)                    {have_dataset=set;};
  void set_mol_clock_3()                          {three_tree_mol_clock=TRUE;};
  void set_clock_type(CLOCK_TYPE type)            {clock_type=type;};
  void fix_p_non_syn(double **p_ns); 
  void fix_p_non_syn(int pns_num, double p_ns);
  void fix_p_non_syn(int pns_num, int rate, double p_ns);
  void fix_p_inter_group(int pns_num, double p_intgrp[]);
  void fix_p_inter_group(int pns_num, int rate_num, double p_intgrp);
  void fix_p_inter_group(int pns_num, int rate_num);
  void fix_zero_brn_lens()						  {fix_zero_brns=TRUE;};
  void set_treefile(const char filename[500])           {strcpy(treefile, filename);};
  void set_datafile(const char filename[500])           {strcpy(datafile, filename);};
  void set_gencode_file(char filename[500])       {strcpy(gencode_file, filename);};
  void set_model(LKMODEL curmodel);
  void set_model(LKMODEL curmodel, int num_states);
  void set_dataformat(DATAFORMAT format)          {current_dataformat=format;};
  void set_dataset(Sequence_dataset *the_data)    {current_dataset=the_data;};
  void set_genetic_code(GENETIC_CODE code)        {curr_genetic_code=code;};
  void disable_property(AA_PROPERTIES prop);
  void set_p_non_syn_used(int pns_num, BOOL val);
  void set_p_inter_group_used(int pns_num, BOOL val);
  void set_aa_property_used(int pns_num, AA_PROPERTIES prop, BOOL val);
  void set_aa_property_used(int pns_num, int prop, BOOL val);
  void set_matrix_coeff_used(int pns_num, int matrix_num, BOOL val);
  void use_prob_site_rates(int num_rates);
  void set_site_rate_prob(int site, int rate, double rate_prob);
  void optimize_rate_probs()					  {opt_rate_probs=TRUE;};
  void set_use_full_conprobs()					  {use_full_conprobs=TRUE;};
  void disallow_LCAP_pos_weights()				  {allow_LCAP_pos_weights=FALSE;};
  void treat_track_gaps_as_missing(BOOL val)	  {track_gaps_as_missing = val;};
  void use_single_property(AA_PROPERTIES prop);
  void use_single_property(int prop_num);
  void use_only_hetero_internal();
  void set_WGX_depth(int d)                     {if (d>1) {WGX_depth=fabs(d);} else {WGX_depth=2;}};
  void set_guess_order()                      {guess_order=TRUE;};
    void set_use_break_pos()                    {use_break_pos=TRUE;};

  ~Exchange();
  void set_use_codon_position_rates();
  void set_use_arbitrary_site_rates(int *initial_rates);
  void set_site_rate(int site, int rate_num);
  void set_use_generic_site_rates();
  void set_generic_site_rate_prob(int rate, double prob);
  void set_use_codon_basefreqs(double newfreqs[3][4]);
  void fix_basefreq()                             {basefreqs_fixed=TRUE;};
    void mark_taxa_unused(int taxa_num);
    
#ifdef MPI_CODE_VERSION
  Message_pass_interface* get_mpi_info() {return(mpi_inter);};
  void set_mpi_info(Message_pass_interface *info) {mpi_inter=info;};
#endif
#ifdef _OPEN_MP_VERSION_
	int get_num_open_mp_threads()					  {return(num_open_mp_threads);}
#endif

 protected:
  char treefile[500], datafile[500], gencode_file[500];
  int num_taxa, num_sites, num_branches, num_rates,  *branch_kaks1,
     num_params, num_codons, *site_rates, Exchange::*num_localities, 
	 condlike_size, num_p_nonsyn, num_matrices, num_live_properties, Exchange::*num_nonsyn_params,
	 Exchange::*num_patterns, num_nonsyn_patterns, standard_model_nonsyn_params, **aa_rate_index, num_unallowed_SNP_states, WGX_depth, max_breaks;
#ifdef _OPEN_MP_VERSION_
	int num_open_mp_threads;
#endif
	int (Exchange::*site_rate_num_func)(int site);  
  double (Exchange::*site_rate_func)(int site);
  double (Exchange::*site_rate_prob_func)(int site, int rate);
  double *rates, basefreqs[4], codon_basefreqs[3][4], obs_trs_trv, kappa, **p_non_syn, **p_inter_group, 
    ***aa_properties, **matrix_coeff, saved_lnL, **site_rate_probs, *rate_probs, 
	dupl_parallel_rate, dupl_fix_rate, strand_switch_prob, fix_rate_scale, 
	loss_rate_scale, fix_loss_rate, snp_gains_to_losses, snp_loss_to_switch;
  BOOL  nonstandard_gencode, model_codons, codon_freqs, basefreqs_fixed, 
    pre_optimize, aa_class_model, branch_freqs, have_dataset, **p_inter_group_fix, **p_inter_group_same, **p_non_syn_fix, 
	three_tree_mol_clock, ks_fixed, is_kaks1, aa_property_allowed[NUM_AA_PROPS+1], *p_non_syn_used, 
	*p_inter_group_used, **aa_prop_used, **matrix_coeff_used, opt_rate_probs, allow_LCAP_pos_weights,
	rooted_tree, fix_zero_brns, track_gaps_as_missing, use_full_conprobs, scale_likelihood, guess_order, use_break_pos, all_taxa_used, *taxa_used;
  SITE_RATE_TYPE rate_type;
  LKMODEL current_model;
  DATAFORMAT current_dataformat;
  GENETIC_CODE curr_genetic_code;
  CLOCK_TYPE clock_type;
	SNP_STATE unallowed_states[5];
  Sequence_dataset *current_dataset;

#ifdef MPI_CODE_VERSION
  Message_pass_interface *mpi_inter;
#endif

  //Functions
  void set_rates();
  int single_rate_site_rate_num(int site)             {return(0);};
  int codon_positions_site_rate_num(int site);
  int arbitrary_site_rate_num(int site);
  double single_rate_site_rates(int site)             {return(rates[0]);};
  double codon_positions_site_rates(int site);
  double arbitrary_site_rates(int site);
  double arbitrary_site_rate_prob(int site, int rate);
  double codon_position_site_rate_prob(int site, int rate);
  double singe_rate_site_rate_prob(int site, int rate);
  double site_specific_site_rate_prob(int site, int rate);
  double generic_rate_prob(int site, int rate);



};


#endif
