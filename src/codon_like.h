//Copyright 1999-2005 Gavin Conant

#include "gen_code.h"
#include "lin_alg.h"
#include "gen_dna_funcs.h"
#include "maxlike.h"
#include "score_matrix.h"

#ifndef ___CODON_LIKE_H___
#define ___CODON_LIKE_H___



class  Codon_model : virtual public Like_model
//Class with implements common features of codon-based models
{
 public:
  double root_freq(int site); 
  void set_expect_subs();
  void set_ratios(Branch *taxa);
  double get_obs_trs_trv_ratio();
  virtual ~Codon_model();

 protected:
  double prop_syn, prop_nsyn, neu_prop_syn, neu_prop_nsyn, stop_excess;
  Genetic_code *curr_code;
  Linear_Algebra *curr_lin_alg;

  void calc_codon_matrices (int rate_num, Branch *taxa);
  void calc_first_derivatives(Branch *taxa, double ut);
  void calc_second_derivatives(Branch *taxa, double ut);
  void calc_first_derivatives(Branch *taxa, double ut, int rate);
  void calc_second_derivatives(Branch *taxa, double ut, int rate);
  void initialize_arrays();
  void calc_transprobs (Branch *taxa, int rate_num);
};

//Implement specific codon-models, but do not specify model
//of nucleotide substitution
class MG_94_model :  virtual public Codon_model
{
 public:
  void describe_results();
  int rate_param_size();
  ~MG_94_model() {};
 protected:
  BOOL change_rate(int rate_num, double *rate_info);  
  double prob_nonsyn(int start[3], int end[3], Branch *taxa, int rate_num);  
  void set_num_pns_used();
  int set_pns_id_indices(double par[], PARAM_TYPE types[]);
};


class C_00_model : virtual public Codon_model
{
 public:
  void describe_results();
  int rate_param_size();
 protected:
  BOOL change_rate(int rate_num, double *rate_info);
  double prob_nonsyn(int start[3], int end[3], Branch *taxa,  int rate_num);
  void set_num_pns_used();
  int set_pns_id_indices(double par[], PARAM_TYPE types[]);
};


class LCAP_model : virtual public Codon_model
{
 public:
  void describe_results();
  int rate_param_size();
 protected:
  BOOL change_rate(int rate_num, double *rate_info);
  double prob_nonsyn(int start[3], int end[3], Branch *taxa,  int rate_num);
  void set_num_aa_props_used();
  int set_aa_prop_id_indices(double par[], PARAM_TYPE types[]);
};



class MultiMatrix_model : virtual public Codon_model
{
public:
  void describe_results();
  int rate_param_size();
 protected:
  BOOL change_rate(int rate_num, double *rate_info);
  double prob_nonsyn(int start[3], int end[3], Branch *taxa, int rate_num);
  AA_matrix **the_matrices;

  void set_num_pns_used();
  int set_pns_id_indices(double par[], PARAM_TYPE types[]);

};


//The following classes use multiple inheritance to create specific
//instances of codon-based models with specific models of nucleotide
//substitution


class MG_94_JC_model : public JC_model, public MG_94_model
{
 public:
  MG_94_JC_model ();
  MG_94_JC_model(Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree, Genetic_code *ccode);
  ~MG_94_JC_model() {delete curr_lin_alg;};
 protected:
  void num_params_model();
  void intialize_parameters (double par[], PARAM_TYPE types[]);
};


class MG_94_K2P_model : public K2P_model, public MG_94_model
{
 public:
  MG_94_K2P_model ();
  MG_94_K2P_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree, Genetic_code *ccode);
  ~MG_94_K2P_model () {delete curr_lin_alg;};
 protected:
  void num_params_model();
  void intialize_parameters (double par[], PARAM_TYPE types[]);
};


class MG_94_HKY_model : public HKY_model, public MG_94_model
{
 public:
  MG_94_HKY_model ();
  MG_94_HKY_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree, Genetic_code *ccode);
  ~MG_94_HKY_model ()   {delete curr_lin_alg;};
 protected:
  void num_params_model();
  void intialize_parameters (double par[], PARAM_TYPE types[]);
};


class MG_94_GG_98_model : public GG_98_model, public MG_94_model
{
 public:
  MG_94_GG_98_model ();
  MG_94_GG_98_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree, Genetic_code *ccode);
  ~MG_94_GG_98_model () {delete curr_lin_alg;};
 protected:
  void num_params_model();
  void intialize_parameters (double par[], PARAM_TYPE types[]);
};


class C_00_JC_model : public JC_model, public C_00_model
{
 public:
  C_00_JC_model ();
  C_00_JC_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree, Genetic_code *ccode);
  ~C_00_JC_model () {delete curr_lin_alg;};
protected:
  void num_params_model();
  void intialize_parameters (double par[], PARAM_TYPE types[]);
};


class C_00_K2P_model : public K2P_model, public C_00_model
{
 public:
  C_00_K2P_model ();
  C_00_K2P_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree, Genetic_code *ccode);
  ~C_00_K2P_model () {delete curr_lin_alg;};
 protected:
  void num_params_model();
  void intialize_parameters (double par[], PARAM_TYPE types[]);
};


class C_00_HKY_model : public HKY_model, public C_00_model
{
 public:
  C_00_HKY_model ();
  C_00_HKY_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree, Genetic_code *ccode);
  ~C_00_HKY_model () {delete curr_lin_alg;};
 protected:
  void num_params_model();
  void intialize_parameters (double par[], PARAM_TYPE types[]);
};


class C_00_GG_98_model : public GG_98_model, public C_00_model
{
 public:
  C_00_GG_98_model ();
  C_00_GG_98_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree, Genetic_code *ccode);
  ~C_00_GG_98_model () {delete curr_lin_alg;};
 protected:
  void num_params_model();
  void intialize_parameters (double par[], PARAM_TYPE types[]);
};


class LCAP_JC_model : public JC_model, public LCAP_model
{
 public:
  LCAP_JC_model ();
  LCAP_JC_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree, Genetic_code *ccode);
  ~LCAP_JC_model () {delete curr_lin_alg;};
 protected:
  void num_params_model();
  void intialize_parameters (double par[], PARAM_TYPE types[]);
};


class LCAP_K2P_model : public K2P_model, public LCAP_model
{
 public:
  LCAP_K2P_model ();
  LCAP_K2P_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree, Genetic_code *ccode);
  ~LCAP_K2P_model () {delete curr_lin_alg;};
 protected:
  void num_params_model();
  void intialize_parameters (double par[], PARAM_TYPE types[]);
};


class LCAP_HKY_model : public HKY_model, public LCAP_model
{
 public:
  LCAP_HKY_model ();
  LCAP_HKY_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree, Genetic_code *ccode);
  ~LCAP_HKY_model () {delete curr_lin_alg;};
 protected:
  void num_params_model();
  void intialize_parameters (double par[], PARAM_TYPE types[]);
 // double find_ut(Branch *taxa);
};


class LCAP_GG_98_model : public GG_98_model, public LCAP_model
{
 public:
  LCAP_GG_98_model ();
  LCAP_GG_98_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree, Genetic_code *ccode);
  ~LCAP_GG_98_model () {delete curr_lin_alg;};
 protected:
  void num_params_model();
  void intialize_parameters (double par[], PARAM_TYPE types[]);
};




class MultiMatrix_JC_model : public JC_model, public MultiMatrix_model
{
 public:
  MultiMatrix_JC_model ();
  MultiMatrix_JC_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree, Genetic_code *ccode, AA_matrix **the_mats);
  ~MultiMatrix_JC_model () {delete curr_lin_alg;};
protected:
  void num_params_model();
  void intialize_parameters (double par[], PARAM_TYPE types[]);
};


class MultiMatrix_K2P_model : public K2P_model, public MultiMatrix_model
{
 public:
  MultiMatrix_K2P_model ();
  MultiMatrix_K2P_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree, Genetic_code *ccode, AA_matrix **the_mats);
  ~MultiMatrix_K2P_model () {delete curr_lin_alg;};
 protected:
  void num_params_model();
  void intialize_parameters (double par[], PARAM_TYPE types[]);
};


class MultiMatrix_HKY_model : public HKY_model, public MultiMatrix_model
{
 public:
  MultiMatrix_HKY_model ();
  MultiMatrix_HKY_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree, Genetic_code *ccode, AA_matrix **the_mats);
  ~MultiMatrix_HKY_model () {delete curr_lin_alg;};
 protected:
  void num_params_model();
  void intialize_parameters (double par[], PARAM_TYPE types[]);
};





#endif





