//Copyright 1999-2002 Gavin Conant


#include "gen_dna_funcs.h"
#include "maxlike.h"

#ifndef ___NUCLEOTIDE_LIKE_H___
#define ___NUCLEOTIDE_LIKE_H___


class  Nucleotide_model : virtual public Like_model
//This class implements the basic features characteristic of a 
//nucleotide-based model of evolution
{
 public:
  void describe_results();
  double root_freq(int site);
  int rate_param_size();
 protected:
  void initialize_arrays();
  void calc_transprobs(Branch *taxa, int rate_num);
  BOOL change_rate(int rate_num, double *rate_info);
  double prob_nonsyn(int start[3], int end[3], Branch *taxa, int rate_num);
};


//Multiple inheritance is now used to create specific nucleotide models derived from
//the generic models and the Nucleotide_model class

class  JC_Nucleotide_model :  public JC_model, public Nucleotide_model
{
 public:
  JC_Nucleotide_model () {};
  JC_Nucleotide_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree)
    {assemble (cexchange, cdata, ctree);};
  void num_params_model() {curr_exchange->set_num_params(0);};
  void calc_first_derivatives(Branch *taxa, double ut);
  void calc_second_derivatives(Branch *taxa, double ut);
  void calc_first_derivatives(Branch *taxa, double ut, int rate);
  void calc_second_derivatives(Branch *taxa, double ut, int rate);
  double get_obs_trs_trv_ratio();
 protected:
  void intialize_parameters (double par[], PARAM_TYPE types[]);
};


class  K2P_Nucleotide_model :  public K2P_model, public Nucleotide_model
{
 public:
  K2P_Nucleotide_model () {};
  K2P_Nucleotide_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree)
    {assemble (cexchange, cdata, ctree);};
  void num_params_model() {curr_exchange->set_num_params(1);};
  void calc_first_derivatives(Branch *taxa, double ut);
  void calc_second_derivatives(Branch *taxa, double ut);
  void calc_first_derivatives(Branch *taxa, double ut, int rate);
  void calc_second_derivatives(Branch *taxa, double ut, int rate);
  double get_obs_trs_trv_ratio();
 protected:
  void intialize_parameters (double par[], PARAM_TYPE types[]);
};



class  HKY_Nucleotide_model :  public HKY_model, public Nucleotide_model
{
 public:
  HKY_Nucleotide_model () {};
  HKY_Nucleotide_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree)
    {assemble (cexchange, cdata, ctree);};
  void num_params_model();
  void calc_first_derivatives(Branch *taxa, double ut);
  void calc_second_derivatives(Branch *taxa, double ut);
   void calc_first_derivatives(Branch *taxa, double ut, int rate);
  void calc_second_derivatives(Branch *taxa, double ut, int rate);
  double get_obs_trs_trv_ratio();
 protected:
  void intialize_parameters (double par[], PARAM_TYPE types[]);
  double PurPyrfreqs(Branch *taxa, int base);
};


class  GG_98_Nucleotide_model :  public GG_98_model, public Nucleotide_model
{
 public:
  GG_98_Nucleotide_model () {};
  GG_98_Nucleotide_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree)
    {assemble (cexchange, cdata, ctree);};
  void num_params_model();
 protected:
  void intialize_parameters (double par[], PARAM_TYPE types[]);
};



class  HKY_Nucleotide_12_3_model :  public HKY_Nucleotide_model
{
 public:
  HKY_Nucleotide_12_3_model () {};
  HKY_Nucleotide_12_3_model (Exchange *cexchange, Sequence_dataset *cdata, Tree *ctree) 
    {assemble (cexchange, cdata, ctree);};
  void num_params_model() {curr_exchange->set_num_params(5);};
  double find_lnL_on_tree_pos_diff(Tree *ctree);
  double get_obs_trs_trv_ratio();
 protected:
  void intialize_parameters (double par[], PARAM_TYPE types[]);
};



#endif





