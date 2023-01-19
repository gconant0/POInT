#include <fstream>
#include "gen_dna_funcs.h"
#include "tree.h"
#include "libranlib.h"
#include "exchange.h"
#include "gen_code.h"
#include "maxlike.h"
#include "codon_like.h"
#include "nucleotide_like.h"
#include "phylo_model_matrix.h"


#ifndef ___SIM_DATA_H___
#define ___SIM_DATA_H___



class Simulate_data
{
 public:
  Sequence_dataset *dataset;

  //Functions
  Simulate_data();
  Simulate_data(Tree *ctree, Exchange *cexchange, Like_model *cmodel);
  void make_dataset();
	void make_dataset(double **rate_probs, int num_probs, int *used_rates);

 protected:
  ifstream infile;
  Tree *curr_tree;
  Like_model *curr_model;
  Exchange *curr_exchange;
  
   

  //Functions
  void get_next_site(int site_num, int rate);
  virtual void ascend_tree(int curr_site, int site_num, Branch *curr_branch, int rate);
  virtual int get_root_site()=0;
  virtual void write_data(int site, int site_num, Branch *taxa)=0;
};



class Simulate_Nucleotide : public Simulate_data 
{
 public:
   Simulate_Nucleotide();
   Simulate_Nucleotide(Tree *ctree, Exchange *cexchange, Like_model *cmodel);
 protected:
  void write_data(int site, int site_num, Branch *taxa);
  int get_root_site();
};


class Simulate_Amino_Acid : public Simulate_data
{
 public:
  Simulate_Amino_Acid();
  Simulate_Amino_Acid(Tree *ctree, Exchange *cexchange, Like_model *cmodel);

 protected:
  void write_data(int site, int site_num, Branch *taxa);
  int get_root_site();
};


class Simulate_Codon : public Simulate_data
{
 public:
  Simulate_Codon();
  Simulate_Codon(Tree *ctree, Exchange *cexchange, Like_model *cmodel, Genetic_code *ccode);
 
 protected:
  Genetic_code *curr_code;
  void write_data(int site, int site_num, Branch *taxa);
  int get_root_site();
};


class Simulate_Dupl : public Simulate_data
{
 public:
  Simulate_Dupl();
  Simulate_Dupl(Tree *ctree, Exchange *cexchange, Like_model *cmodel);
 
 protected:
  void write_data(int site, int site_num, Branch *taxa);
  int get_root_site();
};

class Simulate_from_PhyloMatrix : public Simulate_data
{
public:
    Simulate_from_PhyloMatrix ();
    Simulate_from_PhyloMatrix (Tree *ctree, Exchange *cexchange, Like_model *cmodel, Phylo_Matrix *cmatrix);
protected:
    Phylo_Matrix *curr_model_matrix;
    
    void ascend_tree(int curr_site, int site_num, Branch *curr_branch, int rate);
    void write_data(int site, int site_num, Branch *taxa);
    int get_root_site();
};

#endif





