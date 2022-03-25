#ifndef ___GEN_CODE_H___
#define ___GEN_CODE_H___

#include <iostream>
#include <fstream>
#include "gen_dna_funcs.h"
#include "exchange.h"

using namespace::std;

//To make my life easier, I have made the base class "Genetic_code" equivilent to the 
//derived Univ_Genetic_code (so I don't have to re-write a bunch of other programs)
class Genetic_code 
{
 public:
  Genetic_code ();
  Genetic_code (Exchange *cexchange, BOOL standard_code, const char *code_file);
  Genetic_code (Exchange *cexchange);


  int get_pos_n (int codon_num, int position);
  int get_codon_num(int *codon);
  Sequence_dataset * make_codon_sequence (Sequence_dataset *nucleotide_seqs);
  Sequence_dataset * make_codon_sequence_arbitrary_size (Sequence_dataset *nucleotide_seqs);
  BOOL is_non_synon(int *codon1, int *codon2);
  BOOL is_stop (int codon_num);
  int one_diff(int codon_num, int position_num, int iteration);
  int multiple_subs(int *codon1, int *codon2);
  void diff_pos(int *codon1, int *codon2, BOOL is_diff[3]);
  double synon_paths(int num_subs, int *codon1, int *codon2);
  double non_synon_paths(int num_subs, int *codon1, int *codon2);
  void count_by_degen(int *codon1, int *codon2, int diff_pos, double &trs_paths, double &trv_paths);
  void Comeron_levels(int *codon1, int *codon2, int diff_pos, double trs_paths[5], double trv_paths[5], double divisor, int degen);
  void use_zero_cys_cc();
  int Li_93_degen(int codon[3], int pos);
  int Comeron_95_degen(int codon[3], int pos);
  double possible_synon_subs(int *codon);
  double possible_non_synon_subs(int *codon);
  double aa_pair_diff(int codon1[3], int codon2[3], AA_PROPERTIES prop);
  double aa_pair_diff(int aa1, int aa2, AA_PROPERTIES prop);
  double aa_pair_diff(int codon1[3], int codon2[3], int prop);
  double aa_pair_diff(int aa1, int aa2, int prop);
  double aa_prop_val(char aa, int prop) {return individual_aa_props(aa, curr_exchange->get_prop_num_n(prop));};
  char return_AA(int *codon);
  char return_AA(int codon_num);
  int get_non_stops()                {return(non_stops);};
  int get_num_stops()                {return(num_stops);}; 
  int get_stop_num(int stop);
  ~Genetic_code();


 protected:
  int num_stops, non_stops, stop_pos[12];
  char gen_code[4][4][4];
  double ***aa_diffs;
  BOOL cys_cc_zero;
  Exchange *curr_exchange;
 
  virtual void make_code();
  void get_external_code(const char *filename);
  void setup_aa_diffs();
  double individual_aa_props(char aa, AA_PROPERTIES prop);
  void locate_stops();
};


class Univ_genetic_code : public Genetic_code 
{
 public:
    Univ_genetic_code() : Genetic_code () {cerr<<"ERROR: Default constructor of Univ_genetic_code\n";};
    Univ_genetic_code(Exchange *cexchange, BOOL standard_code, char *code_file) :
    Genetic_code(cexchange) {};
    Univ_genetic_code(Exchange *cexchange) :
    Genetic_code(cexchange) {};
};

class Vert_mito_genetic_code : public Genetic_code
{
 public:
  Vert_mito_genetic_code() : Genetic_code () {cerr<<"ERROR: Default constructor of Vert_genetic_code\n";};
  Vert_mito_genetic_code(Exchange *cexchange);
 protected:
  void make_code();
}; 

class Yeast_mito_genetic_code : public Genetic_code
{
 public:
  Yeast_mito_genetic_code() : Genetic_code () {cerr<<"ERROR: Default constructor of Yeast_mito_genetic_code\n";};
  Yeast_mito_genetic_code(Exchange *cexchange);
 protected:
  void make_code();
};

class Mold_mito_genetic_code : public Genetic_code
{
 public:
  Mold_mito_genetic_code() : Genetic_code () {cerr<<"ERROR: Default constructor of Mold_mito_genetic_code\n";};
  Mold_mito_genetic_code(Exchange *cexchange);
 protected:
  void make_code();
};

class Invert_mito_genetic_code : public Genetic_code
{
 public:
  Invert_mito_genetic_code() : Genetic_code () {cerr<<"ERROR: Default constructor of Invert_mito_genetic_code\n";};
  Invert_mito_genetic_code(Exchange *cexchange);
 protected:
  void make_code();
};

class Ciliate_nuc_genetic_code : public Genetic_code
{
 public:
  Ciliate_nuc_genetic_code() : Genetic_code () {cerr<<"ERROR: Default constructor of Ciliate_nuc_genetic_code\n";};
  Ciliate_nuc_genetic_code(Exchange *cexchange);
 protected:
  void make_code();
};

class Echino_mito_genetic_code : public Genetic_code
{
 public:
  Echino_mito_genetic_code() : Genetic_code () {cerr<<"ERROR: Default constructor of Echino_mito_genetic_code\n";};
  Echino_mito_genetic_code(Exchange *cexchange);
 protected:
  void make_code();
};

class Mycoplasma_genetic_code : public Genetic_code
{
public:
	Mycoplasma_genetic_code() : Genetic_code() {cerr<<"ERROR: Default constructor of Mycoplasma_genetic_code\n";};
	Mycoplasma_genetic_code(Exchange *cexchange);
protected:
	void make_code();
};


Genetic_code * create_genetic_code(Exchange *curr_exchange);

#endif







