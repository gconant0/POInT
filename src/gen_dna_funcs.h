//General functions of use with DNA sequence work
//Copyright 1999-2002 Gavin Conant

#ifndef ___GEN_DNA_FUNCS_H___
#define ___GEN_DNA_FUNCS_H___
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <string>
#include <algorithm>



#ifndef LN_2
#define LN_2 0.693147180559945

#endif


#define MIN_FLOAT -1e300


//Provides datatypes and classes that make handling DNA and protein sequences easier.
//Functions here define the mapping between amino acid and and nucleotide characters
//and the integer representation of them used elsewhere (i.e. A->0, C->1 ...)

enum LKMODEL {JC_NUCLEOTIDE, K2P_NUCLEOTIDE, HKY_NUCLEOTIDE, GG_98_NUCLEOTIDE, MG_94_JC, MG_94_K2P, 
	      MG_94_HKY, MG_94_GG_98, C_00_JC, C_00_K2P, C_00_HKY, C_00_GG_98,  LCAP_JC, LCAP_K2P, 
	      LCAP_HKY, LCAP_GG_98, AAMULMAT_JC, AAMULMAT_K2P, AAMULMAT_HKY, 
		  BPS_JC, BPS_K2P, BPS_HKY, DUPL, DUPL_FIX, DUPL_PARALLEL, DUPL_PARALLEL_2_RATE, 
		  DUPL_PARALLEL_FIX, DUPL_PARALLEL_FIX_SUBF,
		  DUPL_SUBF_3_RATE, DUPL_SUBF_ONLY, DUPL_2_RATE_NOSUBF,
		  DUPL_NOSTATE, DUPL_FIX_NOSTATE, DUPL_PARALLEL_NOSTATE, DUPL_PARALLEL_2_RATE_NOSTATE, 
		  DUPL_PARALLEL_FIX_NOSTATE, DUPL_SUBF_ALL_STATES,  DUPL_PARALLEL_FIX_SUBF_NOSTATE, 
		  DUPL_SUBF_ALL_STATES_NOSTATE, DUPL_SUBF_3_RATE_NOSTATE, DUPL_SUBF_ONLY_NOSTATE, 
		  DUPL_SUBF_3_RATE_ALLSTATES_NOSTATE, DUPL_2_RATE_NOSUBF_NOSTATE, 
		  DUPL_2_RATE_NOSUBF_ALLSTATES_NOSTATE, DUPL_SLOW_LOSS_CONV_FIX, 
		  DUPL_SLOW_LOSS_CONV_FIX_NOSTATE, DUPL_ARBITRARY, SNP_EVOL_MODEL, ORTHO_PARSIMONY};
enum PARAM_TYPE {PUR_PYR_SPLIT, BRANCH, A_G_SPLIT, C_T_SPLIT, TRS_TRV, CLASS_SUB, SUB_TYPE,
		 AA_PROP_FAC, SCALING_FAC, THIRD_POS_RATE, THREE_TREE_COMMON, 
		 THREE_TREE_SPLIT, THREE_TREE_COMMON_KA, THREE_TREE_SPLIT_KA, MATRIX_COEFF, RATE_PROB, 
		 DUPL_PARALLEL_RATE, DUPL_FIX_RATE, STRAND_SWITCH_PROB, CON_REL_FIX_RATE, CON_FIX_LOSS_RATE, 
		CON_REL_LOSS_RATE, SNP_GAIN_LOSS_RATIO, SNP_LOSS_SWITCH_RATIO, ARBITRARY_PARAM};
#ifndef TRUE
	enum BOOL {FALSE, TRUE};
#else
#define BOOL bool
#endif

enum DATATYPE {NUCLEIC, PROTEIN, CODON, DUPL_STATUS, SNP_DATA, ARBITRARY};
enum DATAFORMAT {NEXUS, PIR, PHYLIP, FASTA, DUPL_DATA};
enum ARG_TYPE {NONE, TREE, DATA, MODEL, NUM_RATES, CFGFILE, 
               SET_BASEFREQS, SET_TRS_TRV, GEN_CODE,  
               PRE_OPTIM, EXCLUDE_STOP,DFORMAT, DTYPE, START, END, HELP, BAD};
enum AA_PROPERTIES {CHEM_COMP, POLARITY, VOLUME, ISO_ELEC, HYDROPATHY, SCALING};
enum GENETIC_CODE {UNIVERSAL, VERT_MITO, YEAST_MITO, MOLD_MITO, INVERT_MITO, CILIATE_NUC, ECHINO_MITO, MYCOPLASMA};
enum CLOCK_TYPE {KS_CLOCK, KA_CLOCK};
enum DUPL_LOSS_STATES {BOTH_PRESENT, COPY1, COPY2, BOTH_PRESENT_FIXED, BOTH_1_BIAS, BOTH_2_BIAS, 
		GENERIC_SINGLE_COPY, MISSING, LOST, COPY1_OR_BOTH, COPY2_OR_BOTH, BOTH_FIXED_SUBF, COPY1_BIAS, COPY2_BIAS};

enum SNP_STATE {TYPE_A, TYPE_B, HETEROZYGOUS, SNP_ABSENT, DATA_MISSING};

struct Entry {
  int position;
  double value;
};

struct Ties {
  int num;
  double value;
};

class Molecule_Sequence 
//Holds a sequence of amino acids or nucleotides and a name for that sequence.
//= operator allows copying of such sequences
{
 public:
  Molecule_Sequence ();
  Molecule_Sequence (int sequence_length);
  Molecule_Sequence (int sequence_length, DATATYPE tp);
  
  int operator==(Molecule_Sequence  & compareto);
  Molecule_Sequence  & operator=(Molecule_Sequence  & assign_from);
  int operator[](int element);
  virtual BOOL compare_elements(int element, Molecule_Sequence * other_seq, int other_element);
  int Sequence_size() { return(size);};
  char * Sequence_name() {return(name);};
  void Assign_site(int site_num, int assignment);
  void Assign_name(const char *assigned_name)   {strcpy(name, assigned_name);};
  DATATYPE get_datatype()		{return(type);};
 
  virtual ~Molecule_Sequence () {if (sequence !=0) delete[] sequence;};

 protected:
  int size, *sequence;
  char name[700];
  DATATYPE type;
};


class Sequence_dataset
//This class allows the creation of objects that act like Molecule_Sequence
//arrays, but which perform array index checking etc.
{
 public:
  Molecule_Sequence **actual_sequences, *dummy;

  
  Sequence_dataset();
  Sequence_dataset(int num_taxa, int sequence_length);
  Sequence_dataset(int num_taxa, int *lenghts);
  Sequence_dataset(int num_taxa, int sequence_length, DATATYPE tp);
  Sequence_dataset(int num_taxa, int *lenghts, DATATYPE tp);
  Sequence_dataset  & operator=(Sequence_dataset  & assign_from);  
  Molecule_Sequence& operator[] (int element);
  int  Num_sequences()  {return(num_sequences);};
  DATATYPE get_datatype() {return(type);};
  ~Sequence_dataset();

 protected:
  int num_sequences;
  DATATYPE type;
};

class DNA_Ambig_Molecule_Sequence : public Molecule_Sequence {
public:
	DNA_Ambig_Molecule_Sequence ();
	DNA_Ambig_Molecule_Sequence (int sequence_length);
	BOOL compare_elements(int element, Molecule_Sequence * other_seq, int other_element);
    ~DNA_Ambig_Molecule_Sequence () {};

};

class List_Molecule_Sequence
{
 public:
  List_Molecule_Sequence() {current=0; last=0; next=0;};
  Molecule_Sequence *current;
  List_Molecule_Sequence *last, *next;
};


class Amino_acid_group
{
 public:
  Amino_acid_group();
  Amino_acid_group(int ngroups);
  Amino_acid_group(int ngroups, const char *filename);
  Amino_acid_group & operator= (Amino_acid_group  & assign_from);  
  void assign_to_group(int aa, int group);
  int get_group(int aa)  {return(groups[aa]);};
  int get_num_groups()   {return(num_groups);};
  ~Amino_acid_group()  {};
 private:
  int num_groups, groups[20];

};


struct Non_gap_sites {
  int site_num;
  Non_gap_sites *next, *last;
}; 


class Btree_node {
	//Holds a simple binary search tree for use with the Constrain_Param_Lookup class
public:
	Btree_node();
	int internal_node_num, tip_num, left_max;
	Btree_node *parent, *children[2];
};

class Constrain_Param_Lookup {
	//If we have a number of parameters n whose sum is fixed at 1, this class allows to have n-1
	//independantly variable parameters which convert back to the original n parameters
public:
	Constrain_Param_Lookup();
	Constrain_Param_Lookup(int num_p);
	double get_prob_value(int category_num);
	double get_param_value(int param_num);
	void set_param_proportion (int param_num, double param_val);
	int get_num_cat_params()			{return(num_categories-1);};
	~Constrain_Param_Lookup();
protected:
	int num_categories, pow2_masks[32], internal_node_pos, internal_node_count;
	double *cat_proportions;
	Btree_node *search_tree, *root;

	Btree_node * build_tree(int min_id, int max_id);
	Btree_node * get_next_free_node();
	double traverse_tree(int index, Btree_node *curr_node);
	
};


#ifndef ONE_MASK
#define ONE_MASK 1
#endif

class BinaryString {
public:
	BinaryString ();
	BinaryString (int l);
	BinaryString& operator= (BinaryString &assign_from);
	int get_len()         {return(len);};
	int element_n(int n);
	int operator[] (int n);
	void set_element_n(int n, int val);
	int string_int_value();
	~BinaryString();
private:
	int *the_string, len, intval;
	BOOL int_val_correct;
};



extern void initialize_next_list_element (Non_gap_sites *&the_list, BOOL first);
extern char num_to_base(int inbase);
extern char num_to_aa(int inaa);
extern char num_to_dupl_data(int indupl);
extern char num_to_snp_state(int instate);
extern int readchar_to_base(char inbase);
extern int readchar_to_aa(char inaa);
extern int readchar_to_dupl(char indupl);
extern int readsnp_to_snpstate(char instate);
extern int get_num_ambig_states(int inbase);
extern BOOL bases_equal(int inbase1, int inbase2);
extern int get_ambig_state_n(int inbase, int ambig_num);
extern int loss_state_to_dupl(DUPL_LOSS_STATES lossstate);
extern DUPL_LOSS_STATES dupl_to_loss_state (int dupl);
extern DUPL_LOSS_STATES get_observable_dupl_state(DUPL_LOSS_STATES inval);
extern SNP_STATE snp_to_snpstate(int state);
extern int snpstate_to_snp(SNP_STATE inval);
extern BOOL is_base(char inbase);
extern BOOL is_aa(char inaa);
extern BOOL is_dupl_data(char indupl);
extern BOOL is_snpstate(char instate);
extern BOOL base_is_ambig(int base);
extern BOOL valid_ambig_spec(int ambig_base, int test_base);
extern BOOL observable_dupl_state(DUPL_LOSS_STATES the_state);
extern double observed_basefreqs (Sequence_dataset *curr_data, int freq);
extern Sequence_dataset * remove_gaps (Sequence_dataset *orig_data, BOOL protein_seq);
extern double observed_codon_basefreqs(Sequence_dataset *curr_data, int codon_pos, int freq);
extern char to_upper(char inchar);
extern char to_lower(char inchar);
extern int exact_match (char string1[], int stringlen, char word[], int wordlen);
extern int word_match(char string1[], int stringlen, char word[], int wordlen);
extern int word_match(std::string my_string, std::string my_word);
extern int exact_match (Molecule_Sequence *seq, int pos1, Molecule_Sequence *motif, int pos2);
extern int loc_word_match (char string1[], int stringlen, char word[], int wordlen);
extern int string_to_int(const char *instring);
//extern double string_to_float(char *instring);
extern double string_to_float(const char *instring);
extern void double_to_string (char *the_string, int string_len, int precision, double val);
extern void int_to_string (char *the_string, int string_len, int val);
extern char to_ucase(char inchar);
extern void to_ucase(char *instring);
extern DATAFORMAT guess_dataformat (const char *filename, int len);
extern BOOL is_transition (int base_1, int base_2);
extern double integer_power(double probablity, int affects);
extern void enum_permutations(int n, int &num_perm, int *ids, int **&perms);
extern int recurse_factorial(int i);
extern double long float_factorial(int i);
extern double long log_factorial(int i);
extern double long table_prob(int vals[2][2]);
extern int n_choose_k(int n, int k);
extern long double float_n_choose_k(int n, int k);
double logadd(double val1, double val2);
double logsubtract(double val1, double val2);
double calc_pearson_correl (int num_elements, double *val1, double *val2);
double calc_spearman_correl (int num_elements, double *val1, double *val2);
double calc_shannon_entropy(int *counts, int num_items);
double calc_mutual_information (double *indep_freqs, double **joint_freqs, int num_items);
double calc_variance(double *vals, int num_items);
double calc_fishers_exact_p(int vals[2][2], BOOL left);
int rank_cmp(const void * a, const void* b);
int delta (int a, int b);
extern int get_dna_comp(int inbase);
extern Sequence_dataset * complement_dna (Sequence_dataset *orig_data);
extern Sequence_dataset * make_writable_dupl_dataset(Sequence_dataset *orig_data);

#endif




