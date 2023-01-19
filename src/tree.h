//Contains data structures and class definitions for Tree operations
//Used in tree.cpp
//Copyright 1999-2002 Gavin Conant


#ifndef ___TREE_H___
#define ___TREE_H___

#ifndef MIN_BRLEN
#define MIN_BRLEN 1.0e-8
#endif



#include <fstream>
#include "gen_dna_funcs.h"
#include "exchange.h"

using namespace::std;

struct EPattern_list {
	EPattern_list *last, *next;
	int *pattern, pattern_branch;
};


class Branch
{
 public:
  Branch();
  Branch(Exchange *cexchange);
  Branch(Exchange *cexchange, int b_num);
  Branch  & operator=(Branch  & assign_from);

  int get_taxa_id()                                 {return(taxa_id);};
  int get_brn_num()                                 {return(brn_num);};
  void set_taxa_id(int id)                          {taxa_id=id;};
  void set_brn_id(int id)                          {brn_num=id;};
  void set_brnlen(double brln)                      {brlen=brln;};
  double get_brnlen()                               {return(brlen);};
  double expect_subs_site()                         {return(expect_numsubs_site);};
  double expect_nsyn_subs_site()                    {return(nsyn_expect_numsubs_site);};
  void set_expect_subs_site(double expect)          {expect_numsubs_site=expect;};
  void set_expect_nsyn_subs_site(double expect)     {nsyn_expect_numsubs_site=expect;};
  double get_brlen_ci(BOOL upper);
  void set_brlen_ci(BOOL upper, double val);
  void set_p_nonsyn_num(int pns_n)                  {pns_num=pns_n;};
  void set_p_inter_group_num(int pitg_n)	        {pitg_num=pitg_n;};
  void set_aa_prop_num(int aanum, AA_PROPERTIES prop);
  void set_matrix_coeff_num(int matrix, int matrix_num);
  void set_syn_ratio(double ratio)                  {syn_ratio=ratio;};
  void set_nonsyn_ratio(double ratio)               {nonsyn_ratio=ratio;};
  void set_nonsyn_pattern(int pat);
  int get_p_nonsyn_num()                             {return(pns_num);};
  int get_p_intergroup_num()						{return(pitg_num);};
  int get_aa_prop_num(AA_PROPERTIES prop);
  int get_matrix_coeff_num(int matrix);
  int get_aa_prop_num(int prop);
  int get_nonsyn_pattern()                          {return(this->*which_nonsyn_patt);};
  double get_syn_ratio()                            {return(syn_ratio);};
  double get_nonsyn_ratio()                         {return(nonsyn_ratio);};
  char* get_name()                                  {return(name);};
  void set_name(const char * new_name)                    {strcpy(name, new_name);};
  void name_this_branch();
	BOOL branch_has_name()                           {return(has_name);};
    BOOL is_fixed_zero()                            {return(zero_fixed);};
    void make_fixed_zero()                            {zero_fixed=TRUE;};
    void set_pruned()                                    {pruned=TRUE;};
  virtual double get_trpb(int rate, int start, int end);
  virtual void set_trpb(int rate, int start, int end,
		double trpb)                        {transprobs[rate][start][end]=trpb;};
  virtual long double get_cond_prob (int site)                   {return(condprobs[site]);};
  virtual void set_cond_prob (int site, long double prob)        {condprobs[site]=prob;};
#ifdef _OPEN_MP_VERSION_
	long double get_cond_prob_locale (int locale, int site)                   {return(condprobs_locale[locale][site]);};
	void set_cond_prob_locale (int locale, int site, long double prob)        {condprobs_locale[locale][site]=prob;};
#endif
  double get_trpb_prime(int start, int end)         {return(transprobs_prime[0][start][end]);};
  double get_trpb_prime(int start, int end, int rate)         {return(transprobs_prime[rate][start][end]);};

  void set_trpb_prime(int start, int end, double trpb_prime)
    {transprobs_prime[0][start][end]=trpb_prime;};
  void set_trpb_prime(int start, int end, int rate, double trpb_prime)
    {transprobs_prime[rate][start][end]=trpb_prime;};


  double get_trpb_double_prime(int start, int end)  {return(transprobs_double_prime[0][start][end]);};
  double get_trpb_double_prime(int start, int end, int rate)  {return(transprobs_double_prime[rate][start][end]);};
  void set_trpb_double_prime(int start, int end, double trpb_double_prime)
    {transprobs_double_prime[0][start][end]=trpb_double_prime;};
  void set_trpb_double_prime(int start, int end, int rate, double trpb_double_prime)
    {transprobs_double_prime[rate][start][end]=trpb_double_prime;};

  void set_save_conprob_index(int val)    {save_conprobs_index=val;};
  int get_save_conprob_index()				{return(save_conprobs_index);};

  double get_branch_basefreq(int base)              {return(basefreqs[base]);};
  void set_branch_basefreq(int base, double freq)   {basefreqs[base]=freq;};
  BOOL is_tip()                                     {return(tip);};
  void set_tip(BOOL set)                            {tip=set;};
  BOOL is_uninitialized()                           {return(uninitialized);};
    BOOL has_extern_trpb()                          {return(extern_trpb);};
    BOOL is_pruned()                                {return(pruned);};
  void initialized()                                {uninitialized=FALSE;};
  Branch * get_parent()                             {return(parent);};
  Branch * get_child(int childnum);
  Branch * get_sibling()                            {return(sibling);};
  void set_parent(Branch *newpar)                   {parent=newpar;};
  void null_parent()                                {parent=0;};
  void set_child(Branch *newchild, int childnum);
  void null_child(int childnum);
  void set_sibling(Branch *newsib)                  {sibling=newsib;};
  void null_sibling()                               {sibling=0;};
    void set_dist_to_root();
	double get_dist_to_root()						{return(dist_to_root);};
    void set_param_set(int n);
    int get_param_set()           {return(param_set_id);};

  ~Branch();
  
 protected:
	 int taxa_id, brn_num, pns_num, pitg_num, aa_prop_num[NUM_AA_PROPS+1], nonsyn_pattern, Branch::*which_nonsyn_patt,
		 *matrix_coeff_nums, save_conprobs_index, param_set_id;
  BOOL tip, uninitialized, has_name, zero_fixed, extern_trpb, pruned;
  char name[700];
  double expect_numsubs_site, nsyn_expect_numsubs_site, syn_ratio, nonsyn_ratio, brlen,brlen_upper, brlen_lower,  ***transprobs, ***transprobs_prime,
    ***transprobs_double_prime, basefreqs[4], dist_to_root;
  long double *condprobs;
#ifdef _OPEN_MP_VERSION_
	long double **condprobs_locale;
#endif
  Branch *parent,  *children[2], *sibling;
  Exchange *curr_exchange;
};


struct Branch_list {
	Branch_list *next;
	Branch *element;
};

void add_to_branch_list(Branch_list *&list, Branch *new_ele);

class Tree
{
 public:
	Branch **tree;
	//Functions
	Tree();
	Tree(Exchange *cexchange);
	Tree(Exchange *cexchange, BOOL rooted);
	Tree(Exchange *cexchange, Branch *new_root); 
	Branch * operator[] (int element);
	Tree & operator= (Tree & assign_from);
	void describe_tree();
	void name_branches()                {find_root()->name_this_branch();};
	virtual Branch * initialize_branch();
	Branch * get_pott_branch_num(Branch *start, int &curr_num, int target_num);
	int count_branches_above(Branch *start);
	void re_root_tree(Branch * root_from);
	Branch * tree_bisect_reconnect(Branch * split_loc, int tree_1_join, int tree_2_join, int base_tree);
	void prune_and_reconnect(Branch *subtree_root, Branch *reconnect_loc, double split_fac);
	void set_root(Branch *start);
	void set_rooted_tree(BOOL rooted);
	void set_constrain_brn(Branch *the_brn)             {constrain_brn=the_brn;};
	Branch * find_root()                                {return(root);};
	Branch * find_null_branch_id()                      {return(find_null_branch());};
	Branch * find_left_tip(Branch * start);
	Branch * get_leftmost_tip();
	Branch * get_constrain_brn()                        {return(constrain_brn);};
	void diganose_tree(Branch * start);
	void prune_from_tree(Branch *prune_branch);
	void add_back_prune(Branch *add_to);
	void remove_tip(Branch * taxa);
	void run_dist_to_root()                            {find_root()->set_dist_to_root();};
	BOOL is_root(Branch * taxa);
	BOOL rooted_tree()                                 {return(is_rooted_tree);};
	void set_as_siblings(Branch *sib1, Branch *sib2);
	void set_as_parent_child(Branch *parent, Branch *child, int child_num);
	int other_child(int child_num);
	int find_max_tree_depth();
	void copy_tree(Tree *intree);
	void set_num_nonsyn_params();
	~Tree();

 protected:
    
  //Variables
  int num_patterns, max_p_nonsyn, *max_evol_pattern;
  EPattern_list *start_patt, *the_patterns;
  Branch *single_branch, *pair_branch1, *pair_branch2, 
    *old_sibling, *old_parent, *root, *current,
    *left_tip, *dummy, *prune_root, *extra_branch, *constrain_brn;
  Branch_list *start_extra_brn, *extra_brn, *start_tips, *tips;
  BOOL is_rooted_tree, tree_is_local_mem, three_taxa_tree, null_root_brlns;
  Exchange *curr_exchange;
  


  //Functions
	void get_tips(Branch *start, Branch *myroot);
	void build_tree(Branch *tree_loc, Branch *myroot, Branch *&child1, Branch *&child2);
  Tree * split_tree(Branch *split_point); 
  Branch * join_trees(Tree *tree1, Tree *tree2, Branch *join1, Branch *join2, int choosetree);
  Branch* invert_clade(Branch *new_parent, Branch *old_parent, Branch *old_sibling);
  Branch* find_null_branch()                       {return(root->get_child(0));};
  void get_num_distinct_evol_patterns();
	void ascend_assign_clade (Branch *other_branch, Branch *my_pos);
};

void merge_trees (Tree *tree1, Tree *tree2, Exchange *curr_exchange1, Exchange *curr_exchange2, Tree *&new_tree, Exchange *&new_exchange);

void double_tree (Tree *tree, Exchange *curr_exchange, Tree *&new_tree, Exchange *&new_exchange, int level);

#endif






