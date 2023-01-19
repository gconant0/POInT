#ifndef ___READ_TREE_H___
#define ___READ_TREE_H___

#include "tree.h"
#include "gen_dna_funcs.h"


struct Groups_list {
  Groups_list *last, *next;
  int aa[20];
};




class Read_Tree
{
 protected:
  Tree *curr_tree;
  Exchange *curr_exchange;

  //Functions
  Branch * initialize_branch ();  
  Branch * add_to_branch(Branch *add_to, Branch *add);
};

class Read_PAUP_Tree : public Read_Tree
{
 public:

  //Functions
  Tree * create_tree_from_file(Exchange *cexchange, Sequence_dataset *curr_data);
  virtual Tree * create_tree_from_file(Exchange *cexchange, Sequence_dataset *curr_data, BOOL rooted);
  virtual void set_num_nonsyn_params()                          {};

 protected:
    
  //Variables
  char rel, rel1;
  int nest_level;

  //Functions
  virtual Branch* get_tip(ifstream &treein);
  virtual Branch* get_interior(ifstream &treein);
  virtual Branch* make_subtree(ifstream &treein);
  Branch* make_base_unrooted(ifstream &treein); 
 };


class Read_PAUP_w_Settings_Tree : public Read_PAUP_Tree
{
public:
	Read_PAUP_w_Settings_Tree();
	Read_PAUP_w_Settings_Tree(Exchange *cexchange);
	~Read_PAUP_w_Settings_Tree();
protected:
	int *evol_pattern;
	Branch* get_tip(ifstream &treein);
	Branch* get_interior(ifstream &treein);
	
};


class Read_PAUP_w_brlen_CI_Tree : public Read_PAUP_Tree
{
public:
	Read_PAUP_w_brlen_CI_Tree();
	Read_PAUP_w_brlen_CI_Tree(Exchange *cexchange);
	~Read_PAUP_w_brlen_CI_Tree();
protected:
	double con_int[2];
	Branch* get_tip(ifstream &treein);
	Branch* get_interior(ifstream &treein);
	Branch* make_subtree(ifstream &treein);
	
};

class Two_Taxa_Tree
{
 public:
  Tree * create_two_taxa_tree(Exchange *cexchange, double dist);
 private:
  Tree *curr_tree;
};

class Three_Taxa_Tree
{
 public:
  Tree * create_three_taxa_tree(Exchange *cexchange, Sequence_dataset *curr_data, double dist[3]);
 private:
  Tree *curr_tree;
};



#endif

