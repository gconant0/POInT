#ifndef ___WRITE_TREE_H___
#define ___WRITE_TREE_H___
#include <iostream>
#include <fstream>
#include <string>
#include "gen_dna_funcs.h"
#include "score_matrix.h"
#include "tree.h"


enum TREE_TYPE {NEXUS_TREE, PHYLIP_TREE};

class Write_Tree
{
 public:
  //Functions
	Write_Tree() {};
	void write_tree(const char *output_file, const char *prog_name,
			  Tree *ctree, Exchange *cexchange);
    void write_tree(string output_file, string prog_name, Tree *ctree, Exchange *cexchange);
    void write_tree_to_string(string &tree_string, string prog_name, Tree *ctree, Exchange *cexchange);
	void write_tree(const char *output_file, const char *prog_name, AA_matrix **the_mats,
			  Tree *ctree, Exchange *cexchange);
	

 protected:
	int dlines;
	char **descript, program_name[100], matrix_file[100], filename[500];
	//ofstream treefile;
    ostream *treefile;
	Exchange *curr_exchange;
	Tree *tree_object;
	AA_matrix **the_matrices;

	virtual void write_descript_string();
	virtual void write_tree_string()=0;
	virtual void decend_clade(Branch *current);
	virtual void get_sibling(Branch *sib);
};



class Write_Nexus_Tree :  public Write_Tree 
{
public:
	Write_Nexus_Tree() : Write_Tree() {};
 protected:
  void write_tree_string();
  void decend_clade(Branch *current);
  void get_sibling(Branch *sib);
};


class Write_GCC_Tree : public Write_Nexus_Tree
{
public:
	Write_GCC_Tree() : Write_Nexus_Tree() {};
protected:
	void decend_clade(Branch *current);
	void get_sibling(Branch *sib);
	void write_branch_params(Branch *current);
};

class Write_Phylip_Tree : public Write_Tree
{
public:
    Write_Phylip_Tree() : Write_Tree() {};
protected:
    void write_descript_string() {};
    void write_tree_string();
    void decend_clade(Branch *current);
    void get_sibling(Branch *sib);
};


#endif









