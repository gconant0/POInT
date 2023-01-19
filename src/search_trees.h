#include <iostream>
#include <math.h>
#include "tree.h"
#include "maxlike.h"
#include "codon_like.h"
#include "other_like.h"
#include "gen_dna_funcs.h"
#include "exchange.h"
#include "powell.h"
#include "write_tree.h"
#include "gen_code.h"
#include "genome_list.h"
#include "genome_tripl_list.h"
#include "genome_ploidy_like.h"
#include "phylo_model_matrix.h"

class Tree_searcher {
public:
	Tree **best_trees;
	Exchange **best_exchanges;
	Dupl_NoState_Base_model *local_model;

	Tree_searcher() {save_num=0; local_model=0; cerr<<"Error: call to default constructor of class Tree_searcher\n";};
	Tree_searcher(Exchange *cexchange);
	Tree_searcher(Exchange *cexchange, int n);
	Tree_searcher(Exchange *cexchange, Genetic_code *ccode);
	Tree_searcher(Exchange *cexchange, Genetic_code *ccode, int n);
	Tree_searcher(Exchange *cexchange, Clade *cgenomes, WGD_Data *chomologs);
    Tree_searcher(Exchange *cexchange, Clade *cgenomes, WGD_Data *chomologs, char *post_file);
    Tree_searcher(Exchange *cexchange, Clade *cgenomes, WGX_Data *chomologs, Phylo_Matrix *cmatrix);
    Tree_searcher(Exchange *cexchange, Clade *cgenomes, WGX_Data *chomologs, Phylo_Matrix *cmatrix, int n);
	virtual void start_search()=0;
	int get_num_saved()                             {return(num_saved);};
	virtual void write_out_trees(int num, char *partial_file_name, const char *prog_name);
	void set_start_end(int s, int e)   {start=s; end=e;};
	void write_post_probs();
	~Tree_searcher();

protected:
	int save_num, num_saved, *best_num, treenum, start, end;
	char post_prob_file[100];
	double *best_scores, **best_params;
	Exchange *eval_exchange, *curr_exchange;
	Tree *eval_tree;
	Like_model *eval_model;
	Genetic_code *curr_code;
	Powell *eval_powell;
	Write_Nexus_Tree *writeout_tree;
	Clade *the_genomes;
	WGD_Data *the_homologs;
    WGX_Data *the_WGX_homologs;
    Phylo_Matrix *the_matrix;

	void allocate();
    virtual void save_tree(int save_num, Tree * source_tree);
    virtual void setup() {};
	void evaluate_tree(Tree *tree_to_run, Exchange *exchange_to_run);

};

class Exhaustive_Tree_searcher : public Tree_searcher
{
public:
	Exhaustive_Tree_searcher() {cerr<<"Error: call to default constructor of class Exhaustive_Tree_searcher\n";};
	Exhaustive_Tree_searcher(Exchange *cexchange);
	Exhaustive_Tree_searcher(Exchange *cexchange, int n);
	Exhaustive_Tree_searcher(Exchange *cexchange, Genetic_code *ccode);
	Exhaustive_Tree_searcher(Exchange *cexchange, Genetic_code *ccode, int n);
	Exhaustive_Tree_searcher(Exchange *cexchange, Clade *cgenomes, WGD_Data *chomologs);
	Exhaustive_Tree_searcher(Exchange *cexchange, Clade *cgenomes, WGD_Data *chomologs, char *post_file);
	virtual void start_search();
	~Exhaustive_Tree_searcher()  {};
protected:

	void assemble_trees(Tree *partial_tree, Exchange *partial_exchange, int num_added);

};

class Exhaustive_Tree_PhyloMat_searcher : public Tree_searcher
{
public:
    Exhaustive_Tree_PhyloMat_searcher(Exchange *cexchange, Clade *cgenomes, WGX_Data *chomologs, Phylo_Matrix *cmatrix);
    Exhaustive_Tree_PhyloMat_searcher(Exchange *cexchange, Clade *cgenomes, WGX_Data *chomologs, Phylo_Matrix *cmatrix, int n);
    virtual void start_search();
    void write_out_trees(int num, char *partial_file_name, const char *prog_name);
    ~Exhaustive_Tree_PhyloMat_searcher()  {};
protected:
    Tree_Ex *full_tree, **best_trees_Ex;
    Ploidy_Like_model *my_model;
    
    void PhyloMat_allocate();
    void assemble_trees(Tree_Ex *partial_tree, Exchange *partial_exchange, int num_added);
    virtual void setup();
    virtual void save_tree(int save_num, Tree * source_tree);
};
