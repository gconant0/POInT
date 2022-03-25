#include "gen_dna_funcs.h"
#include "exchange.h"
#include "lin_alg.h"
#include "tree.h"
#include "write_tree.h"
#include <iostream>
#include <string>
#include <fstream>
#include <map>

#ifndef ___PHYLO_MODEL_MATRIX_H____
#define ___PHYLO_MODEL_MATRIX_H____



enum Numerical_Param_Type {MINUS_INF_TO_INF, ZERO_TO_INF, ONE_TO_INF, ONE_PLUS_TO_INF, ZERO_TO_ONE};


enum Entry_Type {ZERO, ONE, PARAM};
class TransProb_Matrix {
public:
    TransProb_Matrix () {transprobs=0; num_states=0; num_rates=0; condprobs=0; std::cerr<<"Error: call to default constructor of TransProb_Matrix\n";};
    TransProb_Matrix (int n_states);
    TransProb_Matrix (int n_states, int n_rates);
    TransProb_Matrix& operator=(TransProb_Matrix& assign_from);
    double get_transprob (int start, int end);
    double get_transprob (int start, int end, int rate);
    void set_transprob(int start, int end, double val);
    void set_transprob(int start, int end, int rate, double val);
    long double get_cond_prob(int state) {return(condprobs[state]);};
    void set_cond_prob(int state, long double val) {condprobs[state]=val;};
    
    ~TransProb_Matrix();
protected:
    int num_rates, num_states;
    double ***transprobs;
    long double *condprobs;
    
};

class Model_State {
public:
    Model_State();
    Model_State(int s_id, string new_name);
    Model_State(int s_id, string new_name, int s_level, int bin_mask, int tlevels);
    int get_state_id();
    string get_state_name();
    int get_state_level()  {return(state_level);};
    int get_binary_rep()     {return(binary_mask);};
    int num_state_redunds();
    int get_ith_redund_state(int redund_num);
    void assign_redunds(std::map<int,int> redund_hash);
    void make_root_state()   {root_state=TRUE;};
    BOOL is_observed_state()  {return(observed_state);};
    BOOL is_root_state()        {return(root_state);};
    BOOL position_present(int pos)  {return(positions_present[pos]);};
    void initialize_cross_ref(int nc);
    void set_cross_ref(int crnum, Model_State *other_state);
    void set_level_state_id(int i)  {level_state_id=i;};
    int get_level_state_id()        {return(level_state_id);};
    Model_State * get_cross_ref(int num) {return(cross_refs[num]);};
    ~Model_State();
protected:
    int state_num, num_redund, *redund_assigns, num_redund_from, state_level, num_cross_refs, binary_mask, level_state_id, total_levels;
    string name;
    BOOL observed_state, root_state, *positions_present;
    Model_State **cross_refs;
    
};

class Param_Set {
public:
    Param_Set();
    Param_Set(string name, Numerical_Param_Type type, double init_val, int set_size, BOOL g);
    Param_Set(string name, Numerical_Param_Type type, double* init_vals, int set_size, BOOL g);
    string get_global_param_name();
    string get_param_set_name(int set_id);
    BOOL is_global_parameter()         {return(is_global);};
    Numerical_Param_Type get_type()    {return(param_type);};
    double get_parameter(int set_num);
    double get_scaled_parameter(int my_set);
    int get_real_set_size()         {return(real_set_size);};
    void set_parameter(double scaled_val, int set_id);
    ~Param_Set();
    
protected:
    int real_set_size;
    BOOL is_global;
    string global_name;
    Numerical_Param_Type param_type;
    double *vals;
};

class Phylo_Matrix {
public:
    double **rate_matrix, **receive_matrix;
    
    Phylo_Matrix();
    Phylo_Matrix(Exchange *cexchange, string model_file);
    
    TransProb_Matrix* get_tp_matrix_num(int matrix_num);
    Model_State * get_nth_state(int n)  {return(the_states[n]);};
    Model_State * get_nth_level_ith_state(int level, int position);
    Model_State * get_ith_full_state_level(int level, int position);
    Model_State * get_masked_state(int mask);
    Model_State * get_first_root_state()  {return(first_root_state);};
    Model_State * get_nth_root_state(int n) {return(root_states[n]);};
    double get_param(int param_num);
    double get_scaled_param(int param_num);
    int get_num_params() {return(num_parameters);};
    double find_ut (Branch *taxa)  {return(taxa->get_brnlen());};
    void set_param(int param_num, double scaled_val);
    void initialize_rate_matrix()       {initialize_rate_matrix(0);};
    void initialize_rate_matrix(int set_id);
    void calc_transprobs(Branch *my_branch);
    void calc_all_transprobs(Tree *curr_tree);
    int get_num_states()    {return(num_states);};
    int num_state_redunds(int state_num);
    int get_ith_redund_state(int state_num, int redund_num);
    int get_num_root_states()  {return(num_root_states);};
    BOOL model_is_hierarchical()   {return(hierarchical_model);};
    BOOL model_is_branch_specific()  {if (num_param_sets ==1) return(FALSE); else return(TRUE);};
    BOOL model_is_symmetric() {return(is_symm);};
    void set_force_brnlen_unity()  {brnlen_one=TRUE;}
    int get_num_levels()  {return(num_levels);};
    int get_num_level_states (int level);
    int get_num_full_level_states (int level);
    string get_model_name() {return(model_name);};
    string get_param_name(int param_num);
    void initialize_all_cross_refs(int level, int nc);
    Entry_Type get_trans_type(int state1, int state2);
    int num_trans_params(int state1, int state2);
    Param_Set* get_nth_trans_param(int state1, int state2, int param_num);
    
    ~Phylo_Matrix();
    
protected:
    int num_states, num_parameters, **num_params_per_entry,
        num_root_states, num_levels, mask_array_size, *num_states_by_level, *total_states_by_level,
        num_param_sets, num_params_per_set, *param_set_lookup, *param_set_id;
    //double *parameters;
    
    Entry_Type **matrix_descript;
    //string *parameter_names, model_name;
    string  model_name;
    BOOL rate_matrix_invalid, hierarchical_model, brnlen_one, is_symm;
    //Numerical_Param_Type *param_types;
    Model_State **the_states, ***level_states, ***full_level_states, **masked_states, *first_root_state, **root_states;
    Exchange *curr_exchange;
    TransProb_Matrix **the_matrices;
    Generic_Linear_Algebra *curr_lin_alg;
    Param_Set **the_parameters, ****matrix_params;
    
    
    Numerical_Param_Type get_param_type(string type_name);
    Param_Set * get_param(string param_name);
    void setup_masked_states();
    void initialize_linear_param_list();
};


class Branch_Ex: public Branch
{
public:
    Branch_Ex();
    Branch_Ex(Exchange *cexchange, TransProb_Matrix *new_mat, int b_num);
    Branch_Ex  & operator=(Branch_Ex  & assign_from);
    double get_trpb(int rate, int start, int end);
    void set_trpb(int rate, int start, int end, double trpb);
    long double get_cond_prob (int site);
    void set_cond_prob (int site, long double prob);
    TransProb_Matrix * get_tb_matrix()  {return my_transprobs;};
    
protected:
    TransProb_Matrix * my_transprobs;
};


class Tree_Ex : public Tree
{
public:
    Tree_Ex();
    Tree_Ex(Exchange *cexchange, Phylo_Matrix *my_mat, BOOL rooted);
    Tree_Ex  & operator=(Tree_Ex  & assign_from);
    Branch_Ex * get_nth_branch(int n);
    ~Tree_Ex()  {if (tree_Ex !=0) delete[] tree_Ex;};
protected:
    Branch_Ex **tree_Ex;
    
};

class Read_PAUP_Tree_Ex
{
public:
    Tree_Ex *curr_tree;
    Exchange *curr_exchange;
    
    //Functions
    Branch * initialize_branch ();
    
    Tree_Ex * create_tree_from_file(Exchange *cexchange, Sequence_dataset *curr_data, Phylo_Matrix *curr_matrix, BOOL rooted);
protected:
    char rel, rel1;
    int nest_level;
    
    virtual Branch* get_tip(ifstream &treein);
    virtual Branch* get_interior(ifstream &treein);
    virtual Branch* make_subtree(ifstream &treein);
    Branch* make_base_unrooted(ifstream &treein);
};

class Read_PAUP_Settings_Tree_Ex : public Read_PAUP_Tree_Ex
{
protected:
    Branch* make_subtree(ifstream &treein);
    Branch* get_tip(ifstream &treein);
    Branch* get_interior(ifstream &treein);
};

class Write_Tree_Arb_Model : public Write_Nexus_Tree
{
public:
    void write_tree(string output_file, string prog_name, Phylo_Matrix *cmatrix, Phylo_Matrix *rmatrix, Tree *ctree, Exchange *cexchange);
    void write_descript_string();
protected:
    Phylo_Matrix *the_matrix, *root_matrix;
};
#endif
