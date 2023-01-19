#ifndef ___SCAFFOLD_WGX_ANNEAL_H___
#define ___SCAFFOLD_WGX_ANNEAL_H___

#include "anneal_template.h"
#include "gen_dna_funcs.h"
#include "exchange.h"


class Gene_block
{
public:
    Gene_block & operator=(Gene_block &assign_from)  {end=assign_from.end; start=assign_from.start; return(*this);};
    int start, end;
};


class Scaffold_Gene
{
public:
    Scaffold_Gene();
    Scaffold_Gene (int num, string new_name, int cnum);
    Scaffold_Gene& operator=(Scaffold_Gene &assign_from);
    
    string get_name()							{return(name);};
    string get_name_string()                    {return(name);};
    int get_gene_num()							{return(gene_num);};
    Scaffold_Gene * get_neighbor(int index);
    void set_neighbor(Scaffold_Gene * new_neighbor, int index);
    void set_gene_num(int n)                    {gene_num=n;}
    int get_num_neighbors()						{return(num_neighbors);};
    BOOL keep_gene()							{return(keep);};
    void discard_gene()							{keep=FALSE;};
    int get_contig_num()                {return(contig_num);};
    int get_num_homologs()              {return(num_homologs);};
    BOOL gene_used()                    {return(used);};
    void set_used(int pnum)                     {pillar_num=pnum; used=TRUE;};
    void set_unused()                   {used=FALSE; pillar_num=-1;};
    int get_pillar()                    {return(pillar_num);};
    void set_num_homologs(int numh)     {num_homologs=numh;}
    void add_tandem(string new_tandem);
    BOOL is_tandem()                    {return(tandem);};
    int get_num_tandems()               {return(num_tandem);};
    string get_nth_tandem(int n)        {return(tandem_names[n]);};
    ~Scaffold_Gene();
    
protected:
    int gene_num, num_neighbors, contig_num, num_tandem, num_homologs, pillar_num;
    string name, *tandem_names;
    BOOL keep, used, tandem;
    Scaffold_Gene **neighbors;
};


BOOL are_neighbors(Scaffold_Gene* gene1, Scaffold_Gene* gene2);

class Scaffold_Genome
{
public:
    Scaffold_Genome();
    Scaffold_Genome(string gname, int ngenes);
    Scaffold_Genome(string gname, int ngenes, Scaffold_Gene **in_genelist);
    Scaffold_Gene& operator[] (int index_num);
    Scaffold_Genome& operator= (Scaffold_Genome &assign_from);
    
    //void use_gene(int gene_num);
    void omit_gene(int gene_num);
    int get_num_genes()             {return(num_genes);};
    string get_genome_name()        {return(genome_name);};
    void reset_neighbors();
    void print_tandems(string filename);
    
    ~Scaffold_Genome();
protected:
    int num_genes;
    string genome_name;
    Scaffold_Gene **gene_list, *null_gene;
};

Scaffold_Genome* read_genome(string genome_file);
Scaffold_Genome* collapse_tandems (string tandem_file, Scaffold_Genome *orig_genome);

class Gene_homologs
{
public:
    Gene_homologs();
    Gene_homologs(Scaffold_Gene *to_me);
    Gene_homologs(int nhomologs, Scaffold_Gene* to_me, Scaffold_Gene** homologs, double *vals);
    
    int get_num_homologs()          {return(num_homologs);};
    
    void add_homolog(Scaffold_Gene *new_homolog, double val);
    Scaffold_Gene* get_me()                  {return(me);};
    Scaffold_Gene* get_nth_homolog(int num)  {return(the_homologs[num]);};
    Scaffold_Gene* get_closest_homolog();
    double get_nth_val(int num)     {return(the_vals[num]);};
    ~Gene_homologs();
    
protected:
    int num_homologs;
    double *the_vals;
    Scaffold_Gene *me, **the_homologs;
};


class Homolog_Set
{
public:
    Homolog_Set();
    Homolog_Set(Scaffold_Genome *agenome, Scaffold_Genome *dgenome);
    
    Homolog_Set& operator= (Homolog_Set &assign_from);
    Gene_homologs& operator[] (int index_num);
    
    int get_num_ances_genes()                   {return(ances_genome->get_num_genes());};
    void fix_genomes()                          {genomes_fixed=TRUE;};
    BOOL are_genomes_fixed()                      {return(genomes_fixed);};
    
    Scaffold_Genome* get_ancestral_genome()      {return(ances_genome);};
    Scaffold_Genome* get_dupl_genome()          {return(dupl_genome);};
    BOOL stored_order()                         {return(has_order);};
    void get_new_order(string order_file);
    int get_order_pos(int orig_pos);
    
    ~Homolog_Set();
protected:
    BOOL genomes_fixed, has_order;
    int *assigned_order;
    Gene_homologs **the_homologs, *null_homolog;
    Scaffold_Genome *ances_genome, *dupl_genome;
};

Homolog_Set* read_homolog_set (string homolog_file, Scaffold_Genome *ances_genome, Scaffold_Genome *dupl_genome, double cutoff, int count_cutoff);


class Pillar
{
public:
    Pillar();
    Pillar(Exchange *cexchange, Scaffold_Gene *out, Gene_homologs *homos, int id);
    int get_WGX_depth()                                                 {return(depth);};
    int get_pillar_id()                                                 {return(pillar_id);};
    Scaffold_Gene* get_outgroup_gene()                                   {return(outgroup);};
    Scaffold_Gene* get_WGX_gene(int index_num)                          {return(WGX_genes[index_num]);};
    
    Scaffold_Gene* get_outgroup_neighbor(int index_num)                 {return(out_neighbors[index_num]);};
    Scaffold_Gene* get_WGX_neighbor(int WGX_index, int neighbor_index)  {return(WGX_neighbors[WGX_index][neighbor_index]);}
    
    void set_WGX_gene(Scaffold_Gene* the_gene, int index_num);
    void set_outgoup_neighbor(Scaffold_Gene* new_neighbor, int index_num);
    void set_WGX_neighbor(Scaffold_Gene* new_neighbor, int WGX_index, int neighbor_index);
    Gene_homologs* get_my_homologs()                {return(my_homologs);};
    ~Pillar();
                                                                         
protected:
    int depth, pillar_id;
    Scaffold_Gene *outgroup, **WGX_genes, **out_neighbors, ***WGX_neighbors;
    Gene_homologs *my_homologs;
    Exchange *curr_exchange;
};


class WGX_scaffold
{
public:
    WGX_scaffold();
    WGX_scaffold(Homolog_Set *new_homologs, Exchange *cexchange);
    
    Pillar& operator[] (int element);
    WGX_scaffold& operator= (WGX_scaffold &assign_from);
    
    void change_pillar_order (int *new_order);
    int get_num_pillars()                       {return(num_pillars);};
    int get_num_blocks()                        {return(num_blocks);};
    int get_block_start(int bl);
    int get_block_end(int bl);
    void set_blocks();
    int get_num_empty_pillars()                 {return(num_empty);};
    int get_nth_empty_pillar_index(int index)   {return(empty_index[index]);};
    
    int get_order(int pt)                       {return(pillar_order[pt]);};
    void update_neighbors();
    void optimize_assigns();
    double get_neighbor_count();
    int get_num_open_pillars()                  {return(num_open_pillars);};
    Pillar* get_nth_open_pillar(int n)          {return(open_pillar_index[n]);};
    Pillar* get_nth_ordered_pillar(int n)        {return(the_pillars[pillar_order[n]]);};
    Homolog_Set* get_homolog_set()              {return(the_homologs);};
    BOOL order_preset()                         {return(have_init_order);};
    
    ~WGX_scaffold();
protected:
    int num_pillars, num_blocks, *pillar_order, num_open_pillars, num_empty, *empty_index,*order_backref;
    BOOL *pillar_gene_fixed, have_init_order;
    Pillar **the_pillars, *null_pillar, **open_pillar_index;
    Homolog_Set *the_homologs;
    Exchange *curr_exchange;
    Gene_block *block_set, *curr_block_data;
    
    
    void check_pillars();
};




class Scaffold_point
{
 public:

    int *pillar_order, *new_order;
    Exchange *curr_exchange;
    WGX_scaffold *my_scaffold;
    
    Scaffold_point()    {cerr<<"Invalid call to base constructor of Track_point\n";};
    Scaffold_point(Exchange *cexchange);
    Scaffold_point& operator=(Scaffold_point & assign_from);
    BOOL check_order();

    ~Scaffold_point ();
protected:
    
};



class Scaffold_anneal : public Anneal<Scaffold_point>
{
 public:
    Scaffold_anneal() : Anneal<Scaffold_point>() {};
    Scaffold_anneal(Exchange *cexchange, int nwalkers) : Anneal<Scaffold_point>(cexchange, nwalkers) {};
 private:
    void move(int walker);
    void after_move(int walker_num);
    void set_params(Space_point<Scaffold_point> *values) {};
    void initialize_params(Space_point<Scaffold_point> *curr_condition);
    void output();
};


#endif









