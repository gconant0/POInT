#ifndef ___LIKELIHOOD_ANNEAL_RECOMBO_H___
#define ___LIKELIHOOD_ANNEAL_RECOMBO_H___

#include "anneal_template.h"
#include "gen_dna_funcs.h"
#include "exchange.h"
#include "genome_tripl_list.h"




//Allows the space point to see the data--hack
extern Clade *global_genomes;
extern WGX_Data *global_homologs;

class Gene_block
{
public:
    Gene_block & operator=(Gene_block &assign_from)  {end=assign_from.end; start=assign_from.start; return(*this);};
    int start, end;
};



class Track_point
{
 public:

    int *homolog_order, move1, move2, pattern, rev, loc, num_blocks, density_size, *location_density, max_breaks;
    BOOL reorg, real_nums, **in_valid_break;
    Clade *the_genomes;
    WGX_Data *the_homologs;
    Exchange *curr_exchange;
    WGX_Tracks *the_tracks;
    Gene_block *block_set;

    Track_point()    {cerr<<"Invalid call to base constructor of Track_point\n";};
    Track_point(Exchange *cexchange);
    Track_point& operator=(Track_point & assign_from); 

    int get_num_blocks()                        {return(num_blocks);};
    int get_block_start(int bl);
    int get_block_end(int bl);
    void set_blocks();
    void set_blocks_v2();
    BOOL check_order();
    int get_score();
    void set_density();
    
    
    ~Track_point ();
protected:
    int cutoff, *pillar_scores, *pillar_lookup, **pillar_neighbor_dirs, search_array_size;
    BOOL score_valid, *used_pillars, **has_ends;
    
    void initialize_pillar_nums();
    void guess_order();
    
   
    void score_neighbors(int current_pillar);
    int get_pillar_search_num(int pillar_num);
    int get_best_pillar();
    int get_compl_pillar(int pillar_lookup_id);
    
 
};

class Track_anneal : public Anneal<Track_point>
{
 public:
    Track_anneal() : Anneal<Track_point>() {};
    Track_anneal(Exchange *cexchange, int nwalkers);
 private:
#ifdef _OPEN_MP_VERSION_
    double *pre_scores;
    Track_point **test_points;
#endif

    virtual int choose_move_start(int walker);
    virtual int choose_far_loc(int walker, int block_start, int block_length);
    double tune_brn;
    void move(int walker);
    void after_move(int walker_num);
    void set_params(Space_point<Track_point> *values) {};
    void initialize_params(Space_point<Track_point> *curr_condition);
    void output();
};


class Track_anneal_density : public Track_anneal
{
public:
    Track_anneal_density();
    Track_anneal_density(Exchange *cexchange, int nwalkers);
                         
    ~Track_anneal_density();
protected:
    
    int choose_move_start(int walker);
    int choose_far_loc(int walker, int block_start, int block_length);
};


#endif









