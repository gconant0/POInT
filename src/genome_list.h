#include <iostream>
#include "gen_dna_funcs.h"
#include <fstream>
#include <string>

using namespace::std;

#ifndef ___GENOME_LIST_H___
#define ___GENOME_LIST_H___


class Gene {
public:
    int my_pillar;
    
    
	Gene();
	Gene(int num, char *new_name);
    Gene(int num,  string new_name);
	Gene& operator= (Gene &assign_from);
	const char* get_name()							{return(name.c_str());};
    string get_name_string()                    {return(name);};
	int get_gene_num()							{return(gene_num);};
    int get_start_pos()                         {return(start_pos);};
    int get_end_pos()                           {return(end_pos);};
    string get_chrom_name()                     {return(chrom_string);};
	Gene * get_neighbor(int index);
	void set_neighbor(Gene * new_neighbor, int index);
    void set_location(string ch_string, int s, int e);
	int get_num_neighbors()						{return(num_neighbors);};
	BOOL keep_gene()							{return(keep);};
	void discard_gene()							{keep=FALSE;};
	~Gene();
	
protected:
	int gene_num, num_neighbors, start_pos, end_pos;
	BOOL keep;
	string name, chrom_string;
	Gene **neighbors;
	
};



class Contig {
public:
    
	Contig();
	Contig(int ngenes, char **gene_names);
    Contig(int ngenes,  string *gene_names);
	Gene& operator[] (int element);
	Contig& operator=(Contig & assign_from);
	int get_num_genes()							{return(num_genes);};

	~Contig();

protected:
	int num_genes;
	Gene **the_genes, *null_gene;

	void assign_neighbors();
};



class List_Contig {
public:
	List_Contig();
	List_Contig(List_Contig *lst, Contig *contig);
	List_Contig *last, *next;
	Contig *the_contig;
};


class Genome {
public:
	Genome();
	Genome(int ncontigs, char *name, List_Contig *start_contigs);
    Genome(int ncontigs, string name, List_Contig *start_contigs);
	Contig& operator[] (int element);
	Genome& operator= (Genome &assign_from);
	int get_num_contigs()						{return(num_contigs);};
	const char* get_name()							{return(genome_name.c_str());};
    string get_name_string()                    {return(genome_name);};
    string get_web_link()                       {return(web_link);};
    void set_web_link(string inlink)               {web_link=inlink;};
    void print_genome(string filename);
	~Genome();
protected:
	int num_contigs;
	string genome_name, web_link;
	Contig **the_contigs, *null_contig;
};


class Clade {
public:
	Clade();
	Clade(int ngenomes, Genome **genomes);
	Genome& operator[] (int element);
	int get_num_genomes()						{return(num_genomes);};
	~Clade();
protected:
	int num_genomes;
	Genome **the_genomes, *null_genome;
};



class List_Names {
public:
	List_Names();
	List_Names(List_Names *lst, char *new_name);
    List_Names(List_Names *lst, string new_name);
	List_Names *last, *next;
	string name;
};



class Read_Genome {
public:
	Read_Genome();
	Genome* get_genome(const char *filename);
    Genome* get_genome(string filename);

protected:
	int num_contigs, contig_num, last_contig_num;
	string *name_list, next_contig_name1;
	ifstream infile;
	List_Names *start_names, *list_names;
	List_Contig *start_contigs, *list_contigs;
	

	Contig* get_contig();
};



class WGD_Locus {
public:
	WGD_Locus();
	WGD_Locus(int sp_num, char *n1, char *n2, Genome *genome);
	int get_contig(int index)						{return(contig_index[index]);};
	int get_gene(int index)						{return(gene_index[index]);};
	Gene * get_gene_obj(int index);				
	BOOL has_duplicate()					{return(has_second);};
	Genome * get_genome()					{return(the_genome);};
	

protected:
	int species_num, contig_index[2], gene_index[2];
	char name1[50], name2[50];
	BOOL has_second;
	Genome *the_genome;


	void find_id(char *name, int &contig_id, int &gene_id);
	
};


//Note: We assume the input arrays are indexed in the same way as the genome objects in Clade
class Homologs {
public:
	Homologs();
	Homologs(char **first_orthos, char **second_orthos, Clade *genomes);
	WGD_Locus& operator[] (int element);
	~Homologs();

protected:
	WGD_Locus **the_loci, *null_locus;
	Clade *the_genomes;

};



class WGD_Data {
public:
	WGD_Data();
	WGD_Data(int nhomologs, char ***first_orthos, char ***second_orthos, Clade *genomes);
	Homologs & operator[] (int element);
	int get_num_homologs()					{return(num_homologs);};
	~WGD_Data();

protected:
	int num_homologs;
	Homologs **the_homologs, *null_homolog;
	Clade *the_genomes;

	int get_species_index(int read_index);
};



class Read_WGD_Data {
public:
	Read_WGD_Data();
	WGD_Data * get_data(char *filename, int &num_sites, Clade *genomes);
	~Read_WGD_Data();

protected:
	int *species_indexes, sites, num_genomes;
	char **genome_names, ***first_orthos, ***second_orthos;
	ifstream infile;
	Clade *the_genomes;
	
	void set_indexes();

};


class Gene_Track_List
{
public:
	Gene_Track_List();
	Gene_Track_List& operator= (Gene_Track_List &assign_from);

	WGD_Locus *my_locus;
	int index_num, to_track_num, dist_to_last, partition_num, mark_num, num_joins;
	Gene_Track_List *last, *next;
};

class Tracking_List 
{
public:
	Tracking_List();
	Tracking_List(Gene_Track_List **ele, Tracking_List *lst, int n);
	~Tracking_List();

	int num;
	Tracking_List *last, *next;
	Gene_Track_List **element;

};


class Track_List_List
{
	//Note--this class should only be used within the Track_Stack class
public:
	Track_List_List();
	Track_List_List(Gene_Track_List *new_ele, Track_List_List *old_top);
	Gene_Track_List *element;
	Track_List_List *next, *last;
};


class Track_Stack {
public:
	Track_Stack();
	Gene_Track_List * get_bottom();
	Gene_Track_List * get_next();
	void delete_element(Gene_Track_List *element);
	void pop();
	void push(Gene_Track_List *entry);
	int get_size()							{return(size);};
	void reset()							{curr_pos=bottom; pos_num=0;};
	void restore_orig();
	void add_back(Gene_Track_List *entry);
	void clear();
	~Track_Stack();
protected:
	int size, full_size, orig_size, pos_num;
	Track_List_List *top, *bottom, *orig_bottom, *orig_top, *curr_pos, *full_top;
};


class Reorg_Pattern {
public:
	Reorg_Pattern();
	int get_num_pieces()				{return(num_pieces);};
	int get_insert_loc()				{return(insert_loc);};
	int	get_piece_num(int n);
	int get_reversal(int n);
	BOOL is_last_join_same(int n)			{return(last_join_same[n]);};
	BOOL is_next_join_same(int n)			{return(next_join_same[n]);};
	BOOL is_valid_w_o_add()					{return(valid_without_add);};
	void set_num_pieces(int n)				{num_pieces = n;};
	void set_insert_loc(int val)			{insert_loc=val;};
	void set_piece_num(int n, int val);
	void set_reversal(int n, int val);
	void set_last_same(int n)				{last_join_same[n]=TRUE;};
	void set_next_same(int n)				{next_join_same[n]=TRUE;};
	void not_valid_wo_add()					{valid_without_add=FALSE;};
protected:
	int num_pieces, piece_num[3], insert_loc;
	BOOL last_join_same[2], next_join_same[2], valid_without_add;
	int reversal[3];
};

class All_Reorg_Patterns {
public:
	All_Reorg_Patterns();
	All_Reorg_Patterns(BOOL no_new);
	Reorg_Pattern & operator[] (int element);
	int get_num_patterns()				{return(num_patterns);};
	~All_Reorg_Patterns();
protected:
	int num_patterns;
	Reorg_Pattern *the_patterns;
};

class Genome_Track {
public:
	Genome_Track();    	Genome_Track (int id, Genome *genome, WGD_Data *homologs);
	void number_contigs(int *order);
	void number_contigs(int *order, int size);
	void print_tracking();
	void print_tracking(BOOL use_file);
	BOOL diagnose_list();
	Gene_Track_List* get_gene_track(int locus_num, int track_num);
	BOOL has_back_link(int locus_num, int track_num);
	int get_dist_to_last(int locus_num, int track_num);
	int joins_after(int n);
	int joins_after_list(int n);
	BOOL all_poss_joins_after(int n);
	BOOL has_double_break(int n );
	BOOL has_tracked_double_break(int n );
	int count_num_breaks();
	int count_num_breaks(int size);
	int count_num_list_track_breaks();
	int get_list_position_number(int n);
	void set_null_lasts();
	Tracking_List * create_list_element(int homolog_num);
	int get_nth_tracking_list_pos(int n);
	void try_position_reorgs(int n, int **scores, BOOL *fully_connected);
	void try_position_inserts(int start_pos, int end_pos, int **scores, BOOL *fully_connected);
	int do_section_insert(int start_pos, int end_pos, int loc, int rev, BOOL do_link);
	int do_reorg(int n, int m, int pattern_num, BOOL keep);
	int break_and_rejoin(int break_before, int size);
	int break_and_rejoin(int break_before);
	int break_and_rejoin(Tracking_List *right_end, Tracking_List *start, Tracking_List *end, 
		int &save_index, int start_index, int num_left_sections, int *left_mins, int *left_maxs);
	int break_and_rejoin(Tracking_List *right_end, Tracking_List *start, Tracking_List *end, 
		int &save_index, int start_index, int num_left_sections, int *left_mins, int *left_maxs, BOOL do_link);
	void reset_list_partition_numbers();
	void reset_partition_numbers();

	void start_track_list(int n);
	int add_to_position_n(int n, Tracking_List *new_locus_start, Tracking_List *new_locus_end, int size, BOOL do_link);
	//BOOL has_forward_link(int locus_num, int track_num);
	Genome * get_genome()               {return (the_genome);};
	void reset_internal_track_list()    {last_list_pos=0; partial_track_pos=partial_track_start;};
	~Genome_Track();

protected:
	int taxa_id, pow2[32], list_len, last_list_pos, reversal_pos, num_old, section_mins[4], section_maxs[4], precut_max_partition, postcut_min_partition;
	Gene_Track_List **inferred_tracking, **new_inferred_tracking, 
		*broken_lefts[6], *broken_rights[6], *new_lasts[500], *new_nexts[500], **reversal_pointers,
		*old_lasts[500], *old_nexts[500],  *dumele;
	Tracking_List *partial_track_start, *partial_track_pos, *partial_track_end, 
	 **hold_reversal_list, *section_starts[3][2], *section_ends[3][2], *check_break_pos[2][2][3];
	Track_Stack section_stacks[3][2][2], new_piece, new_left, new_middle_left, new_middle_right, new_right;
	All_Reorg_Patterns possible_patterns, *reorg_patterns;

	Genome *the_genome;
	WGD_Data *the_homologs;

	void recurse_assemble(Track_Stack *&lefts, Track_Stack *&rights, 
		int left_end, int size);
	
	void set_stacks(Track_Stack *new_stack, int end, int stop_point, BOOL left);
	void set_stacks(Track_Stack *new_stack, Tracking_List *end, Tracking_List *stop_point, BOOL left);
	int make_joins(Track_Stack *lefts, Track_Stack *rights, Gene_Track_List *last_left, Gene_Track_List *new_right, BOOL do_link);
	int make_joins(Track_Stack *lefts, Track_Stack *rights, Gene_Track_List *last_left, Gene_Track_List *new_right, BOOL do_link, int &save_index);
	void set_link_counts(Track_Stack *lefts, Track_Stack *rights);
	int set_section_min_maxs(int pattern, int *revs, int *section_pos, int *pos, int break_num, int break_side,
		int break_sec_num[2][3]);
	int set_section_min_maxs(int break_num, int break_side, int rev, int break_sec_num[2][3]);
	int find_partition_breaks(Tracking_List *left_end, Tracking_List *start, int num_sections_left, int *left_mins, int *left_maxs);
	int find_partition_breaks(int break_after);
	void find_possible_internal_rejoins(int num_lost_breaks, int old_lost_breaks, int break_num, int break_sec_num[2][3], 
		int num_check_breaks[2][3], Tracking_List *start_pos);
	void find_reversal_internal_rejoins(int break_num, int num_check_breaks[2][3], Tracking_List *sectionA_end, Tracking_List *sectionB_start);
	Tracking_List * get_list_pos(int pos);
	void reverse_section(Tracking_List *start, Tracking_List *end, Tracking_List *&new_start, Tracking_List *&new_end);
	void reset_reversals()    {reversal_pos=0;};
	int assign_adjacency(Gene_Track_List *locus1, Gene_Track_List *locus2);
	int assign_adjacency(Gene_Track_List *locus1, Gene_Track_List *locus2, BOOL do_link);
	int assign_adjacency(Gene_Track_List *locus1, Gene_Track_List *locus2, BOOL do_link, int &save_index);
	void store_tracking(int size);
	Gene_Track_List ** create_new_track_locus(WGD_Locus *the_locus);

	BOOL check_list();

};



class WGD_Tracks {
public:
	WGD_Tracks()      {cerr<<"Error: Call to default constructor of class WGD_Tracks\n";};
	WGD_Tracks(WGD_Data *homologs, Clade *genomes);
	int get_homolog_num(int num)				{return(order[num]);};
	void change_order(int *new_order);
	void change_order();
	void update_tracking();
	void print_all_tracks();
	void print_all_tracks(BOOL use_file);
	int get_num_breaks();
	int get_num_list_breaks();
	int get_num_move_locations()				{return(num_move_locations);};
	void get_move_section_start_end(int move_pos, int &start, int &end);
	BOOL tracking_current()						{return(tracking_correct);};
	Genome_Track & operator[] (int index);
	void swap_elements (int index1, int index2);
	void set_all_null_lasts();
	void search_for_best();
	int section_insert(int start_pos, int end_pos, int rev, int loc, BOOL do_link, BOOL real_nums);
	void check_for_fully_connected_positions();
	

	~WGD_Tracks();
protected:
	int *order, *move_locations, *real_to_move_locs, *move_section_starts, *move_section_ends, num_move_locations, num_pseudo_chroms;
	BOOL tracking_correct, *fully_connected;
	WGD_Data *the_homologs;
	Clade *the_genomes;
	Genome_Track **the_trackings;

	void set_array_for_pos(int pos, int val, int size, int *array);

	
};



#endif
