#include <iostream>
#include <string>
#include "gen_dna_funcs.h"
#include <fstream>
#include "genome_list.h"

using namespace::std;

#ifndef ___GENOME_TRIPL_LIST_H___
#define ___GENOME_TRIPL_LIST_H___





class WGX_Locus {
public:
	WGX_Locus();
	WGX_Locus(int sp_num, int d_level,  string *nms, Genome *genome);
	int get_contig(int index)						{return(contig_index[index]);};
	int get_gene(int index)						{return(gene_index[index]);};
    int get_dupl_level()                        {return(dupl_level);};
    int get_dupl_count()                        {return(dupl_count);};
	Gene * get_gene_obj(int index);				
	BOOL has_duplicate(int index)					{return(has_ith[index]);};
    BOOL has_all_duplicates();
	Genome * get_genome()					{return(the_genome);};
	
    ~WGX_Locus();

protected:
	int species_num, dupl_level, dupl_count, *contig_index, *gene_index;
	string *names;
	BOOL *has_ith;
	Genome *the_genome;

    void alloc_arrays();
	void find_id(string name, int &contig_id, int &gene_id);
	
};


//Note: We assume the input arrays are indexed in the same way as the genome objects in Clade
class Homologs_DX {
public:
	Homologs_DX();
	Homologs_DX(int d_level, string **orthos_n, Clade *genomes);
	WGX_Locus& operator[] (int element);
    int get_dupl_level()   {return(dupl_level);};
	~Homologs_DX();

protected:
    int dupl_level;
	WGX_Locus **the_loci, *null_locus;
	Clade *the_genomes;

};



class WGX_Data {
public:
	WGX_Data();
	WGX_Data(int nhomologs, int d_level, string ***orthos, Clade *genomes);
	Homologs_DX & operator[] (int element);
	int get_num_homologs()					{return(num_homologs);};
    int get_dupl_level()    {return(dupl_level);};
	~WGX_Data();

protected:
	int num_homologs, dupl_level;
	Homologs_DX **the_homologs, *null_homolog;
	Clade *the_genomes;

	int get_species_index(int read_index);
};



class Read_WGX_Data {
public:
	Read_WGX_Data();
	WGX_Data * get_data(string filename, int d_level, int &num_sites, Clade *genomes);
	~Read_WGX_Data();

protected:
	int *species_indexes, sites, num_genomes, dupl_level;
	string *genome_names, ***orthos;
	ifstream infile;
	Clade *the_genomes;
	
	void set_indexes();

};


class Gene_Track_List_DX
{
public:
	Gene_Track_List_DX();
	Gene_Track_List_DX& operator= (Gene_Track_List_DX &assign_from);

	WGX_Locus *my_locus;
    BOOL in_track_line;
	int index_num, to_track_num, dist_to_last, partition_num, mark_num, num_joins, track_pos;
	Gene_Track_List_DX *last, *next;
};

class Tracking_List_DX
{
public:
	Tracking_List_DX();
	Tracking_List_DX(int d_level, Gene_Track_List_DX **ele, Tracking_List_DX *lst, int n);
	~Tracking_List_DX();

    int num, dupl_level;
	Tracking_List_DX *last, *next;
	class Gene_Track_List_DX **element;

};


class Track_List_List_DX
{
	//Note--this class should only be used within the Track_Stack class
public:
	Track_List_List_DX();
	Track_List_List_DX(Gene_Track_List_DX *new_ele, Track_List_List_DX *old_top);
	Gene_Track_List_DX *element;
	Track_List_List_DX *next, *last;
};


class Track_Stack_DX {
public:
	Track_Stack_DX();
	Gene_Track_List_DX * get_bottom();
	Gene_Track_List_DX * get_next();
	void delete_element(Gene_Track_List_DX *element);
	void pop();
	void push(Gene_Track_List_DX *entry);
	int get_size()							{return(size);};
	void reset()							{curr_pos=bottom; pos_num=0;};
	//void restore_orig();
	//void add_back(Gene_Track_List_DX *entry);
	//void clear();
	~Track_Stack_DX();
protected:
	//int size, full_size, orig_size, pos_num;
    int size, pos_num;
	//Track_List_List_DX *top, *bottom, *orig_bottom, *orig_top, *curr_pos, *full_top;
    Track_List_List_DX *top, *bottom,  *curr_pos;
};


class Genome_Track_DX {
public:
	Genome_Track_DX();
    Genome_Track_DX (int id, Genome *genome, WGX_Data *homologs);
	void number_contigs(int *order);
	void number_contigs(int *order, int size);
	void print_tracking();
	void print_tracking(BOOL use_file);
	BOOL diagnose_list();
	Gene_Track_List_DX* get_gene_track(int locus_num, int track_num);
	BOOL has_back_link(int locus_num, int track_num);
    BOOL has_forward_link(int locus_num, int track_num);
	int get_dist_to_last(int locus_num, int track_num);
	int joins_after(int n);
	BOOL all_poss_joins_after(int n);
	BOOL has_tracked_full_break(int n );
	int count_num_breaks();
    int count_num_full_breaks();
	int count_num_breaks(int size);
	int count_num_list_track_breaks();
	int get_list_position_number(int n);
	void set_null_lasts();
	//Tracking_List_DX * create_list_element(int homolog_num);
	int get_nth_tracking_list_pos(int n);
	int break_and_rejoin(int break_before, int size);
	void reset_list_partition_numbers();
	void reset_partition_numbers();

	//void start_track_list(int n);
	//BOOL has_forward_link(int locus_num, int track_num);
	Genome * get_genome()               {return (the_genome);};
	void reset_internal_track_list()    {last_list_pos=0; partial_track_pos=partial_track_start;};
	~Genome_Track_DX();

protected:
    int pow_N[32],taxa_id, list_len, last_list_pos, reversal_pos, num_old, precut_max_partition, postcut_min_partition, *tracks;
    Gene_Track_List_DX **inferred_tracking, **new_inferred_tracking,
		**broken_lefts, **broken_rights, *new_lasts[500], *new_nexts[500],
		*old_lasts[500], *old_nexts[500],  *dumele, **to_ends;
	Tracking_List_DX *partial_track_start, *partial_track_pos, *partial_track_end,
	 **hold_reversal_list;
	Track_Stack_DX  new_left, new_right, track_starts;
	
	Genome *the_genome;
	WGX_Data *the_homologs;

	void recurse_assemble(Track_Stack_DX *&lefts, Track_Stack_DX *&rights,
		int left_end, int size);
    void recurse_assemble_V2(int left_end, int size);
    
    
	void set_stacks(Track_Stack_DX *&new_stack, int end, int stop_point, BOOL left);
	
	int make_joins(Track_Stack_DX *lefts, Track_Stack_DX *rights, Gene_Track_List_DX *last_left, Gene_Track_List_DX *new_right, BOOL do_link, int &save_index);
	void set_link_counts(Track_Stack_DX *lefts, Track_Stack_DX *rights);
	int find_partition_breaks(Tracking_List_DX *left_end, Tracking_List_DX *start, int num_sections_left, int *left_mins, int *left_maxs);
	int find_partition_breaks(int break_after);
	
	
	Tracking_List_DX * get_list_pos(int pos);
	int assign_adjacency(Gene_Track_List_DX *locus1, Gene_Track_List_DX *locus2);
	int assign_adjacency(Gene_Track_List_DX *locus1, Gene_Track_List_DX *locus2, BOOL do_link);
	int assign_adjacency(Gene_Track_List_DX *locus1, Gene_Track_List_DX *locus2, BOOL do_link, int &save_index);
    void mark_used(int left_pos, int right_pos);
    int num_pass_throughs (int pos);
	void store_tracking();
    //Gene_Track_List_DX ** create_new_track_locus(WGX_Locus *the_locus);
    void add_track_starts();
    void set_track_nums();
    void extend_track(Gene_Track_List_DX *my_pos, Gene_Track_List_DX *next_pos, int track_num);

	BOOL check_list();
};



class WGX_Tracks {
public:
    int *order;
	WGX_Tracks()      {cerr<<"Error: Call to default constructor of class WGD_Tracks\n";};
	WGX_Tracks(WGX_Data *homologs, Clade *genomes);
    WGX_Tracks(WGX_Data *homologs, Clade *genomes, int *new_order);
	int get_homolog_num(int num)				{return(order[num]);};
    int get_num_possible_track_orders()     {return(num_orders);};
    int get_single_genome_num_trackings()   {return(single_genome_num_trackings);};
	void change_order(int *new_order);
	void change_order();
	void update_tracking();
	void print_all_tracks();
	void print_all_tracks(BOOL use_file);
    void print_homolog_file(string outfile);
	int get_num_breaks();
    int get_num_positions_w_breaks();
    int get_num_full_breaks();
	int get_num_list_breaks();
	int get_num_move_locations()				{return(num_move_locations);};
	void get_move_section_start_end(int move_pos, int &start, int &end);
	BOOL tracking_current()						{return(tracking_correct);};
	Genome_Track_DX & operator[] (int index);
	void swap_elements (int index1, int index2);
	void set_all_null_lasts();
	void check_for_fully_connected_positions();
	

	~WGX_Tracks();
protected:
	int num_orders, single_genome_num_trackings, *move_locations, *real_to_move_locs, *move_section_starts, *move_section_ends, num_move_locations, num_pseudo_chroms;
	BOOL tracking_correct, *fully_connected;
	WGX_Data *the_homologs;
	Clade *the_genomes;
	Genome_Track_DX **the_trackings;

	void set_array_for_pos(int pos, int val, int size, int *array);

	
};



#endif
