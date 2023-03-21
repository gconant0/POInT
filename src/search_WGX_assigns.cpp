#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>
#include "stdio.h"
#include "maxlike.h"
#include "genome_ploidy_like.h"
#include "exchange.h"
#include "tree.h"
#include "genome_tripl_list.h"
#include "powell.h"
#include "gen_dna_funcs.h"
#include "write_tree.h"
#include "phylo_model_matrix.h"
#include <string>
#include <sstream>
#include "search_trees.h"

//#define _DO_PLOT_
#ifdef _DO_PLOT_
#include <plot.h>
#endif

#ifndef POInT_version
#define POInT_version "v1.61"
#endif

void parse_args(int argc, char** argv, Exchange *curr_exchange, int &num_genomes, std::string *&genome_files,
                std::string &ortho_file, std::string &post_probs_file, std::string &model_file, std::string &root_model_file, string &cond_probs_file, std::string &synteny_file, BOOL &have_treefile, int &start, int & end, int &save_num, int &diag_size, BOOL &no_opt, BOOL &get_blocks, BOOL &do_tracking, double &block_thresh, BOOL &degen, BOOL &find_perfect, BOOL &draw_frame, std::string &socketid, std::string *&full_genome_files, std::string *&prefixes,  BOOL &print_all_trees, BOOL &single_only, BOOL &brn_CI, BOOL &draw_model, BOOL &have_loc_data);

void optimize_single_arrange(Exchange *curr_exchange, Clade *the_genomes, WGX_Data *the_homologs, Phylo_Matrix *the_matrix, Phylo_Matrix *root_matrix, std::string post_probs_file, std::string cond_probs_file, int diag_size, BOOL no_opt, BOOL get_blocks, BOOL do_tracking, BOOL draw_frame, std::string socketid, double block_thresh, BOOL degen, BOOL find_perfect, std::string synteny_file, std::string *&full_genome_files, std::string *&prefixes,BOOL &print_all_trees, BOOL &single_only, BOOL brn_CI, BOOL draw_model, BOOL have_loc_data);

#ifdef _DO_PLOT_
void print_blocks(Exchange *curr_exchange, Clade *the_genomes, WGX_Data *the_homologs, Phylo_Matrix *the_matrix, Ploidy_Like_model *the_model, Tree *current_tree, int **&block_ends, int **&track_breaks, double block_thresh, BOOL degen);
void find_multi_spp_blocks(Exchange *curr_exchange, Clade *the_genomes, WGX_Data *the_homologs, Phylo_Matrix *the_matrix, Ploidy_Like_model *the_model, Tree *current_tree, double thresh, int depth,
                           int *&multi_block_ends, int &num_blocks, BOOL degen);

int draw_block_diagram( Clade *the_genomes, WGX_Data *the_homologs,  int **block_ends, int **track_breaks, int block_list_size, int **multi_block_ends, int *num_multi_blocks);
extern void draw_tracking(Exchange *curr_exchange, Tree *the_tree, Phylo_Matrix *the_matrix, Ploidy_Like_model *the_model,  Clade *the_genomes, WGX_Data *the_homologs, int diag_size);
extern void generate_frames(Exchange *curr_exchange, Tree *the_tree, Phylo_Matrix *the_matrix, Ploidy_Like_model *the_model,  Clade *the_genomes, WGX_Data *the_homologs, std::string *&full_genome_files, std::string *&prefixes, BOOL have_loc_data);

extern int draw_model_diag(Phylo_Matrix *the_matrix, string plotfile, BOOL bitmap, BOOL IPC_call, std::stringstream *plot_ss, int matrix_num);
#endif
#ifdef POInT_daemon
extern int communicate_plots(Exchange *curr_exchange, Tree *the_tree, Phylo_Matrix *the_matrix, Ploidy_Like_model *the_model,  Clade *the_genomes, WGX_Data *the_homologs, std::string socketname, std::string *&full_genome_files, std::string *&prefixes, BOOL have_loc_data);
#endif



void find_perfect_ohnos( Clade *the_genomes, WGX_Data *the_homologs, Ploidy_Like_model *the_model);


void optimize_all_trees(Exchange *curr_exchange, Clade *the_genomes, WGX_Data *the_homologs, Phylo_Matrix *the_matrix, int start, int end, int save_num, string order_file);

void calc_brn_second_dervs(Ploidy_Like_model *the_model, Tree *the_tree, Exchange *curr_exchange);
#if 0
void run_branch_checks(BOOL two_sides, int num_checks, Ploidy_Like_model *the_model, Tree *the_tree, Exchange *curr_exchange);
#endif
void get_max_prob_pattern(Exchange *curr_exchange, Phylo_Matrix *the_matrix, Ploidy_Like_model *the_model, int site, int &opt_track_id, int *taxa_track_ids, double &prob);

string make_gene_tree(Exchange *curr_exchange,  Tree *the_tree, Phylo_Matrix *the_matrix, TREE_TYPE my_type, WGX_Data *the_homologs, Ploidy_Like_model *the_model,  Clade *the_genomes, int pillar_num);


#if 0
void write_probs(std::string prob_file, int taxa_id, WGX_Data *the_homologs, Exchange *curr_exchange, Ploidy_Like_model *the_model);
#endif

int main(int argc, char *argv[])
{
	int i, num_genomes, num_homologs, start, end, **block_ends, **track_breaks, save_num, diag_size;
    double block_thresh;
	std::string *genome_files, ortho_file, post_probs_file, model_file, root_model_file, cond_probs_file, socketid,	  synteny_file, *full_genome_files, *prefixes;
	BOOL have_treefile, no_opt, get_blocks, do_tracking, degen, find_perfect, draw_frame, print_all_trees, single_only, brn_CI, draw_model, have_loc_data;
	Exchange current_exchange;
	Sequence_dataset *dummy_seqs;
	Genome **list_of_genomes;
	Read_Genome *read_genome;
	Clade *the_genomes;
	WGX_Data *the_homologs;
	Read_WGX_Data read_homologs;
    Phylo_Matrix *curr_model_matrix, *root_matrix=0;
    
	

	if (argc > 3) {
		parse_args(argc, argv, &current_exchange, num_genomes, genome_files, ortho_file, post_probs_file, model_file, root_model_file, cond_probs_file, synteny_file,
			have_treefile, start, end, save_num, diag_size, no_opt, get_blocks, do_tracking, block_thresh, degen, find_perfect, draw_frame, socketid, full_genome_files, prefixes,  print_all_trees, single_only, brn_CI, draw_model, have_loc_data);
        //cout<<"Got arguements"<<flush<<endl;

		list_of_genomes=new Genome * [num_genomes];
       	for(i=0; i<num_genomes; i++)
		{
			read_genome=new Read_Genome();
			list_of_genomes[i]=read_genome->get_genome(genome_files[i]);
			delete read_genome;
        }
        
        
		the_genomes=new Clade(num_genomes, list_of_genomes);
		current_exchange.set_num_taxa(num_genomes);
		current_exchange.set_num_rates(1);
        
		//Clear the list we made, as the Clade object has its own copy
		//for(i=0; i<num_genomes; i++)
		//	delete list_of_genomes[i];
		//delete[] list_of_genomes;
        
        curr_model_matrix=new Phylo_Matrix(&current_exchange, model_file);
        
        if (curr_model_matrix->model_is_symmetric() == TRUE) cout<<"Model is symmetric\n";
        
        if (root_model_file != "NONE") {
            root_matrix=new Phylo_Matrix(&current_exchange, root_model_file);
            //Branch length is meaningless for the root model
            root_matrix->set_force_brnlen_unity();
            cout<<"Read root model from "<<root_model_file<<endl;
        }
        
		the_homologs=read_homologs.get_data(ortho_file, curr_model_matrix->get_num_levels(), num_homologs, the_genomes);
        
       // for(i=0; i<current_exchange.get_num_taxa(); i++)
        //    cout<<"Main Site 9073 taxa "<<i<<" num dupls: "<<(*the_homologs)[9073][i].get_dupl_count()<<endl;

#if 0
        for(i=0; i<the_genomes->get_num_genomes(); i++) {
            for(j=0; j<the_homologs->get_num_homologs(); j++) {
                for(k=0; k<the_homologs->get_dupl_level(); k++) {
                    if ((*the_homologs)[j][i].get_gene_obj(k) !=0)
                        cout<<(*the_homologs)[j][i].get_gene_obj(k)->get_name()<<"\t";
                    else cout<<"NONE\t";
                }
                cout<<endl;
            }
        }
#endif

		dummy_seqs=new Sequence_dataset(the_genomes->get_num_genomes(), 1, ARBITRARY);


		for(i=0; i<the_genomes->get_num_genomes(); i++)
			(*dummy_seqs)[i].Assign_name((*the_genomes)[i].get_name());

		current_exchange.set_dataset(dummy_seqs);
		current_exchange.set_num_taxa(the_genomes->get_num_genomes());
        
	    current_exchange.set_num_sites(the_homologs->get_num_homologs());
        current_exchange.set_model(DUPL_ARBITRARY, curr_model_matrix->get_num_states());
		
		
		
		if (have_treefile == TRUE)
			optimize_single_arrange(&current_exchange, the_genomes, the_homologs, curr_model_matrix, root_matrix, post_probs_file, cond_probs_file, diag_size, no_opt, get_blocks, do_tracking,draw_frame, socketid, block_thresh, degen, find_perfect, synteny_file, full_genome_files, prefixes, print_all_trees, single_only, brn_CI, draw_model, have_loc_data);
		else
			optimize_all_trees(&current_exchange, the_genomes, the_homologs, curr_model_matrix, start, end, save_num, ortho_file);

        
		delete dummy_seqs;
        delete curr_model_matrix;
		delete the_homologs;
		delete the_genomes;
        for(i=0; i<num_genomes; i++)
        	delete list_of_genomes[i];
        delete[] list_of_genomes;
		return(0);

	}
	else
	{
#ifdef _DO_PLOT_
        cerr<<"POInT "<<POInT_version<<"\nUsage: POInT -g:<genome file> -g:<genome file> -o:<ortholog file> -m:<Model file>  (-r:<Root model file>) (-t:treefile)\n"<<"(-p:<posteriortrackprobs file>) (-b:#) (-B:#) (-i:#) (-w) (-c:<conditional probabilities file> (-no_opt) (-s:<start>:<end>) (-zerolengthfixed) (-x:#TreestoSave) (-a:<complete genome file>) (-estimateBrnCI) (-h:save all single copy gene trees) (-H:Save all gene trees) (-q:FullSyntenyGenesFile) (-draw_models)\n";
#else
        cerr<<"POInT "<<POInT_version<<"\nUsage: POInT -g:<genome file> -g:<genome file> -o:<ortholog file> -m:<Model file>  (-r:<Root model file>) (-t:treefile)\n"<<"(-p:<posteriortrackprobs file>) (-c:<conditional probabilities file> (-no_opt) (-s:<start>:<end>) (-zerolengthfixed) (-x:#TreestoSave) (-estimateBrnCI) (-h:save all single copy gene trees) (-H:Save all gene trees) (-q:FullSyntenyGenesFile)\n";
#endif
		return(-1);
	}

}

void parse_args(int argc, char** argv, Exchange *curr_exchange, int &num_genomes, std::string *&genome_files,
                std::string &ortho_file, std::string &post_probs_file, std::string &model_file, std::string &root_model_file, string &cond_probs_file, std::string &synteny_file, BOOL &have_treefile, int &start, int & end,  int &save_num, int &diag_size, BOOL &no_opt, BOOL &get_blocks, BOOL &do_tracking, double &block_thresh, BOOL &degen,  BOOL &find_perfect, BOOL &draw_frame, std::string &socketid, std::string *&full_genome_files, std::string *&prefixes, BOOL &print_all_trees, BOOL &single_only, BOOL &brn_CI, BOOL &draw_model, BOOL &have_loc_data)
{
    int i, j, cnt_genome=0, cnt_full=0, cnt_pre=0, pos;
    std::string treefile, block_string, size_string;
    char dump[30];
	ifstream treein;
    BOOL one_partial_loc=FALSE;
    
    block_thresh=0.9;
 
    full_genome_files=0;
    prefixes=0;
    
    have_loc_data=FALSE;
    degen=FALSE;
	have_treefile=FALSE;
    no_opt=FALSE;
    get_blocks=FALSE;
    do_tracking=FALSE;
    find_perfect=FALSE;
    draw_frame=FALSE;
    print_all_trees= FALSE;
    single_only=FALSE;
    brn_CI=FALSE;
    draw_model=FALSE;
	treefile="out.tre";
    root_model_file="NONE";
    cond_probs_file="NONE";
    synteny_file="NONE";
    save_num=1;
    diag_size=20;
    
	num_genomes=0;
	start=end=-1;
    
    socketid="";

	for (i=1; i<argc; i++)
		if ((argv[i][1] == 'g') || (argv[i][1] == 'G')) num_genomes++;

	genome_files = new std::string [num_genomes];
	

	post_probs_file="NONE";
	

	for(i=1; i<argc; i++) {
        cout<<"Looking at arg "<<i<<" = "<<argv[i]<<flush<<endl;
        
		switch(argv[i][1]) {
		case 'g':
		case 'G':
            genome_files[cnt_genome]=argv[i];
            genome_files[cnt_genome]=genome_files[cnt_genome].substr(3, genome_files[cnt_genome].length()-3);
            //cout<<"Genome "<<cnt_genome<<" is "<<genome_files[cnt_genome]<<endl;
			cnt_genome++;
		break;
#ifdef _DO_PLOT_
        case 'i':
        case 'I':
                do_tracking=TRUE;
                size_string=argv[i];
                size_string=size_string.substr(3, size_string.length()-3);
                diag_size=string_to_int(size_string.c_str());
                cout<<"Tracking diagrams of size "<<diag_size<<" requested\n";
                
                break;
        case 'w':
        case 'W':
                draw_frame=TRUE;
                socketid=argv[i];
                if (socketid.length() >3) {
                    socketid=socketid.substr(3, socketid.length()-3);
                    cout<<"Running POInT as plotting daemon: "<<socketid<<"\n";
                }
                else socketid="";
                cout<<"Interactive tracking windows requested\n";
                
                break;
        case 'd':
        case 'D':
                draw_model=TRUE;
                break;
#endif
        case 'm':
		case 'M':
            model_file=argv[i];
            model_file=model_file.substr(3, model_file.length()-3);
			break;
                
		case 't':
		case 'T':
            treefile=argv[i];
            treefile=treefile.substr(3, treefile.length()-3);
			curr_exchange->set_treefile(treefile.c_str());
			treein.open(treefile.c_str());
			
			if (!treein.fail()) {
				treein>>dump;
				if (strcmp(dump, "#NEXUS") == 0)
					have_treefile=TRUE;
				treein.close();
				treein.clear();
                    }
			break;
                
		case 'o':
		case 'O':
            ortho_file=argv[i];
            ortho_file=ortho_file.substr(3, ortho_file.length()-3);
			break;
		case 'p':
		case 'P':
            post_probs_file=argv[i];
            post_probs_file=post_probs_file.substr(3, post_probs_file.length()-3);
			break;
					break;
		case 'r':
		case 'R':
                root_model_file=argv[i];
                root_model_file=root_model_file.substr(3, root_model_file.length()-3);
                break;
        case 'c':
        case 'C':
                cond_probs_file=argv[i];
                cond_probs_file=cond_probs_file.substr(3, cond_probs_file.length()-3);
                break;
        case 'q':
        case 'Q':
                synteny_file=argv[i];
                synteny_file=synteny_file.substr(3, synteny_file.length()-3);
                break;
        case 'n':
        case 'N':
                no_opt=TRUE;
                break;
#ifdef _DO_PLOT_
        case 'B':
                degen=TRUE;
        case 'b':
                get_blocks=TRUE;
                cout<<"Parental block diagram requested\n";
                block_string=argv[i];
                block_string=block_string.substr(3, block_string.length()-3);
                block_thresh=string_to_float(block_string.c_str());
                break;
#endif
        case 'x':
        case 'X':
                j=3;
                while(j<strlen(argv[i])) {
                    dump[j-3]=argv[i][j];
                    j++;
                }
                save_num=string_to_int(dump);
                break;
        case 's':
        case 'S':
				j=3;
				while(argv[i][j] != ':') {
					dump[j-3]=argv[i][j];
					j++;
				}
				j++;
				start=string_to_int(dump);
				pos=j;
				while(j<strlen(argv[i])) {
					dump[j-pos]=argv[i][j];
					j++;
				}
				end=string_to_int(dump);
				break;
        case 'z':
        case 'Z':
                cout<<"Forcing zero-length input branches to remain at 0\n";
                curr_exchange->fix_zero_brn_lens();
                break;
        case 'F':
        case 'f':
                find_perfect=TRUE;
                break;

        case 'A':
        case 'a':
                if (argv[i][1] == 'a') one_partial_loc=TRUE;
                if (argv[i][1] == 'A') have_loc_data=TRUE;
                if (full_genome_files==0) {
                    full_genome_files=new string [num_genomes];
                    for(j=0; j<num_genomes; j++) full_genome_files[j]="NONE";
                }
                full_genome_files[cnt_full]=argv[i];
                full_genome_files[cnt_full]=full_genome_files[cnt_full].substr(3, full_genome_files[cnt_full].length()-3);
                cnt_full++;
                cout<<"Have genome order data: "<<argv[i]<<endl;
            break;

        case 'k':
        case 'K':
                if (prefixes ==0){
                    prefixes =new string[num_genomes];
                    for(j=0; j<num_genomes; j++) prefixes[j]="NONE";
                }
                prefixes[cnt_pre]=argv[i];
                prefixes[cnt_pre]=prefixes[cnt_pre].substr(3, prefixes[cnt_pre].length()-3);
                cnt_pre++;
            break;
        case 'e':
        case 'E':
                brn_CI=TRUE;
                break;
                
        case 'h':
                single_only=TRUE;
        case 'H':
                print_all_trees=TRUE;
                cout<<"Printing trees with arguement "<<argv[i]<<endl;
                break;
                
		}

    }
    
    if ((have_loc_data==TRUE) &&( one_partial_loc==TRUE)) {
        have_loc_data=FALSE;
        delete full_genome_files;
        full_genome_files=0;
        cerr<<"ERROR: Mismatch in type of provided full genome lookup files"<<endl;
    }
  
}




void optimize_single_arrange(Exchange *curr_exchange, Clade *the_genomes, WGX_Data *the_homologs, Phylo_Matrix *the_matrix, Phylo_Matrix *root_matrix, std::string post_probs_file,
                             std::string cond_probs_file, int diag_size, BOOL no_opt, BOOL get_blocks, BOOL do_tracking, BOOL draw_frame, std::string socketid, double block_thresh, BOOL degen, BOOL find_perfect, std::string synteny_file, std::string *&full_genome_files, std::string *&prefixes, BOOL &print_all_trees, BOOL &single_only, BOOL brn_CI, BOOL draw_model, BOOL have_loc_data)
//This function will optimize the model parameters given a single tree and ancestral arrangement
{
    int i, j,k,**block_ends, **track_breaks, **multi_block_ends, *num_multi_blocks, stop_depth, num_genes, my_track;
    std::string outtreefile, gene_tree_file, prog_name, tree_string, plotfile;
    BOOL bidir_syn, do_pillar, all_single;
	Ploidy_Like_model *the_model;
	Tree *current_tree;
	Read_PAUP_Tree_Ex get_tree;
    Read_PAUP_Settings_Tree_Ex get_settings_tree;
	Powell current_powell;
	Write_Tree_Arb_Model writeout_tree;
    std::ofstream fout;
    

    outtreefile =  std::string (curr_exchange->get_treefile()) + ".out";
    
   if (the_matrix->model_is_branch_specific() ==FALSE)
        current_tree=get_tree.create_tree_from_file(curr_exchange, curr_exchange->get_dataset(), the_matrix, TRUE);
   else {
       cout<<"Reading branch-specific model tree: "<<curr_exchange->get_treefile()<<endl;
       current_tree=get_settings_tree.create_tree_from_file(curr_exchange, curr_exchange->get_dataset(), the_matrix, TRUE);
   }
	 
    if (root_matrix==0)
        the_model=new Ploidy_Like_model(curr_exchange, current_tree, the_genomes, the_homologs, the_matrix);
    else
        the_model = new Ploidy_Like_RootDelta_model(curr_exchange, current_tree, the_genomes, the_homologs, the_matrix, root_matrix);
    
        
    //the_model->the_tracks->print_all_tracks();
    
    the_model->the_tracks->update_tracking();
	the_model->the_tracks->print_all_tracks();
    
#if 1
    if (synteny_file != "NONE") {
        fout.open(synteny_file.c_str());
        
        if (!fout.fail()) {
            for(i=0; i<the_genomes->get_num_genomes(); i++) {
                for(j=0; j<the_homologs->get_num_homologs()-1; j++) {
                    bidir_syn=TRUE;
                    for(k=0; k<the_homologs->get_dupl_level(); k++) {
                            if (((*the_model->the_tracks)[i].has_back_link(j,k) == FALSE) ||
                                ((*the_model->the_tracks)[i].has_back_link(j+1,k) == FALSE))
                                    bidir_syn=FALSE;
                    }
                    
                    if (bidir_syn==TRUE) {
                        for(k=0; k<the_homologs->get_dupl_level(); k++) {
                            if ((*the_homologs)[j][i].get_gene_obj(k) !=0)
                                fout<<(*the_homologs)[j][i].get_gene_obj(k)->get_name()<<endl;
                        }
                    }
                    
                }
            }
            fout.close();
        }
    }
#endif
    
    
    if(find_perfect==TRUE) find_perfect_ohnos( the_genomes, the_homologs,   the_model);
    cout<<"Total breaks in dataset: "<<the_model->the_tracks->get_num_breaks()<<endl;
   
    std::cout<<the_model->find_appropriate_ln_like()<<endl;

  //  the_model->the_tracks->update_tracking();
  //  the_model->the_tracks->print_all_tracks();
   // the_model->print_tracking_probs("test.txt");
   
    std::cout<<"Second lnl: "<<the_model->find_appropriate_ln_like()<<endl;
    
    if (no_opt == FALSE) {
        cout<<"Begining global optimization: there are "<<curr_exchange->get_num_params()<<" parameters\n"<<flush;

        //the_model->recalculate_transprobs();
        curr_exchange->set_saved_lnL(the_model->find_appropriate_ln_like());
        curr_exchange->set_saved_lnL(current_powell.Init_min(the_model, curr_exchange, FALSE));
    }
    else {
        curr_exchange->set_saved_lnL(the_model->find_appropriate_ln_like());
    }

    if (brn_CI==TRUE)   calc_brn_second_dervs(the_model, current_tree, curr_exchange);
    
	
    if (print_all_trees == TRUE) {
        the_model->get_gene_conditional_probs();
        
        for (i=0; i<the_homologs->get_num_homologs(); i++) {
            if (single_only==FALSE) do_pillar=TRUE;
            else {
                all_single=TRUE;
                do_pillar=TRUE;
                for(j=0; j<the_genomes->get_num_genomes(); j++) {
                    num_genes=0;
                    for(k=0; k<the_homologs->get_dupl_level(); k++)
                        if( (*the_model->the_tracks)[j].get_gene_track(i, k)->my_locus != 0) num_genes++;
                    
                    if (num_genes>1) all_single=FALSE;
                }
                if (all_single==FALSE) do_pillar=FALSE;
                
            }
            
            if (do_pillar==TRUE) {
                std::ostringstream ss;
                ss<<"Pillar"<<i<<".tre";
                gene_tree_file=ss.str();
                
                fout.open(gene_tree_file.c_str());
                tree_string=make_gene_tree(curr_exchange, current_tree, the_matrix, NEXUS_TREE, the_homologs, the_model,  the_genomes, i);
                fout<<tree_string;
                fout.close();
            }
        }
    }
    
	for(i=0; i<curr_exchange->get_num_branches(); i++)
		(*current_tree)[i]->set_expect_subs_site(the_model->get_expect_sub_site((*current_tree)[i]->get_brnlen()));

#ifdef _DO_PLOT_
    if (get_blocks==TRUE) {
        print_blocks(curr_exchange, the_genomes, the_homologs, the_matrix, the_model, current_tree, block_ends, track_breaks, block_thresh,degen);
        
        stop_depth=curr_exchange->get_num_taxa()/2;
        multi_block_ends=new int*[stop_depth+1];
        num_multi_blocks=new int [stop_depth+1];
        
        for(i=1; i<=stop_depth; i++)
                find_multi_spp_blocks(curr_exchange, the_genomes, the_homologs, the_matrix, the_model, current_tree, block_thresh, i, multi_block_ends[i], num_multi_blocks[i], degen);
        draw_block_diagram(the_genomes, the_homologs,  block_ends, track_breaks, stop_depth, multi_block_ends, num_multi_blocks);
        
        for(i=0; i<the_genomes->get_num_genomes(); i++) {
            delete[] block_ends[i];
            delete[] track_breaks[i];
        }
        delete[] block_ends;
        delete[] track_breaks;
    }
#endif
	if (post_probs_file != "NONE")
        the_model->print_tracking_probs(post_probs_file);

	if (cond_probs_file != "NONE")
        the_model->print_transpoint_branch_probs(cond_probs_file);
#ifdef _DO_PLOT_
    if (do_tracking == TRUE) {
        cout<<"Drawing tracking\n";
        draw_tracking(curr_exchange, current_tree, the_matrix, the_model,  the_genomes, the_homologs, diag_size);
    }
    if (draw_model == TRUE) {
        //plotfile="ModelPlot.png";
        plotfile="ModelPlot.ps";
        draw_model_diag(the_matrix, plotfile, FALSE, FALSE, 0, 0);
        if (the_matrix->model_is_branch_specific()==TRUE) {
            plotfile="ModelPlot1.ps";
            draw_model_diag(the_matrix, plotfile, FALSE, FALSE, 0, 1);
        }
    }
    
    
    if (draw_frame == TRUE){
#ifdef POInT_daemon
        cout<<"Daemon plot enabled: socket id is "<<socketid<<endl;
        if (socketid != "") {
            communicate_plots(curr_exchange, current_tree, the_matrix, the_model,  the_genomes, the_homologs, socketid, full_genome_files, prefixes, have_loc_data);
        }
        else {
            generate_frames(curr_exchange, current_tree, the_matrix, the_model,  the_genomes, the_homologs, full_genome_files, prefixes, have_loc_data);
        }
#else
        generate_frames(curr_exchange, current_tree, the_matrix, the_model,  the_genomes, the_homologs, full_genome_files, prefixes, have_loc_data);
#endif
    }
#endif

	
    prog_name = "Tree inferred with POInT ";
    stringstream ss2;
    ss2 << POInT_version;
    prog_name = prog_name + ss2.str();
    prog_name = prog_name + " (fixed tree)";
   
   
    writeout_tree.write_tree(outtreefile, prog_name, the_matrix, root_matrix, current_tree, curr_exchange);
	
	delete the_model;
	delete current_tree;

}

#ifdef _DO_PLOT_
void find_multi_spp_blocks(Exchange *curr_exchange, Clade *the_genomes, WGX_Data *the_homologs, Phylo_Matrix *the_matrix, Ploidy_Like_model *the_model, Tree *current_tree, double thresh,
                           int depth, int *&multi_block_ends, int &num_blocks, BOOL degen)
{
    int taxa, track_id, pillar_id, genome, pillar_start, k,max_track_id, prob_pos, *best_track_states, *my_track_states;
    double last_sum, max_sum, best_prob;
    BOOL *my_omits, *null_omits, match;
    
    //my_omits=new BOOL[curr_exchange->get_num_taxa()];
    
    multi_block_ends=new int [the_homologs->get_num_homologs()];
    
    best_track_states=new int [curr_exchange->get_num_taxa()];
    my_track_states=new int [curr_exchange->get_num_taxa()];
    
    multi_block_ends[0]=0;
    num_blocks=1;
    
    
    
    for(pillar_id=1; pillar_id<the_homologs->get_num_homologs(); pillar_id++) {
        max_track_id=0;
        best_prob=the_model->get_post_prob(pillar_id, 0);
        
        for(track_id=1; track_id<the_model->the_tracks->get_num_possible_track_orders(); track_id++) {
            if (the_model->get_post_prob(pillar_id, track_id) > best_prob) {
                best_prob=the_model->get_post_prob(pillar_id, track_id);
                max_track_id=track_id;
            }
        }
        
        if (degen == TRUE) best_prob*=2.0;
        
        max_sum= the_model->recurse_shared_tracking(max_track_id, pillar_id, 0, depth, my_omits, null_omits);
        if (degen == TRUE) max_sum*=2.0;
        the_model->set_track_states(max_track_id, best_track_states);

        last_sum=0;
        for(track_id=0; track_id<the_model->the_tracks->get_num_possible_track_orders(); track_id++) {
            the_model->set_track_states(track_id, my_track_states);
            match=TRUE;
            
            for(taxa=0; taxa<curr_exchange->get_num_taxa(); taxa++) {
                if ((my_omits[taxa] ==FALSE) && (my_track_states[taxa] != best_track_states[taxa])) match=FALSE;
            }
            
           
            
            if (match == TRUE)
                last_sum+=the_model->get_post_prob(pillar_id-1, track_id);
        }
        if (degen == TRUE) last_sum*=2;
        
        //cout<<"Pillar "<<pillar_id<<": "<<depth<<": best track: "<<max_track_id<<" solo probe "<<best_prob<<" sum : "<<max_sum<< " last: "<<last_sum<<" ";
       // for(taxa=0; taxa<curr_exchange->get_num_taxa(); taxa++) {
         //   if (my_omits[taxa] == TRUE)
           //     cout<<"X";
           // else cout<<"O";
       // }
       // cout<<endl;
        
        //if ( (((max_track_ids[1] == last_max_track_ids[0]) && (max_track_ids[0] == last_max_track_ids[1]))
       //       || ((max_track_ids[0]==last_max_track_ids[0]) && (max_track_ids[1] == last_max_track_ids[1]))) &&
        //     ((max_sum >= thresh) && (last_sum >= thresh)) ) {
        if ((max_sum >= thresh) && (last_sum >= thresh)) {
            multi_block_ends[pillar_id]=0;
            
        }
        else {
             multi_block_ends[pillar_id]=1;
            num_blocks++;
            
        }
        //cout<<"\t"<<num_blocks<<endl;
        
    }
    if (multi_block_ends[the_homologs->get_num_homologs()-1] == 0) {
        multi_block_ends[the_homologs->get_num_homologs()-1]=1;
        num_blocks++;
    }
    
    delete[] best_track_states;
    delete[] my_track_states;
    delete[] my_omits;
    
    
    cout<<"Found a total of "<<num_blocks<<" multi-species blocks at depth "<<depth<<" and prob. threshold: "<<thresh<<"\n";
}
#endif

void find_perfect_ohnos( Clade *the_genomes, WGX_Data *the_homologs, Ploidy_Like_model *the_model)
{
    int num_perfect=0, pillar_id, taxa_id, track_id, *perfect_locs;
    BOOL all_present, all_tracked;
    
    for(pillar_id=1; pillar_id<the_homologs->get_num_homologs()-1; pillar_id++) {
        all_present=TRUE;
        for(taxa_id=0; taxa_id<the_genomes->get_num_genomes(); taxa_id++) {
            for(track_id=0; track_id<the_homologs->get_dupl_level(); track_id++) {
                if ((*the_model->the_tracks)[taxa_id].get_gene_track(pillar_id, track_id)->my_locus == 0) all_present=FALSE;
            }
        }
        
        if (all_present == TRUE) {
             num_perfect++;
            all_tracked=TRUE;
            for(taxa_id=0; taxa_id<the_genomes->get_num_genomes(); taxa_id++) {
                for(track_id=0; track_id<the_homologs->get_dupl_level(); track_id++) {
                    if ((*the_model->the_tracks)[taxa_id].has_back_link(pillar_id, track_id) == FALSE) all_tracked=FALSE;
                    if ((*the_model->the_tracks)[taxa_id].has_back_link(pillar_id+1, track_id) == FALSE) all_tracked=FALSE;
                }
            }
            
            if (all_tracked ==TRUE) {
               
                
                cout<<pillar_id;
                 for(track_id=0; track_id<the_homologs->get_dupl_level(); track_id++) {
                     for(taxa_id=0; taxa_id<the_genomes->get_num_genomes(); taxa_id++) {
                         cout<<"\t"<<(*the_model->the_tracks)[taxa_id].get_gene_track(pillar_id, track_id)->my_locus->get_gene_obj((*the_model->the_tracks)[taxa_id].get_gene_track(pillar_id, track_id)->index_num)->get_name();
                     }
                 }
                cout<<endl;
            
            }
        }
        
    }
        
    cout<<num_perfect<<" pillars were fully duplicated\n";
    
}

#ifdef _DO_PLOT_

void print_blocks(Exchange *curr_exchange, Clade *the_genomes, WGX_Data *the_homologs, Phylo_Matrix *the_matrix, Ploidy_Like_model *the_model, Tree *current_tree, int **&block_ends, int **&track_breaks, double thresh, BOOL degen)
{
    int taxa, pillar_id, genome, pillar_start, k, track_id, num_nonblocked, num_blocks, block_size, num_break,
        *taxa_track_ids, last_track, my_track, max_id, last_max, my_breaks, total_breaks=0;
    double max_prob, *last_probs, *curr_probs;
    
    if (degen == TRUE) {thresh=thresh/2.0;}
    
    last_probs=new double[the_homologs->get_dupl_level()];
    curr_probs=new double[the_homologs->get_dupl_level()];
    
    block_ends=new int *[the_genomes->get_num_genomes()];
    track_breaks=new int *[the_genomes->get_num_genomes()];
    for(taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) {
        block_ends[taxa]=new int [the_homologs->get_num_homologs()];
        track_breaks[taxa]=new int [the_homologs->get_num_homologs()];
    }
    
    the_model->get_gene_conditional_probs();
    taxa_track_ids=new int[the_genomes->get_num_genomes()];
    
    cout<<"Total breaks in dataset: "<<the_model->the_tracks->get_num_breaks()<<endl;
    for(taxa=0; taxa<the_genomes->get_num_genomes(); taxa++)
        cout<<"Taxa "<<taxa<<" has "<<(*the_model->the_tracks)[taxa].count_num_full_breaks()<<" double breaks\n";
    cout<<"Break positions\n";
    cout<<"TAXAID\tPILLAR#\n";
    for(taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) {
        num_nonblocked=0;
        num_blocks=0;
        block_size=0;
        pillar_id=0;
        
        track_breaks[taxa][0]=0;
        block_ends[taxa][0]=0;
        my_breaks=0;
        
        for(my_track=0; my_track<the_homologs->get_dupl_level(); my_track++) last_probs[my_track]=0.0;
        
        for(track_id=0; track_id<the_model->the_tracks->get_num_possible_track_orders(); track_id++) {
            the_model->set_track_states(track_id, taxa_track_ids);
            last_probs[taxa_track_ids[taxa]]+=the_model->get_post_prob(pillar_id, track_id);
        }
        
        last_max=0;
        for(my_track=1; my_track<the_homologs->get_dupl_level(); my_track++) {
            if (last_probs[my_track] >last_probs[last_max]) last_max = my_track;
        }
        
        for(pillar_id=1; pillar_id<the_homologs->get_num_homologs(); pillar_id++) {
            for(my_track=0; my_track<the_homologs->get_dupl_level(); my_track++) curr_probs[my_track]=0.0;
            
            for(track_id=0; track_id<the_model->the_tracks->get_num_possible_track_orders(); track_id++) {
                the_model->set_track_states(track_id, taxa_track_ids);
                curr_probs[taxa_track_ids[taxa]]+=the_model->get_post_prob(pillar_id, track_id);
            }
            
            num_break=0;
            for(k=0; k<the_homologs->get_dupl_level(); k++) {
                if ((*the_model->the_tracks)[taxa].has_back_link(pillar_id+1, k) == FALSE)
                    num_break++;
            }
        
            if (num_break < the_homologs->get_dupl_level())
                track_breaks[taxa][pillar_id]=0;
            else {
                track_breaks[taxa][pillar_id]=1;
                cout<<taxa<<"\t"<<pillar_id<<endl;
                my_breaks++;
            }
            
            max_id=0;
            for(my_track=1; my_track<the_homologs->get_dupl_level(); my_track++) {
                if (curr_probs[my_track] > curr_probs[max_id]) max_id=my_track;
            }
            
            //cout<<taxa<<"\t"<<pillar_id<<"\t"<<max_id<<"\t"<<curr_probs[max_id]<<endl;
            
            if (((max_id != last_max)&& (curr_probs[max_id] > thresh) && (last_probs[last_max] >thresh)) || (curr_probs[max_id] < thresh)) {
           // if (((max_id != last_max)&& (curr_probs[max_id] > 0.7) && (last_probs[last_max] >0.7)) || (curr_probs[max_id] < 0.7)) {
                block_ends[taxa][pillar_id]=1;
            }
            else {
                block_ends[taxa][pillar_id]=0;
            }
            last_max=max_id;
            
            for (my_track=0; my_track<the_homologs->get_dupl_level(); my_track++) last_probs[my_track]=curr_probs[my_track];

        }
        cout<<"Taxa "<<taxa<<" has "<<my_breaks<<" double strand breaks\n";
        total_breaks+=my_breaks;
    }
    cout<<"Total double breaks in dataset: "<<total_breaks<<endl;
    
}


int draw_block_diagram( Clade *the_genomes, WGX_Data *the_homologs,  int **block_ends, int **track_breaks, int block_list_size, int **multi_block_ends, int *num_multi_blocks)
{
    int taxa, taxa2, pillar, thandle, *mblock_starts, *mblock_ends, *mblock_sizes, block_cnt, block_set, last_id, *all_breaks;
    double spacex, spacey, box_len, box_width, border, y_pos, x_pos, line_x, pillar_ratio, label_space;
    char attrib[100], value[100], loc[50];
    BOOL last, all;
    FILE *outfile;
    
    spacex=10000.0;
    spacey=200.0*(block_list_size+the_genomes->get_num_genomes());
    
    border=500.0;
    
    label_space = 0.05*spacex;
    
    box_len = spacex-(2.0*border+label_space);
    box_width = (double)((spacey-2.0*border)/(2.0*(the_genomes->get_num_genomes()+0.5)));
    
    strcpy(attrib, "PAGESIZE");
    strcpy(value, "letter");
    pl_parampl (attrib, value);
    outfile=fopen("outbreaks.ps", "w");
    
    /* create a Postscript Plotter that writes to standard output */
    if ((thandle = pl_newpl ("ps", stdin, outfile, stderr)) < 0)
    {
        fprintf (stderr, "Couldn't create Plotter\n");
        return 1;
    }

    pl_selectpl (thandle);       /* select the Plotter for use */
    
    if (pl_openpl () < 0)       /* open Plotter */
    {
        fprintf (stderr, "Couldn't open Plotter\n");
        return 1;
    }
    pl_erase();
    pl_fspace (0.0, 0.0, spacex, spacey); /* specify user coor system */
    pl_flinewidth (1.0);       /* line thickness in user coordinates */
    pl_pencolorname ("black");    /* path will be drawn in red */
    pl_erase ();                /* erase Plotter's graphics display */
    
    pl_fontname("Helvetica");
    pl_ftextangle (0);         /* text inclination angle (degrees) */
    
    
    y_pos=spacey-border;
    x_pos=label_space+border;
    pl_flinewidth (1.0);
    pl_fline(x_pos, y_pos, x_pos+box_len, y_pos);
    
    for(pillar=0; pillar<the_homologs->get_num_homologs(); pillar++) {
        pillar_ratio=(1.0*pillar)/(1.0*the_homologs->get_num_homologs());
        if ((pillar%500)==0) {
            line_x= pillar_ratio*box_len+x_pos;
            pl_pencolorname ("black");
            
            pl_fline(line_x, y_pos, line_x, y_pos-(0.5*box_width));
            pl_fmove (line_x, y_pos-(0.5*box_width));
            
            pl_ffontsize (12.0);
            
            int_to_string(loc, 49, pillar);
            pl_alabel('c', 'c', loc);
            
        }

    }
    
    pl_ffontsize (16.0);
    
    
    all_breaks=new int[the_homologs->get_num_homologs()];
    for(pillar=0; pillar<the_homologs->get_num_homologs(); pillar++) {
        all_breaks[pillar]=1;
        for(taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) {
            if (track_breaks[taxa][pillar] ==0) all_breaks[pillar]=0;
        }
    
    }
    
    x_pos=label_space+border;
    pl_flinewidth (0.2);
    
    for(pillar=0; pillar<the_homologs->get_num_homologs(); pillar++) {
        pillar_ratio=(1.0*pillar)/(1.0*the_homologs->get_num_homologs());
        if (all_breaks[pillar]==1) {
            line_x= pillar_ratio*box_len+x_pos;
            pl_pencolorname ("orange");
            pl_fline(line_x, border, line_x, spacey-border);
        }
    }
    delete[] all_breaks;
    
    for(taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) {
        y_pos = border+(2*(taxa*box_width));
        pl_pencolorname ("black");
        pl_fmove (border, y_pos);
        pl_alabel('c', 'c', (*the_genomes)[taxa].get_name());
        x_pos=label_space+border;
        
        pl_flinewidth (1.0);       /* line thickness in user coordinates */
        pl_fbox (x_pos, y_pos, x_pos+box_len, y_pos+box_width);
        
        pl_flinewidth (0.2);       /* line thickness in user coordinates */
        
        last=FALSE;
        for(pillar=0; pillar<the_homologs->get_num_homologs(); pillar++) {
            pillar_ratio=(1.0*pillar)/(1.0*the_homologs->get_num_homologs());
            
            if (block_ends[taxa][pillar] ==1) {
                if (last == FALSE) {
                    line_x= pillar_ratio*box_len+x_pos;
                    pl_pencolorname ("red");
                    pl_fline(line_x, y_pos+box_width, line_x, y_pos+(0.7*box_width));
                    last=TRUE;
                }
                else {
                    line_x= pillar_ratio*box_len+x_pos;
                    pl_pencolorname ("gray");
                    pl_fline(line_x, y_pos+box_width, line_x, y_pos+(0.8*box_width));
                }
                
            }
            else last=FALSE;
            
            if (track_breaks[taxa][pillar] ==1) {
                line_x= pillar_ratio*box_len+x_pos;
                pl_pencolorname ("blue");
                pl_fline(line_x, y_pos, line_x, y_pos+(0.3*box_width));
            }
            
            
        }
    }
    
    
    
    
    for(block_set=1; block_set <= block_list_size; block_set++) {
    
        y_pos = border+(2.0*((the_genomes->get_num_genomes()+block_set)*box_width));
        
        
        pl_pencolorname ("black");
        x_pos=border;
        
        pl_fmove (x_pos, y_pos);
        
        pl_ffontsize (12.0);
        
        int_to_string(loc, 49, block_set);
        pl_alabel('c', 'c', loc);
        
        pl_flinewidth (0.0);
        pl_filltype(1);
        
        mblock_starts=new int [num_multi_blocks[block_set]];
        mblock_ends=new int [num_multi_blocks[block_set]];
        mblock_sizes=new int [num_multi_blocks[block_set]];
        
        block_cnt=0;
        mblock_starts[0]=0;
        last_id=1;
        for(pillar=0; pillar<the_homologs->get_num_homologs(); pillar++) {
            if (multi_block_ends[block_set][pillar] ==1) {
                mblock_ends[block_cnt]=pillar-1;
                mblock_sizes[block_cnt]=mblock_ends[block_cnt]-mblock_starts[block_cnt]+1;
                block_cnt++;
                mblock_starts[block_cnt]=pillar;
            }
        }
        
        mblock_ends[block_cnt]=the_homologs->get_num_homologs()-1;
        mblock_sizes[block_cnt]=mblock_ends[block_cnt]-mblock_starts[block_cnt]+1;
        
        if (block_cnt != num_multi_blocks[block_set]-1) cerr<<"ERROR: BLock count error\n";
        
        x_pos=label_space+border;
        
        for(block_cnt=0; block_cnt<num_multi_blocks[block_set]; block_cnt++) {
            //cout<<"Plotting block "<<block_cnt<<" to "<<mblock_starts[block_cnt]<<" to "<<mblock_ends[block_cnt]<<endl;
            cout<<the_genomes->get_num_genomes()-block_set<<"\t"<<block_cnt<<"\t"<<mblock_starts[block_cnt]<<"\t"<<mblock_ends[block_cnt]<<endl;
            if (mblock_sizes[block_cnt] > 5) {
                if (last_id==1) {
                    pl_fillcolorname("blue");
                    last_id=0;
                }
                else {
                    pl_fillcolorname("green");
                    last_id=1;
                }
                
                pl_fbox ((1.0*mblock_starts[block_cnt])/(1.0*the_homologs->get_num_homologs())*box_len+x_pos, y_pos,
                         (1.0*mblock_ends[block_cnt])/(1.0*the_homologs->get_num_homologs())*box_len+x_pos, y_pos+box_width);
            }
            else {
                pl_fillcolorname("gray");
                pl_fbox ((1.0*mblock_starts[block_cnt])/(1.0*the_homologs->get_num_homologs())*box_len+x_pos, y_pos+0.05*box_width,
                         (1.0*mblock_ends[block_cnt])/(1.0*the_homologs->get_num_homologs())*box_len+x_pos, y_pos+0.95*box_width);
            }
            
        }
        delete[] mblock_starts;
        delete[] mblock_ends;
        
    }
    
    
    
    
    if (pl_closepl () < 0)      /* close Plotter */
    {
        fprintf (stderr, "Couldn't close Plotter\n");
        return 1;
    }
    
    pl_selectpl (0);            /* select default Plotter */
    if (pl_deletepl (thandle) < 0) /* delete Plotter we used */
    {
        fprintf (stderr, "Couldn't delete Plotter\n");
        return 1;
    }
    
    return(0);

}
#endif


void optimize_all_trees(Exchange *curr_exchange, Clade *the_genomes, WGX_Data *the_homologs, Phylo_Matrix *the_matrix, int start, int end, int num_save, string order_file)
{
	string prog_name, treefile;
	Exhaustive_Tree_PhyloMat_searcher *the_searcher;
  

	cout<<"Saving a total of "<<num_save<<" tree\n";
    prog_name = "Tree inferred with POInT ";
    stringstream ss2;
    ss2 << POInT_version;
    prog_name = prog_name + ss2.str();
    prog_name = prog_name + " (optimal tree)";
   
    treefile = "searchWGXexhaust_" + order_file;
    
    curr_exchange->set_treefile(treefile.c_str());
    the_searcher=new Exhaustive_Tree_PhyloMat_searcher(curr_exchange, the_genomes, the_homologs, the_matrix, num_save);

	if (start != -1) {
		the_searcher->set_start_end(start, end);
		cout<<"Limiting search to trees "<<start<<" through "<<end<<endl;
	}
	the_searcher->start_search();
	the_searcher->write_out_trees(num_save, curr_exchange->get_treefile(), prog_name.c_str());

	//the_searcher->local_model->the_tracks->update_tracking();
	//the_searcher->local_model->the_tracks->print_all_tracks();




	delete the_searcher;
    
    
}

void calc_brn_second_dervs(Ploidy_Like_model *the_model, Tree *the_tree, Exchange *curr_exchange)
{
    int branch;
    double h_val=0.005, like_double_prime, last_like_double_prime, like, like_right, like_left, orig_brn, curr_brn, stderror;
    
    the_tree->name_branches();
    
    
    cout<<"Branch Length Derivative Data\nBranch#\tBranchName\tBrlen\tlnLdoubleprime\tSTDER\n";
    
    for(branch=0; branch<curr_exchange->get_num_branches(); branch++) {
        h_val=0.0025;
        
        orig_brn=(*the_tree)[branch]->get_brnlen();
        like=the_model->find_appropriate_ln_like();
        
        curr_brn=orig_brn-h_val;
        (*the_tree)[branch]->set_brnlen(curr_brn);
        the_model->calc_transprobs((*the_tree)[branch], 0);
        like_left=the_model->find_appropriate_ln_like();
        
        curr_brn=orig_brn+h_val;
        (*the_tree)[branch]->set_brnlen(curr_brn);
        the_model->calc_transprobs((*the_tree)[branch], 0);
        like_right=the_model->find_appropriate_ln_like();
        
        
        like_double_prime = (like_right - (2.0*like)+ like_left) / (h_val*h_val);
        
        last_like_double_prime=2.0*like_double_prime;
        //like_double_prime=2.0*last_like_double_prime;
        //cout<<setw(14)<<setprecision(10)<<"Initial dervi est: L/R/C"<<like_left<<", "<<like_right<<", "<<like<<" = "<<like_double_prime<<endl;
        while ((fabs(last_like_double_prime-like_double_prime)>0.1) && (like_right <like) && (like_left<like)) {
            h_val*=0.5;
            last_like_double_prime=like_double_prime;
            
            curr_brn=orig_brn-h_val;
            (*the_tree)[branch]->set_brnlen(curr_brn);
            the_model->calc_transprobs((*the_tree)[branch], 0);
            like_left=the_model->find_appropriate_ln_like();
            
            curr_brn=orig_brn+h_val;
            (*the_tree)[branch]->set_brnlen(curr_brn);
            the_model->calc_transprobs((*the_tree)[branch], 0);
            like_right=the_model->find_appropriate_ln_like();
            
            
            like_double_prime = (like_right - (2.0*like)+ like_left) / (h_val*h_val);
            //cout<<setw(14)<<setprecision(10)<<"Dervi est with h= "<<h_val<<" L/R/C"<<like_left<<", "<<like_right<<", "<<like<<" = "" : "<<like_double_prime<<endl;
            
        }
        
        stderror = 1.0/sqrt(-1.0*like_double_prime);
        
        (*the_tree)[branch]->set_brnlen(orig_brn);
        the_model->calc_transprobs((*the_tree)[branch], 0);
        
        cout<<branch<<"\t"<<(*the_tree)[branch]->get_name()<<"\t"<<(*the_tree)[branch]->get_brnlen()<<"\t"<<like_double_prime<<"\t"<<stderror<<endl;
        
        
    }
    
}

#if 0
void run_branch_checks(BOOL two_sides, int num_checks, Ploidy_Like_model *the_model, Tree *the_tree, Exchange *curr_exchange)
{
    int real_check, i, branch;
    double interval, curr_diff, orig_brn, curr_brn;
    
    if (num_checks>20) real_check=20;
    else real_check=num_checks;
    
    if (real_check<3) interval=0.1;
    else {
        if (real_check<5) interval=0.05;
        else interval=0.02;
    }
    the_tree->name_branches();
    
    cout<<"Branch Length CI data\nBranch#\tBranchName\tBrnlen\tlnL\n";
    for(branch=0; branch<curr_exchange->get_num_branches(); branch++) {
        orig_brn=(*the_tree)[branch]->get_brnlen();
        
        cout<<branch<<"\t"<<(*the_tree)[branch]->get_name()<<"\t"<<orig_brn<<"\t"<<setw(12)<<setprecision(8)<<curr_exchange->get_saved_lnL()<<endl;

        curr_brn=orig_brn;
        if (two_sides==TRUE) {
            curr_diff=-1.0*interval;
            for(i=0; i<real_check; i++) {
                curr_brn=orig_brn+(curr_diff*orig_brn);
                (*the_tree)[branch]->set_brnlen(curr_brn);
                the_model->calc_transprobs((*the_tree)[branch], 0);
                cout<<branch<<"\t"<<(*the_tree)[branch]->get_name()<<"\t"<<curr_brn<<"\t"<<setw(12)<<setprecision(8)<<the_model->find_appropriate_ln_like()<<endl;
                curr_diff = curr_diff - interval;
            }
            
        }
        
        curr_diff=1.0*interval;
        for(i=0; i<real_check; i++) {
            curr_brn=orig_brn+(curr_diff*orig_brn);
            (*the_tree)[branch]->set_brnlen(curr_brn);
            the_model->calc_transprobs((*the_tree)[branch], 0);
            cout<<branch<<"\t"<<(*the_tree)[branch]->get_name()<<"\t"<<curr_brn<<"\t"<<setw(12)<<setprecision(8)<<the_model->find_appropriate_ln_like()<<endl;
            curr_diff = curr_diff + interval;
        }
        
        
        (*the_tree)[branch]->set_brnlen(orig_brn);
        the_model->calc_transprobs((*the_tree)[branch], 0);
    }
    
}

#endif
string make_gene_tree(Exchange *curr_exchange,  Tree *the_tree, Phylo_Matrix *the_matrix, TREE_TYPE my_type, WGX_Data *the_homologs, Ploidy_Like_model *the_model,  Clade *the_genomes, int pillar_num)
{
    int i, j,subtree, taxa, level, my_track, max_prob_pattern, *taxa_track_ids, brn, *gene_cnts, num_full, num_zero, num_prune=0, tree_nums[3], taxa_id;
    double max_prob;
    BOOL keep_tree=FALSE, three_two_tree=FALSE, single_tree=FALSE, full_tree=FALSE;
    string tree_string, new_taxa_name;
    Exchange *new_exchange, *pruned_exchange;
    Branch *my_sib, *my_par, *new_root;
    Tree  *new_tree, *pruned_tree;
    Write_Tree *writeout_tree;
    std::map<Branch*, Branch*> tree_map, back_map;
    std::map<Branch*, int> subtree_map;
    taxa_track_ids=new int[the_genomes->get_num_genomes()];
    gene_cnts=new int [the_homologs->get_dupl_level()];
    
    
    get_max_prob_pattern(curr_exchange, the_matrix, the_model, pillar_num, max_prob_pattern, taxa_track_ids, max_prob);
    
    
    for(i=0; i<the_homologs->get_dupl_level(); i++) gene_cnts[i]=0;
    
    for(taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) {
        for (level=0; level<the_homologs->get_dupl_level(); level++) {
            my_track=0;
            while(the_model->tracking_permutes[taxa_track_ids[taxa]][my_track] != level) my_track++;
            
            if ((*the_model->the_tracks)[taxa].get_gene_track(pillar_num, level)->my_locus != 0) gene_cnts[my_track]++;
            else num_prune++;
        }
    }
    
    
    num_full=0;
    num_zero=0;
    
    for(i=0; i<the_homologs->get_dupl_level(); i++) {
        if (gene_cnts[i]==the_genomes->get_num_genomes()) num_full++;
        if (gene_cnts[i]==0) num_zero++;
    }
    
    if (num_full == the_homologs->get_dupl_level()) full_tree=TRUE;
    
    if ((num_full==1) && (num_zero == (the_homologs->get_dupl_level()-1))) {
        new_tree = new_tree = new Tree(curr_exchange, TRUE);
        (*new_tree)=(*the_tree);
        new_exchange=curr_exchange;
        keep_tree=TRUE;
        single_tree=TRUE;
    }
    else {
        //cout<<"Have potential tree: nf: "<<num_full<<" nz: "<<num_zero<<" st: "<<single_tree<<" ft: "<<full_tree<<endl;
        if (num_zero == 1)  {
            double_tree (the_tree,  curr_exchange, new_tree, new_exchange, 2);
            three_two_tree=TRUE;
            if (gene_cnts[0]==0) {
                tree_nums[0]=-1;
                tree_nums[1]=0;
                tree_nums[2]=1;
            }
            else {
                tree_nums[0]=0;
                if (gene_cnts[1] ==0) {
                    tree_nums[1]=-1;
                    tree_nums[2]=1;
                }
                else {
                    tree_nums[1]=1;
                    tree_nums[2]=-1;
                }
            }
        }
        else
            double_tree (the_tree,  curr_exchange, new_tree, new_exchange, the_homologs->get_dupl_level());
    }
    
    if (single_tree == FALSE) {
        for(i=0; i<new_exchange->get_num_branches(); i++) {
            if ((*new_tree)[i]->is_tip() ==TRUE) {
                subtree = (*new_tree)[i]->get_taxa_id()/curr_exchange->get_num_taxa();
                //if (three_two_tree==TRUE) subtree=tree_nums[my_track];
                subtree_map[(*new_tree)[i]]=subtree;
                
                //new_taxa_name=(*new_tree)[i]->get_name();
                std::ostringstream ss;
                ss<<(*new_tree)[i]->get_name()<<subtree;
                new_taxa_name=ss.str();
                
                //cout<<"Renaming taxa "<<(*new_tree)[i]->get_name()<<" to "<<new_taxa_name<<endl;
                (*new_tree)[i]->set_name(new_taxa_name.c_str());
            }
        }
    }
    
    
    if ((single_tree==TRUE)  || (full_tree==TRUE) ){
        pruned_tree=new_tree;
        pruned_exchange=new_exchange;
        
        for(taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) {
            for (level=0; level<the_homologs->get_dupl_level(); level++) {
                my_track=0;
                while(the_model->tracking_permutes[taxa_track_ids[taxa]][my_track] != level) my_track++;
                
                if (three_two_tree==TRUE) my_track=tree_nums[my_track];
                
                if (single_tree == FALSE) {
                    std::ostringstream ss;
                    ss<<(*the_genomes)[taxa].get_name_string()<<my_track;
                    new_taxa_name=ss.str();
                }
                else {
                    std::ostringstream ss;
                    ss<<(*the_genomes)[taxa].get_name_string();
                    new_taxa_name=ss.str();
                }
                
                if ((*the_model->the_tracks)[taxa].get_gene_track(pillar_num, level)->my_locus != 0) {
                    brn=0;
                    while((*new_tree)[brn]->get_name() != new_taxa_name) brn++;
                    
                    (*new_tree)[brn]->set_name((*the_model->the_tracks)[taxa].get_gene_track(pillar_num, level)->my_locus->get_gene_obj((*the_model->the_tracks)[taxa].get_gene_track(pillar_num, level)->index_num)->get_name_string().c_str());
                }
            }
        }
        
        
    }
    else {
        pruned_exchange=new Exchange;
        //cout<<"Alloc exchange"<<endl<<flush;
        (*pruned_exchange)=(*new_exchange);
        //cout<<"Alloc/Copy exchange with np= "<<num_prune<<endl<<flush;
        if (three_two_tree==TRUE)
            pruned_exchange->set_num_taxa(new_exchange->get_num_taxa()-(num_prune-the_genomes->get_num_genomes()));
        else
            pruned_exchange->set_num_taxa(new_exchange->get_num_taxa()-num_prune);
        pruned_tree=new Tree(pruned_exchange, TRUE);
        //cout<<"Alloc tree"<<endl<<flush;
        
        
        
        
        for(taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) {
            for (level=0; level<the_homologs->get_dupl_level(); level++) {
                my_track=0;
                while(the_model->tracking_permutes[taxa_track_ids[taxa]][my_track] != level) my_track++;
                
                if (three_two_tree==TRUE) my_track=tree_nums[my_track];
                
                //cout<<"Taxa "<<taxa<<" Level: "<<level<<" My track: "<<my_track<<endl<<flush;
                
                if (my_track != -1) {
                    std::ostringstream ss;
                    ss<<(*the_genomes)[taxa].get_name_string()<<my_track;
                    new_taxa_name=ss.str();
                    
                    if ((*the_model->the_tracks)[taxa].get_gene_track(pillar_num, level)->my_locus == 0) {
                        //cout<<"Taxa "<<new_taxa_name<<" is empty: pruning\n";
                        j=0;
                        while(strcmp((*new_tree)[j]->get_name(), new_taxa_name.c_str()) !=0) j++;
                        
                        new_tree->remove_tip((*new_tree)[j]);
                    }
                    else {
                        brn=0;
                        while((*new_tree)[brn]->get_name() != new_taxa_name) brn++;
                        //cout<<"Renaming "<<new_taxa_name<<" to "<<(*the_model->the_tracks)[taxa].get_gene_track(pillar_num, level)->my_locus->get_gene_obj((*the_model->the_tracks)[taxa].get_gene_track(pillar_num, level)->index_num)->get_name_string()<<endl;
                        (*new_tree)[brn]->set_name((*the_model->the_tracks)[taxa].get_gene_track(pillar_num, level)->my_locus->get_gene_obj((*the_model->the_tracks)[taxa].get_gene_track(pillar_num, level)->index_num)->get_name_string().c_str());
                    }
                }
            }
        }
    }
    
    if ((num_prune>0) && (single_tree==FALSE)) {
        brn=0;
        
        for(i=0; i<new_exchange->get_num_branches(); i++) {
            if ((*new_tree)[i]->is_pruned() ==FALSE) {
                (*pruned_tree)[brn]->initialized();
                (*(*pruned_tree)[brn])=(*(*new_tree)[i]);
                tree_map[(*pruned_tree)[brn]]=(*new_tree)[i];
                back_map[(*new_tree)[i]]=(*pruned_tree)[brn];
                brn++;
            }
        }
        
        taxa_id=0;
        
        for(brn=0; brn <pruned_exchange->get_num_branches(); brn++) {
            if (tree_map[(*pruned_tree)[brn]]->get_parent() !=0) {
                my_par=back_map[tree_map[(*pruned_tree)[brn]]->get_parent()];
                my_sib=back_map[tree_map[(*pruned_tree)[brn]]->get_sibling()];
                
                
                if ((*pruned_tree)[brn]->is_tip() == TRUE) {
                    (*pruned_tree)[brn]->set_taxa_id(taxa_id);
                    taxa_id++;
                }
                //cout<<"Setting "<<brn<<" as "<<(*pruned_tree)[brn]->get_parent()<<" to sib "<<my_sib->get_name()<<" par: "<<my_par->get_name()<<endl<<flush;
                //if (gene_cnts[subtree_map[tree_map[(*pruned_tree)[brn]]]] >1) {
                if (tree_map[(*pruned_tree)[brn]]->get_parent()->get_child(0)==tree_map[(*pruned_tree)[brn]]->get_sibling()) {
                    pruned_tree->set_as_parent_child(my_par, (*pruned_tree)[brn], 1);
                    pruned_tree->set_as_parent_child(my_par, my_sib, 0);
                    pruned_tree->set_as_siblings((*pruned_tree)[brn], my_sib);
                }
                else {
                    pruned_tree->set_as_parent_child(my_par, (*pruned_tree)[brn], 0);
                    pruned_tree->set_as_parent_child(my_par, my_sib, 1);
                    pruned_tree->set_as_siblings((*pruned_tree)[brn], my_sib);
                }
                //}
                
            }
            else {
                new_root=(*pruned_tree)[brn];
                new_root->null_sibling();
                new_root->null_parent();
            }
        }
        pruned_tree->set_root(new_root);
        
    }
    
    
    
    
    
    for(j=0; j<pruned_exchange->get_num_branches(); j++) {
        if (((*pruned_tree)[j] != pruned_tree->find_root()) && ((*pruned_tree)[j] != pruned_tree->find_null_branch_id())) {
            if ((*pruned_tree)[j]->expect_subs_site() ==0) {
                (*pruned_tree)[j]->set_expect_subs_site(0.1);
            }
        }
        
        //cout<<"Branch "<<j<<" name: "<<(*pruned_tree)[j]->get_name()<<" ID "<<(*pruned_tree)[j]<<" P: "<<(*pruned_tree)[j]->get_parent()<<" C1: "<<(*pruned_tree)[j]->get_child(0)<<" C2: "<<(*pruned_tree)[j]->get_child(1);
        //if ((*pruned_tree)[j]->is_tip())
        //    cout<<" Tip: "<<(*pruned_tree)[j]->get_taxa_id()<<endl;
        //else
        //    cout<<" Internal\n";
        
    }
    if (my_type == NEXUS_TREE)
        writeout_tree=new Write_Nexus_Tree();
    else
        writeout_tree=new Write_Phylip_Tree();
    
    
    writeout_tree->write_tree_to_string(tree_string, "POInT Gene Tree", pruned_tree, pruned_exchange);
    
    delete writeout_tree;
    
    if (keep_tree==FALSE) {
        delete new_exchange;
    }
    delete new_tree;
    delete[] gene_cnts;
    delete[] taxa_track_ids;
    
    
    if ((num_prune>0) && (keep_tree=FALSE)) {
        delete pruned_exchange;
        delete pruned_tree;
    }
    
    return(tree_string);
    
}

void get_max_prob_pattern(Exchange *curr_exchange, Phylo_Matrix *the_matrix, Ploidy_Like_model *the_model, int site, int &opt_track_id, int *taxa_track_ids, double &max_prob)
{
    int  track_id, taxa_id, level_id, cnt_breaks, my_track;
    double new_prob;
    BOOL have_double_breaks, use_track;
    string max_genes;
    
    
    
    if (the_model != 0) {
        max_prob =-1.0;
        for(track_id=0; track_id<the_model->the_tracks->get_num_possible_track_orders(); track_id++) {
            if (max_prob ==-1.0) {
                max_prob=the_model->get_post_prob(site, track_id);
                opt_track_id=track_id;
            }
            else {
                if (the_model->get_post_prob(site, track_id) > max_prob) {
                    max_prob=the_model->get_post_prob(site, track_id);
                    opt_track_id=track_id;
                }
            }
            
        }
        
        if (the_model->get_the_homologs()->get_dupl_level()==3) {
            have_double_breaks=FALSE;
            for(taxa_id=0; taxa_id<the_model->get_the_genomes()->get_num_genomes(); taxa_id++) {
                the_model->has_degen[taxa_id]=FALSE;
                cnt_breaks=0;
                for(level_id=0; level_id<the_model->get_the_homologs()->get_dupl_level(); level_id++) {
                    if (((*the_model->the_tracks)[taxa_id].get_gene_track(site, level_id)->last ==0) &&
                        ((*the_model->the_tracks)[taxa_id].get_gene_track(site, level_id)->next ==0)  &&
                        ((*the_model->the_tracks)[taxa_id].get_gene_track(site, level_id)->my_locus == 0)) cnt_breaks++;
                }
                if (cnt_breaks>1) {
                    the_model->has_degen[taxa_id]=TRUE;
                    have_double_breaks=TRUE;
                }
            }
            
            if (have_double_breaks==TRUE) {
                the_model->set_track_states(opt_track_id, taxa_track_ids);
                
                for(taxa_id=0; taxa_id<the_model->get_the_genomes()->get_num_genomes(); taxa_id++) {
                    if (the_model->has_degen[taxa_id] == TRUE) {
                        for(level_id=0; level_id<the_model->get_the_homologs()->get_dupl_level(); level_id++) {
                            if ((*the_model->the_tracks)[taxa_id].get_gene_track(site, level_id)->my_locus != 0) {
                                my_track=0;
                                while(the_model->tracking_permutes[taxa_track_ids[taxa_id]][my_track] != level_id) my_track++;
                                the_model->gene_pos[taxa_id]=my_track;
                            }
                        }
                    }
                    
                }
                
                new_prob=0;
                
                for(track_id=0; track_id<the_model->the_tracks->get_num_possible_track_orders(); track_id++) {
                    the_model->set_track_states(track_id, taxa_track_ids);
                    use_track=TRUE;
                    for(taxa_id=0; taxa_id<the_model->get_the_genomes()->get_num_genomes(); taxa_id++) {
                        if (the_model->has_degen[taxa_id] == TRUE) {
                            for(level_id=0; level_id<the_model->get_the_homologs()->get_dupl_level(); level_id++) {
                                if ((*the_model->the_tracks)[taxa_id].get_gene_track(site, level_id)->my_locus != 0) {
                                    my_track=0;
                                    while(the_model->tracking_permutes[taxa_track_ids[taxa_id]][my_track] != level_id) my_track++;
                                    
                                    if (my_track!=the_model->gene_pos[taxa_id]) use_track=FALSE;
                    
                                }
                            }
                        }
                    }
                    
                    if (use_track == TRUE) new_prob +=the_model->get_post_prob(site, track_id);
                
                }
                the_model->set_track_states(opt_track_id, taxa_track_ids);
                max_prob = new_prob;
            }
           
        }
        
        if (the_matrix->model_is_symmetric() == TRUE) max_prob*=2.0;
        the_model->set_track_states(opt_track_id, taxa_track_ids);
    }
    
}

#if 0

void write_probs(string prob_file, int taxa_id, WGX_Data *the_homologs, Exchange *curr_exchange, Ploidy_Like_model *the_model)
{
	int i, j;
	ofstream fout;

	fout.open(prob_file);

	fout<<"#\t";

	for(j=0; j<curr_exchange->get_condlike_size(); j++) {
		if (j != curr_exchange->get_condlike_size() -1)
			fout<<num_to_dupl_data(j)<<"\t";
		else
			fout<<num_to_dupl_data(j)<<"\n";
	}


	for(i=0; i<curr_exchange->get_num_localities(); i++) {
		fout<<i<<"\t"<<(*the_homologs)[i][taxa_id].get_gene_obj(0)->get_name()<<"\t";
		if ((*the_homologs)[i][taxa_id].has_duplicate() == TRUE)
			fout<<(*the_homologs)[i][taxa_id].get_gene_obj(1)->get_name()<<"\t";
		else
			fout<<"NONE\t";
		for(j=0; j<curr_exchange->get_condlike_size(); j++) {
			if (j!= curr_exchange->get_condlike_size()-1)
				fout<<the_model->get_site_prob(i, dupl_to_loss_state(j))<<"\t";
			else
				fout<<the_model->get_site_prob(i, dupl_to_loss_state(j))<<"\n";

		}
	}

	fout.close();
}
#endif
