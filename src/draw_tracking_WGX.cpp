#include <iostream>
#include <math.h>
#include <fstream>
#include "exchange.h"
#include "tree.h"
#include "genome_tripl_list.h"
#include "gen_dna_funcs.h"
#include "maxlike.h"
#include "genome_ploidy_like.h"
#include "phylo_model_matrix.h"
#include <string>
#include <sstream>
#include <map>
#include "write_tree.h"
#include "plot.h"

using namespace::std;

enum IMAGE_SIZE {SMALL, MEDIUM, BIGGER, LARGE, MASSIVE};

#define MAX_FRAME 150;
//#define num_per_page 50
void draw_tracking(Exchange *curr_exchange, Tree *the_tree, Phylo_Matrix *the_matrix, Ploidy_Like_model *the_model,  Clade *the_genomes, WGX_Data *the_homologs, int diag_size);

void generate_frames(Exchange *curr_exchange, Tree *the_tree, Phylo_Matrix *the_matrix, Ploidy_Like_model *the_model,  Clade *the_genomes, WGX_Data *the_homologs, std::string *&full_genome_files, std::string *&prefixes);

int draw_region (Exchange *curr_exchange,  Tree *the_tree, Phylo_Matrix *the_matrix, Ploidy_Like_model *the_model,  Clade *the_genomes, WGX_Data *the_homologs, int start, int num_per_page, int focus, string focus_name, string filename, BOOL bitmap, BOOL IPC_call, BOOL rescale, std::stringstream *plot_ss, IMAGE_SIZE mysize, int *coords, int arrow_loc, BOOL arrow_left, int arrow_taxa, string tracked_name, string outname, string *prefix_map);

extern void get_max_prob_pattern(Exchange *curr_exchange, Phylo_Matrix *the_matrix, Ploidy_Like_model *the_model, int site, int &opt_track_id, int *taxa_track_ids, double &prob);

extern string make_gene_tree(Exchange *curr_exchange,  Tree *the_tree, Phylo_Matrix *the_matrix, TREE_TYPE my_type, WGX_Data *the_homologs, Ploidy_Like_model *the_model,  Clade *the_genomes, int pillar_num);


void construct_full_genome_lookup(string genomefile, Clade *the_genomes, std::map<std::string, int> gene_hash, std::map<std::string, int> &left_matches, std::map<std::string, int> &right_matches, std::map<std::string, int> &taxa_lookup, std::map<std::string, std::string> &left_name, std::map<std::string, std::string> &right_name, std::map<std::string, std::string> &aliases);
void make_prefixes (Clade *the_genomes, std::string *prefixes, std::string *&prefix_map);
void prune_tree (Exchange *&curr_exchange, Tree *&current_tree, string prune_name);
void rebuild_tree (Branch *curr_old, Branch **lookup_table, Tree *new_tree, int &tips_so_far);
string extract_name(string input_name);

int draw_model_diag(Phylo_Matrix *the_matrix, string plotfile, BOOL bitmap, BOOL IPC_call, std::stringstream *plot_ss, int matrix_num);

void draw_tracking(Exchange *curr_exchange, Tree *the_tree, Phylo_Matrix *the_matrix, Ploidy_Like_model *the_model,  Clade *the_genomes, WGX_Data *the_homologs, int diag_size)
{
	int i=1, start=0, j, taxa,  dupl_level, stop, num_per_page, arrow_loc=-1, arrow_taxa=0;
    BOOL arrow_left=TRUE;
    string filename, number_str, index_file, tracked_name, outname, *prefix_map=0;
    ofstream fout;

    
    tracked_name="NONE";
    outname="NONE";
    num_per_page=diag_size;
    
    the_model->get_gene_conditional_probs();
    
    index_file="tracking_index.txt";
    
    fout.open(index_file.c_str());
    
    
	do {
        stringstream ss;
        ss << i;
        number_str = ss.str();
        filename = "tracking_section";
        filename += "_";
        filename += number_str;
        filename += ".eps";
	
        if (start+num_per_page >= the_homologs->get_num_homologs())
            stop=the_homologs->get_num_homologs()-1;
            else stop=start+num_per_page;
        for(j=start; j<stop; j++) {
            for(taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) {
                for(dupl_level=0; dupl_level< the_homologs->get_dupl_level(); dupl_level++) {
                    if ((*the_model->the_tracks)[taxa].get_gene_track(j, dupl_level)->my_locus != 0)
                        fout<<(*the_model->the_tracks)[taxa].get_gene_track(j, dupl_level)->my_locus->get_gene_obj((*the_model->the_tracks)[taxa].get_gene_track(j, dupl_level)->index_num)->get_name_string()<<"\t"<<filename<<endl;
                }
            }
                
        }
        
		draw_region (curr_exchange, the_tree, the_matrix, the_model, the_genomes, the_homologs, start, num_per_page, -1, "", filename, FALSE, FALSE, FALSE, 0, LARGE,0, arrow_loc, arrow_left, arrow_taxa,
                     tracked_name, outname, prefix_map);
		start += (int)(0.5*num_per_page);
		i++;
	} while(start < the_homologs->get_num_homologs());
    fout.close();
}


void generate_frames(Exchange *curr_exchange, Tree *the_tree, Phylo_Matrix *the_matrix, Ploidy_Like_model *the_model,  Clade *the_genomes, WGX_Data *the_homologs, std::string *&full_genome_files, std::string *&prefixes)
{
    int i,frame_size, pillar, dupl_level, taxa, start, end, max, coords[7], arrow_loc=-1, arrow_taxa, focus;
    BOOL size_valid, is_bmp, get_tree, gene_error, have_full=FALSE, arrow_left, first;
    string gene_name, quit, size, filename, number_str, type_string, tree_choice, tracked_name, outname, *prefix_map;
    std::map<std::string, int> gene_hash, full_left_matches, full_right_matches, full_taxa_ids;
    std::map<std::string, std::string> left_name, right_name, aliases;
    
    for(pillar=0; pillar<the_homologs->get_num_homologs(); pillar++) {
        for(taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) {
            for(dupl_level=0; dupl_level< the_homologs->get_dupl_level(); dupl_level++) {
                if ((*the_model->the_tracks)[taxa].get_gene_track(pillar, dupl_level)->my_locus != 0)
                gene_hash[(*the_model->the_tracks)[taxa].get_gene_track(pillar, dupl_level)->my_locus->get_gene_obj((*the_model->the_tracks)[taxa].get_gene_track(pillar, dupl_level)->index_num)->get_name_string()]=pillar;
            }
        }
    }
    
    make_prefixes (the_genomes, prefixes, prefix_map);
    
    if (full_genome_files!=0) {
        for(taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) {
            if (full_genome_files[taxa] != "NONE") {
                have_full=TRUE;
                construct_full_genome_lookup(full_genome_files[taxa], the_genomes, gene_hash, full_left_matches, full_right_matches, full_taxa_ids, left_name, right_name, aliases);
            }
        }
    }
    
    gene_name="";
    quit=gene_name;
    the_model->get_gene_conditional_probs();
    while(quit != "exit") {
        
    
        cout<<"Enter the gene to plot:\n";
        //cin>>gene_name;
        //quit=gene_name;
        //std::transform(quit.begin(), quit.end(), quit.begin(),::tolower);
        
        gene_error=TRUE;
        first=TRUE;
        
        while((quit != "exit") && (gene_error==TRUE)) {
            if (first == FALSE)
                cout<<"ERROR: gene "<<gene_name<<" not found in dataset or full genome synteny blocks. Please re-enter:\n";
            
            first=FALSE;
            tracked_name="NONE";
            outname="NONE";
            cin>>gene_name;
            quit=gene_name;
            std::transform(quit.begin(), quit.end(), quit.begin(),::tolower);
            
            if (gene_hash.find(gene_name) == gene_hash.end()) {
                if (have_full ==TRUE) {
                    if (!(full_left_matches.find(gene_name) == full_left_matches.end())) {
                        if (full_left_matches[gene_name] != -1) {
                            pillar=full_left_matches[gene_name];
                            tracked_name=left_name[gene_name];
                            outname=gene_name;
                            arrow_left=TRUE;
                            arrow_taxa=full_taxa_ids[gene_name];
                            arrow_loc=pillar;
                            gene_error=FALSE;
                            focus=-1;
                        }
                        else {
                            if (full_right_matches[gene_name] != -1) {
                                pillar=full_right_matches[gene_name];
                                tracked_name=right_name[gene_name];
                                outname=gene_name;
                                arrow_left=FALSE;
                                arrow_loc=pillar;
                                arrow_taxa=full_taxa_ids[gene_name];
                                gene_error=FALSE;
                                focus=-1;
                            }
                        }
                    }
                    
                }
            }
            else {
                pillar = gene_hash[gene_name];
                arrow_loc=-1;
                focus=pillar;
                gene_error=FALSE;
            }
            
        }
        
        
        
        
        if (quit != "exit") {
            cout<<"Enter the window size:\n";
            cin>>size;
            size_valid=TRUE;
            
            try {
                frame_size=std::stoi(size,nullptr);
            }
            catch (std::invalid_argument const &e) {size_valid=FALSE;}
            catch (std::out_of_range const &e) {size_valid=FALSE;}
            
            while(size_valid==FALSE) {
                cout<<"ERROR: Invalid frame size. Please re-enter:\n";
                cin>>size;
                size_valid=TRUE;
                
                try {
                    frame_size=std::stoi(size,nullptr);
                }
                catch (std::invalid_argument const &e) {size_valid=FALSE;}
                catch (std::out_of_range const &e) {size_valid=FALSE;}
            }
            
            frame_size=abs(frame_size);
            
            max=MAX_FRAME;
            if (frame_size > max) {
                frame_size=MAX_FRAME;
                cout<<"Note: maximum frame size exceeded. Reseting to "<<frame_size<<"\n";
            }
            
            
            cout<<"Bitmap image? (Y/N):\n";
            cin>>type_string;
            
            std::transform(type_string.begin(), type_string.end(), type_string.begin(),::tolower);
            
            if (type_string.front() == 'y') is_bmp=TRUE;
            else is_bmp=FALSE;
            
            cout<<"Generate gene tree? (Y/N):\n";
            cin>>tree_choice;
            
            std::transform(tree_choice.begin(), tree_choice.end(), tree_choice.begin(),::tolower);
            
            if (tree_choice.front() == 'y') {
                get_tree=TRUE;
                cout<<"Tree string is "<<make_gene_tree(curr_exchange,  the_tree, the_matrix, NEXUS_TREE, the_homologs, the_model,  the_genomes, pillar)<<endl;
            }
            else get_tree=FALSE;
            
            start = pillar - (frame_size /2);
            
            if (start < 0) start=0;
            
            end=start+frame_size;
            
            stringstream ss;
            ss << frame_size;
            number_str = ss.str();
            filename = gene_name;
            filename += "tracking_";
            filename += number_str;
            if (is_bmp == TRUE)
                filename += ".png";
            else
                filename += ".eps";
            
            //draw_region (curr_exchange, the_tree, the_matrix, the_model, the_genomes, the_homologs, start, frame_size, focus, gene_name, filename, is_bmp, FALSE, 0, LARGE,coords, arrow_loc, arrow_left, arrow_taxa, tracked_name, outname, prefix_map);
            draw_region (curr_exchange, the_tree, the_matrix, the_model, the_genomes, the_homologs, start, frame_size, focus, gene_name, filename, is_bmp, FALSE, FALSE, 0, MASSIVE,coords, arrow_loc, arrow_left, arrow_taxa, tracked_name, outname, prefix_map);
            
            for(i=0; i<7; i++) cout<<coords[i]<<"\t";
            cout<<endl;
        }
        
    }
}


#if 0
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
        cout<<"Have potential tree: nf: "<<num_full<<" nz: "<<num_zero<<" st: "<<single_tree<<" ft: "<<full_tree<<endl;
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
                
                cout<<"Renaming taxa "<<(*new_tree)[i]->get_name()<<" to "<<new_taxa_name<<endl;
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
        cout<<"Alloc exchange"<<endl<<flush;
        (*pruned_exchange)=(*new_exchange);
        cout<<"Alloc/Copy exchange with np= "<<num_prune<<endl<<flush;
        if (three_two_tree==TRUE)
            pruned_exchange->set_num_taxa(new_exchange->get_num_taxa()-(num_prune-the_genomes->get_num_genomes()));
        else
            pruned_exchange->set_num_taxa(new_exchange->get_num_taxa()-num_prune);
        pruned_tree=new Tree(pruned_exchange, TRUE);
        cout<<"Alloc tree"<<endl<<flush;
       
  
    
        
        for(taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) {
            for (level=0; level<the_homologs->get_dupl_level(); level++) {
                my_track=0;
                while(the_model->tracking_permutes[taxa_track_ids[taxa]][my_track] != level) my_track++;
                
                if (three_two_tree==TRUE) my_track=tree_nums[my_track];
                
                cout<<"Taxa "<<taxa<<" Level: "<<level<<" My track: "<<my_track<<endl<<flush;
                
                if (my_track != -1) {
                    std::ostringstream ss;
                    ss<<(*the_genomes)[taxa].get_name_string()<<my_track;
                    new_taxa_name=ss.str();
                    
                    if ((*the_model->the_tracks)[taxa].get_gene_track(pillar_num, level)->my_locus == 0) {
                        cout<<"Taxa "<<new_taxa_name<<" is empty: pruning\n";
                        j=0;
                        while(strcmp((*new_tree)[j]->get_name(), new_taxa_name.c_str()) !=0) j++;
                        
                        new_tree->remove_tip((*new_tree)[j]);
                    }
                    else {
                        brn=0;
                        while((*new_tree)[brn]->get_name() != new_taxa_name) brn++;
                        cout<<"Renaming "<<new_taxa_name<<" to "<<(*the_model->the_tracks)[taxa].get_gene_track(pillar_num, level)->my_locus->get_gene_obj((*the_model->the_tracks)[taxa].get_gene_track(pillar_num, level)->index_num)->get_name_string()<<endl;
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
                cout<<"Setting "<<brn<<" as "<<(*pruned_tree)[brn]->get_parent()<<" to sib "<<my_sib->get_name()<<" par: "<<my_par->get_name()<<endl<<flush;
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
        
        cout<<"Branch "<<j<<" name: "<<(*pruned_tree)[j]->get_name()<<" ID "<<(*pruned_tree)[j]<<" P: "<<(*pruned_tree)[j]->get_parent()<<" C1: "<<(*pruned_tree)[j]->get_child(0)<<" C2: "<<(*pruned_tree)[j]->get_child(1);
        if ((*pruned_tree)[j]->is_tip())
            cout<<" Tip: "<<(*pruned_tree)[j]->get_taxa_id()<<endl;
        else
            cout<<" Internal\n";
    
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
#endif
#if 0

string make_gene_tree(Exchange *curr_exchange,  Tree *the_tree, Phylo_Matrix *the_matrix, TREE_TYPE my_type, WGX_Data *the_homologs, Ploidy_Like_model *the_model,  Clade *the_genomes, int pillar_num)
{
    int i, subtree, taxa, level, my_track, max_prob_pattern, *taxa_track_ids, brn, *gene_cnts, num_full, num_zero, num_prune=0, tree_nums[3];
    double max_prob;
    BOOL keep_tree=FALSE, three_two_tree=FALSE;
    string tree_string, new_taxa_name;
    Exchange *new_exchange;
    Tree  *new_tree;
    Write_Tree *writeout_tree;
    
    taxa_track_ids=new int[the_genomes->get_num_genomes()];
    gene_cnts=new int [the_homologs->get_dupl_level()];
    
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
    
    if ((num_full==1) && (num_zero == (the_homologs->get_dupl_level()-1))) {
        new_tree = new_tree = new Tree(curr_exchange, TRUE);
        (*new_tree)=(*the_tree);
        keep_tree=TRUE;
    }
    else {
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
    for(i=0; i<new_exchange->get_num_branches(); i++) {
        if ((*new_tree)[i]->is_tip() ==TRUE) {
            subtree = (*new_tree)[i]->get_taxa_id()/curr_exchange->get_num_taxa();
            
            
            //new_taxa_name=(*new_tree)[i]->get_name();
            std::ostringstream ss;
            ss<<(*new_tree)[i]->get_name()<<subtree;
            new_taxa_name=ss.str();
            
            cout<<"Renaming taxa "<<(*new_tree)[i]->get_name()<<" to "<<new_taxa_name<<endl;
            (*new_tree)[i]->set_name(new_taxa_name.c_str());
        }
    }
    
    
    get_max_prob_pattern(curr_exchange, the_matrix, the_model, pillar_num, max_prob_pattern, taxa_track_ids, max_prob);
    
    
    for(taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) {
        for (level=0; level<the_homologs->get_dupl_level(); level++) {
            my_track=0;
            while(the_model->tracking_permutes[taxa_track_ids[taxa]][my_track] != level) my_track++;
            
            if (three_two_tree==TRUE) my_track=tree_nums[my_track];
            
            std::ostringstream ss;
            ss<<(*the_genomes)[taxa].get_name_string()<<my_track;
            new_taxa_name=ss.str();
            
            if ((*the_model->the_tracks)[taxa].get_gene_track(pillar_num, level)->my_locus == 0) {
                cout<<"Taxa "<<new_taxa_name<<" is empty: pruning\n";
                prune_tree (new_exchange, new_tree, new_taxa_name);
            }
            else {
                brn=0;
                while((*new_tree)[brn]->get_name() != new_taxa_name) brn++;
                
                (*new_tree)[brn]->set_name((*the_model->the_tracks)[taxa].get_gene_track(pillar_num, level)->my_locus->get_gene_obj((*the_model->the_tracks)[taxa].get_gene_track(pillar_num, level)->index_num)->get_name_string().c_str());
            }
        }
    }
    
    
    
    if (my_type == NEXUS_TREE)
        writeout_tree=new Write_Nexus_Tree();
    else
        writeout_tree=new Write_Phylip_Tree();

    
    writeout_tree->write_tree_to_string(tree_string, "POInT Gene Tree", new_tree, new_exchange);
    
    delete writeout_tree;
    if (keep_tree==FALSE) {
        delete new_exchange;
        
    }
    delete new_tree;
    delete[] gene_cnts;
    delete[] taxa_track_ids;
    return(tree_string);
}






void prune_tree (Exchange *&curr_exchange, Tree *&current_tree, string prune_name)
{
    int j, k, new_branch_id, init_tips;
    BOOL valid;
    Exchange *new_exchange;
    Branch **lookup_table;
    Tree *new_tree;
    
    new_exchange=new Exchange();
    new_exchange->set_num_taxa(curr_exchange->get_num_taxa()-1);
    new_exchange->set_num_sites(curr_exchange->get_num_sites());
    
    lookup_table=new Branch*[curr_exchange->get_num_branches()];
    
    
    new_tree=new Tree(new_exchange, TRUE);
    cout<<"Pruning "<<prune_name<<endl;
    j=0;
    while((j<curr_exchange->get_num_taxa()) && (strcmp(prune_name.c_str(), (*current_tree)[j]->get_name()) !=0)) j++;
    
    if (j<curr_exchange->get_num_branches()) {
        current_tree->remove_tip((*current_tree)[j]);
        

        cout<<"Creating lookup table\n";
        new_branch_id=0;
        for(k=0; k<curr_exchange->get_num_branches(); k++) {
            valid=TRUE;
            if (((*current_tree)[k]->get_sibling()==0) && ((*current_tree)[k] != current_tree->find_root())) valid=FALSE;
            if (((*current_tree)[k]->is_tip() == FALSE ) && (((*current_tree)[k]->get_child(1) ==0) || ((*current_tree)[k]->get_child(0) ==0))) valid=FALSE;
            
            if (valid==TRUE) {
                cout<<"Assigning "<<(*current_tree)[k]->get_brn_num()<<" to new "<<new_branch_id<<endl;
                lookup_table[(*current_tree)[k]->get_brn_num()] =(*new_tree)[new_branch_id];
                new_branch_id++;
            }
        }
        init_tips=0;
        
        rebuild_tree (current_tree->find_root(), lookup_table, new_tree, init_tips);
        new_tree->set_root((*new_tree)[0]);

        
        new_tree->diganose_tree(new_tree->find_root());
        
        for(k=0; k<new_exchange->get_num_branches(); k++) {
            if (((*new_tree)[k] != new_tree->find_root()) && ((*new_tree)[k] != new_tree->find_null_branch_id())) {
                if ((*new_tree)[k]->expect_subs_site() ==0) {
                    (*new_tree)[k]->set_expect_subs_site(0.1);
                }
            }
            
            cout<<"Branch "<<k<<" name: "<<(*new_tree)[k]->get_name()<<" P: "<<(*new_tree)[k]->get_parent()<<" C1: "<<(*new_tree)[k]->get_child(0)<<" C2: "<<(*new_tree)[k]->get_child(1);
            if ((*new_tree)[k]->is_tip())
            cout<<" Tip: "<<(*new_tree)[k]->get_taxa_id()<<endl;
            else
            cout<<" Internal\n";
            
            delete curr_exchange;
            delete current_tree;
            current_tree=new_tree;
            curr_exchange=new_exchange;
            
        }
    }
    else {
        cerr<<"ERROR: Could not prune "<<prune_name<<" from this tree\n";
        delete new_exchange;
        delete new_tree;
    }
    
    delete[] lookup_table;
       
}

void rebuild_tree (Branch *curr_old, Branch **lookup_table, Tree *new_tree, int &tips_so_far)
{
    Branch *new_me, *new_child1, *new_child2, *old_child1, *old_child2;
    if (curr_old->is_tip() ==FALSE) {
        old_child1=curr_old->get_child(0);
        old_child2=curr_old->get_child(1);
        new_me=lookup_table[curr_old->get_brn_num()];
        new_child1=lookup_table[old_child1->get_brn_num()];
        new_child2=lookup_table[old_child2->get_brn_num()];
        cout<<"Old brnaches"<<curr_old->get_name()<<" C1: "<<old_child1->get_name()<<" C2: "<<old_child2->get_name()<<endl;
        cout<<"New branch ids "<<new_me->get_brn_num()<<" C1 "<<new_child1->get_brn_num()<<" C2: "<<new_child2->get_brn_num()<<endl;
        
        
        
        new_tree->set_as_parent_child(new_me, new_child1, 0);
        new_tree->set_as_parent_child(new_me, new_child2, 1);
        new_tree->set_as_siblings(new_child1, new_child2);
        
        new_me->set_expect_subs_site(curr_old->get_brnlen());
        new_me->set_name(curr_old->get_name());
        new_me->set_tip(FALSE);
        
        if (curr_old->get_parent() ==0)
        new_me->set_parent(0);
        
        
        new_child1->set_expect_subs_site(old_child1->get_brnlen());
        new_child1->set_name(old_child1->get_name());
        new_child1->set_tip(old_child1->is_tip());
        
        new_child2->set_expect_subs_site(old_child2->get_brnlen());
        new_child2->set_name(old_child2->get_name());
        new_child2->set_tip(old_child2->is_tip());
        
        if (old_child1->is_tip() == TRUE) {
            new_child1->set_taxa_id(tips_so_far);
            tips_so_far++;
        }
        if (old_child2->is_tip() == TRUE) {
            new_child2->set_taxa_id(tips_so_far);
            tips_so_far++;
        }
        
        if (old_child1->is_tip() == FALSE)
        rebuild_tree(old_child1, lookup_table, new_tree, tips_so_far);
        
        if (old_child2->is_tip() == FALSE)
        rebuild_tree(old_child2, lookup_table, new_tree, tips_so_far);
        
    }
    
}

#endif
    
#define LINE_SIZE 7

int draw_region (Exchange *curr_exchange,  Tree *the_tree, Phylo_Matrix *the_matrix, Ploidy_Like_model *the_model,  Clade *the_genomes, WGX_Data *the_homologs, int start, int num_per_page, int focus, string focus_name, string filename, BOOL bitmap, BOOL IPC_call, BOOL rescale, std::stringstream *plot_ss, IMAGE_SIZE mysize, int *coords, int arrow_loc, BOOL arrow_left, int arrow_taxa, string tracked_name, string outname, string *prefix_map)
{
    int level, i, j, k, thandle, taxa, stop, num_genes, *max_prob_pattern, *last_max_prob_pattern, track_id, slen, find_loc_pillar, num_empty, num_full, my_line_size, base_line_size,
        find_loc_track, **taxa_track_ids, other_level, y_loc, my_track, name_line, num_lines, name_pos, cnt_dupl, *cnt_state, cnt_1, cnt_2, *color_ids, realy, realx;
    double spacex, spacey, border=10, spacer, name_font_size=1.0, center,
      xcenter, ycenter, label_size, width, boxsize, *tracky, tree_y, tree_o,
    box_x_foot, box_y_foot, text_div, font_size, old_fontsize, my_offset, spacer_x,
    *lastline, offset, offset_y, *max_prob, new_prob, text_width, subtrack_size, font, x_y_aspect, box_dim;
    char write_string[100], new_maxlen, attrib[100], value[100], size_string[25], readchar;
    string fullname, subname, split_name[5], new_taxa_name, colors[10], prefix, name_head, test_string;
    Gene_Track_List_DX ***figure_points, **last_drawn=0;
    BOOL *has_last, *track_first, vertical=TRUE;
    char* buf_pt;
    FILE *outfile;
    std::size_t found, size, mypos;
    
    
    colors[0]="mistyrose";
    colors[1]="lightsteelblue1";
    colors[2]="springgreen3";
    colors[3]="bisque3";
    colors[4]="lightsalmon2";
    colors[5]="lavender";
    colors[6]="palegoldenrod";
    colors[7]="wheat1";
    colors[8]="darkolivegreen3";
    
    
    
    last_drawn=new Gene_Track_List_DX *[the_homologs->get_dupl_level()];
    
    for(level=0; level< the_homologs->get_dupl_level(); level++) last_drawn[level]=0;
    
    taxa_track_ids=new int* [num_per_page];
    
    for(j=0; j<num_per_page; j++)
        taxa_track_ids[j]=new int[the_genomes->get_num_genomes()];
    
    last_max_prob_pattern = new int [num_per_page];
    max_prob_pattern = new int [num_per_page];
    max_prob = new double [num_per_page];
    tracky=new double [the_homologs->get_dupl_level()];
    has_last=new BOOL [the_homologs->get_dupl_level()];
    track_first = new BOOL [the_homologs->get_dupl_level()];
    lastline=new double [the_homologs->get_dupl_level()];
    figure_points = new Gene_Track_List_DX**[the_homologs->get_dupl_level()];
    for (level=0; level<the_homologs->get_dupl_level(); level++)
        figure_points[level]=new Gene_Track_List_DX * [num_per_page];
    
    if (start+num_per_page > the_homologs->get_num_homologs())
      stop=the_homologs->get_num_homologs()-1;
    else
      stop=start+num_per_page-1;

    num_genes=stop-start+1;

    color_ids=new int[num_per_page];
    cnt_state=new int [2*the_homologs->get_dupl_level()-1];
  //cout<<"Printing "<<num_genes<<" in this file"<<" Start: "<<start<<" Stop: "<<stop<<endl;

    if (IPC_call==FALSE)
        outfile=fopen(filename.c_str(), "w");
    else {
        size=0;
        //printf ("At open buf = `%s', size = %zu\n", buf_pt, size);
        outfile= open_memstream (&buf_pt, &size);
        
        //fprintf(outfile, "TESTING");
        fflush(outfile);
        //printf ("At test buf = `%s', size = %zu\n", buf_pt, size);
    }
   
    spacex=1200.0;
    spacey=600.0;
    offset=spacex*0.13;
    offset_y=spacey*0.07;
    //offset=offset_y=0;
    box_x_foot=((spacex-2.0*border-offset)/(double)num_per_page);
    
    boxsize=0.83*box_x_foot;
    spacer_x = (1.0/2.0)*(box_x_foot-boxsize);
    
    cout<<"Box size is "<<boxsize<<" foot is "<<box_x_foot<<" Spacer is "<<spacer_x<<endl;
    
    subtrack_size=0.835/(double)the_homologs->get_dupl_level();
    
    spacer = 0.02/((double)the_homologs->get_dupl_level()-1.0)*spacey;
    
    box_y_foot=((subtrack_size/the_genomes->get_num_genomes())*spacey);
    width=0.82*box_y_foot;
    //tracky[0]=spacey-0.3*spacey;
    tracky[0]=spacey-offset_y-width;
    
    x_y_aspect=boxsize/width;
    
    if (x_y_aspect>1.0) {
        vertical=FALSE;
        box_dim=boxsize;
    }
    else {
        vertical=TRUE;
        box_dim=width;
    }
    
   
    
    
    if (IPC_call==FALSE)
        cout<<"Y for "<<0<<" is "<<tracky[0]<<" :" <<(subtrack_size*level)*spacey<<endl;
    
    for(level=1; level< the_homologs->get_dupl_level(); level++) {
        
        //tracky[level]=0+(subtrack_size*level)*spacey;
        tracky[level]=tracky[level-1]-(subtrack_size*spacey+spacer);
        if (IPC_call==FALSE)
            cout<<"Y for "<<level<<" is "<<tracky[level]<<" :" <<(subtrack_size*level)*spacey<<endl;
    }

    if (bitmap == FALSE) {
        strcpy(attrib, "PAGESIZE");
        strcpy(value, "letter,xsize=10in,ysize=5in");
        pl_parampl (attrib, value);
        /* create a Postscript Plotter that writes to standard output */
        
        if ((thandle = pl_newpl ("ps", stdin, outfile, stderr)) < 0) {
              fprintf (stderr, "Couldn't create Plotter\n");
              return 1;
        }
       
        pl_selectpl (thandle);       /* select the Plotter for use */
    }
    else {
        switch (mysize) {
            case SMALL:
                strcpy(size_string, "900x450");
                realx=900;
                realy=450;
                break;
            case MEDIUM:
                strcpy(size_string, "1200x600");
                realx=1200;
                realy=600;
                break;
            case BIGGER:
                strcpy(size_string, "2400x1200");
                realx=2400;
                realy=1200;
                break;
            case LARGE:
            default:
                strcpy(size_string, "1800x900");
                realx=1800;
                realy=900;
                break;
            case MASSIVE:
                strcpy(size_string, "3600x1800");
                realx=3600;
                realy=1800;
                break;
                
        }
        /* set a Plotter parameter */
        pl_parampl ("BITMAPSIZE", size_string);
       
        if ((thandle = pl_newpl ("png", stdin, outfile, stderr)) < 0) {
            fprintf (stderr, "Couldn't create Plotter\n");
            return 1;
        }
       
        pl_selectpl (thandle);       /* select the Plotter for use */
    }

    
    
    if (pl_openpl () < 0)       /* open Plotter */ {
      fprintf (stderr, "Couldn't open Plotter\n");
      return 1;
    }
    
    pl_fspace (0.0, 0.0, spacex, spacey); /* specify user coor system */
    pl_flinewidth (1.0);       /* line thickness in user coordinates */
    pl_pencolorname ("black");    /* path will be drawn in red */
    pl_erase ();                /* erase Plotter's graphics display */
    /*pl_fmove (600.0, 300.0);*/    /* position the graphics cursor */
    if (bitmap == TRUE) {
        pl_fontname ("HersheySans-BoldOblique"); /* choose a Postscript font */
        pl_ftextangle (0);         /* text inclination angle (degrees) */
        pl_ffontsize (15);
        font_size=15.0;
    }
    else {
        pl_fontname("Helvetica-Oblique");
        pl_ftextangle (0);         /* text inclination angle (degrees) */
        pl_ffontsize (14);
        font_size=14;
    }
  
    for(level=0; level<the_homologs->get_dupl_level(); level++) {
        for(i=0; i<the_genomes->get_num_genomes(); i++) {
            fullname=(*the_genomes)[i].get_name_string();
           
            found=fullname.find("thaliana");
            if (found!=std::string::npos)
                fullname=fullname.substr(0, found+8);
            
            found=fullname.find("_rerun");
            if (found!=std::string::npos)
            fullname=fullname.substr(0, found);
            else new_taxa_name=fullname;
            
            found=fullname.find("_v3");
            if (found!=std::string::npos)
                fullname=fullname.substr(0, found);
            else new_taxa_name=fullname;
            
            found=fullname.find("_V2");
            if (found!=std::string::npos)
                fullname=fullname.substr(0, found);
            else new_taxa_name=fullname;
            
            std::size_t found = fullname.find("_");
            if (found!=std::string::npos) {
                //std::cout << "first 'needle' found at: " << found << '\n';
                new_taxa_name = fullname.front();
                new_taxa_name += ". ";
                new_taxa_name += fullname.substr(found+1, fullname.length()-found);
            }
            
            std::transform(new_taxa_name.begin(), new_taxa_name.begin()+1, new_taxa_name.begin(), ::toupper);
           
            text_width = pl_flabelwidth (new_taxa_name.c_str());
            pl_fmove (0.5*text_width+border, tracky[level]-i*box_y_foot+0.5*width);
            pl_alabel('c', 'c', new_taxa_name.c_str());
#if 1
            if (prefix_map !=0) {
                if (prefix_map[i] != "") {
                    prefix = "Prefix: " + prefix_map[i];
                    
                    if (bitmap == TRUE)
                        pl_fontname ("HersheySans");
                    
                    else
                        pl_fontname("Helvetica");
                    
                    old_fontsize=font_size;
                    font_size=0.65*font_size;
                    pl_ffontsize (font_size);
                    text_width = pl_flabelwidth (prefix.c_str());
                    pl_fmove (0.5*text_width+border+20, tracky[level]-i*box_y_foot-0.8*old_fontsize+0.5*width);
                    pl_alabel('c', 'c', prefix.c_str());
                    font_size=old_fontsize;
                    pl_ffontsize (font_size);
                    
                    if (bitmap == TRUE)
                        pl_fontname ("HersheySans-BoldOblique");
                    else
                        pl_fontname("Helvetica-Oblique");
                    
                }
            }
#endif
            
        }
    }

    my_line_size=LINE_SIZE;
    //if (the_genomes->get_num_genomes()>6) my_line_size -=2;
    
    if (the_genomes->get_num_genomes()<5) my_line_size++;
    
    base_line_size=my_line_size;
    
    test_string="";
    for(i=0; i<my_line_size; i++) {
        if (i%2 ==0)
            test_string+="W";
        else
            test_string+="N";
    }
    if (bitmap == TRUE) {
        pl_fontname ("HersheySans-Bold");
        pl_ffontsize (name_font_size);
    }
    else {
        pl_fontname("Helvetica");
         pl_ffontsize (name_font_size);
    }
    
    text_width = pl_flabelwidth (test_string.c_str());
    
    while(text_width < (0.96*box_dim)) {
        name_font_size+=0.05;
        pl_ffontsize (name_font_size);
        text_width = pl_flabelwidth (test_string.c_str());
    }
    
    cout<<"Final text size: "<<name_font_size<<endl;
    
    for(j=0; j<num_genes; j++) {
        if (bitmap == TRUE) {
            pl_fontname ("HersheySans-Bold");
            //if (the_genomes->get_num_genomes()<8)
                pl_ffontsize (12);
            //else {
            //    if (the_genomes->get_num_genomes()<10)
            //        pl_ffontsize (10);
            //    else
            //        pl_ffontsize (8);
           // }
        }
        else {
            pl_fontname("Helvetica");
            pl_ffontsize (12);
        }
        
        pl_ftextangle(90.0);
        
        get_max_prob_pattern(curr_exchange, the_matrix, the_model, start+j, max_prob_pattern[j], taxa_track_ids[j], max_prob[j]);
        double_to_string(write_string, 49, 3, max_prob[j]);
        text_width = pl_flabelwidth (write_string);
        pl_fmove(j*box_x_foot+offset+0.5*boxsize, spacey-0.6*text_width);
        pl_alabel('c', 'c', write_string);
        
        cnt_dupl=0;
        
        for (level=0; level<the_homologs->get_dupl_level(); level++) cnt_state[level]=0;

        for(taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) {
            for (level=0; level<the_homologs->get_dupl_level(); level++) {
                my_track=0;
                while(the_model->tracking_permutes[taxa_track_ids[j][taxa]][my_track] != level) my_track++;
                if ((*the_model->the_tracks)[taxa].get_gene_track(start+j, level)->my_locus != 0) {
                    cnt_dupl++;
                    cnt_state[my_track]++;
                }
            }
        
        }
        
        if (the_homologs->get_dupl_level()==2) {
        
            if (cnt_dupl == (2*the_genomes->get_num_genomes())) color_ids[j]=0;
            else {
                if ((cnt_state[0] == the_genomes->get_num_genomes()) &&( cnt_state[1]==0)) {
                    if (the_matrix->model_is_symmetric()==FALSE)
                        color_ids[j] =1;
                    else
                        color_ids[j] =8;
                }
                else {
                    if ((cnt_state[1] == the_genomes->get_num_genomes()) &&( cnt_state[0]==0)) {
                        if (the_matrix->model_is_symmetric()==FALSE)
                            color_ids[j] =2;
                        else
                            color_ids[j] =8;
                    }
                    else color_ids[j]=7;
                }
            }
        }
        
        if (the_homologs->get_dupl_level()==3) {
            num_empty=0;
            num_full=0;
            for (level=0; level<the_homologs->get_dupl_level(); level++) {
                if (cnt_state[level]==0) num_empty++;
                if (cnt_state[level]==the_genomes->get_num_genomes()) num_full++;
            }
            
            color_ids[j]=7;
            
            if (num_full == 3)
                color_ids[j]=5;
            else {
                if ((num_empty ==1)  && (num_full==2))
                    color_ids[j]=0;
                else {
                    if (num_empty ==2) {
                        if (cnt_state[0] == the_genomes->get_num_genomes()) color_ids[j]=1;
                        if (cnt_state[1] == the_genomes->get_num_genomes()) color_ids[j]=2;
                        if (cnt_state[2] == the_genomes->get_num_genomes()) color_ids[j]=4;
                    }
                }
            }
        }

    }

    //get_max_prob_pattern(curr_exchange, the_matrix, the_model, start, max_prob_pattern, taxa_track_ids[0], max_prob);
    for(taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) {
        for(j=0; j<num_genes; j++) {
            for(level=0; level<the_homologs->get_dupl_level(); level++)
                figure_points[level][j]=0;
        }
      
        for(level=0; level<the_homologs->get_dupl_level(); level++)
            track_first[level]=TRUE;
      
        for (level=0; level<the_homologs->get_dupl_level(); level++) {
            if ((*the_model->the_tracks)[taxa].get_gene_track(start, the_model->tracking_permutes[taxa_track_ids[0][taxa]][level])->last != 0) {
                has_last[the_model->tracking_permutes[taxa_track_ids[0][taxa]][level]]=TRUE;
                lastline[the_model->tracking_permutes[taxa_track_ids[0][taxa]][level]]=offset-0.1*offset;
            }
            else
                has_last[the_model->tracking_permutes[taxa_track_ids[0][taxa]][level]]=FALSE;
        }
        
        tree_y=9;
        tree_o=8;

        if (coords != 0) {
            if (rescale==FALSE) {
                coords[0]=(int)((offset/spacex)*(double)realx);
                coords[1]=(int)(((boxsize)/spacex)*(double)realx);
                coords[2]= (int)(((box_x_foot-boxsize)/spacex)*(double)realx);
                coords[3] = realy-(int)(((tracky[0]+width)/spacey)*(double)realy);
                coords[4] = realy-(int)(((tree_o+5.0*tree_y)/spacey)*(double)realy);
                coords[5]= realy-(int)(((tree_o+3.0*tree_y)/spacey)*(double)realy);
                coords[6]= realy-(int)(((tree_o)/spacey)*(double)realy);
            }
            else {
                coords[0]=(int)((offset/spacex)*(double)(realx/2.0));
                coords[1]=(int)(((boxsize)/spacex)*(double)(realx/2.0));
                coords[2]= (int)(((box_x_foot-boxsize)/spacex)*(double)(realx/2.0));
                coords[3] = (int)(realy/2.0)-(int)(((tracky[0]+width)/spacey)*(double)(realy/2.0));
                coords[4] = (int)(realy/2.0)-(int)(((tree_o+5.0*tree_y)/spacey)*(double)(realy/2.0));
                coords[5]= (int)(realy/2.0)-(int)(((tree_o+3.0*tree_y)/spacey)*(double)(realy/2.0));
                coords[6]= (int)(realy/2.0)-(int)(((tree_o)/spacey)*(double)(realy/2.0));
            }
        }
        
        
        if ((IPC_call==TRUE) &&(bitmap == TRUE)) {
            for(j=0; j<num_genes; j++) {
                pl_fline(j*box_x_foot+offset+(0.15*boxsize), tree_o+3.0*tree_y, j*box_x_foot+(0.4*boxsize)+offset, tree_o+3.0*tree_y);
                pl_fline(j*box_x_foot+offset+(0.15*boxsize), tree_o+3.0*tree_y, j*box_x_foot+offset+(0.15*boxsize), tree_o+2.0*tree_y);
                pl_fline(j*box_x_foot+offset+(0.15*boxsize), tree_o+2.0*tree_y,  j*box_x_foot+(0.45*boxsize)+offset, tree_o+2.0*tree_y);
                pl_fline(j*box_x_foot+offset+(0.05*boxsize), tree_o+2.5*tree_y, j*box_x_foot+(0.15*boxsize)+offset, tree_o+2.5*tree_y);
                pl_fline(j*box_x_foot+offset+(0.05*boxsize), tree_o+1.0*tree_y, j*box_x_foot+(0.425*boxsize)+offset, tree_o+1.0*tree_y);
                pl_fline(j*box_x_foot+offset, tree_o+1.75*tree_y, j*box_x_foot+(0.05*boxsize)+offset, tree_o+1.75*tree_y);
                pl_fline(j*box_x_foot+offset+(0.05*boxsize), tree_o+2.5*tree_y, j*box_x_foot+offset+(0.05*boxsize),tree_o+1.0*tree_y);
                
                if (bitmap == TRUE) {
                    pl_fontname ("HersheySans-Bold");
                    pl_ffontsize (13);
                }
                else {
                    pl_fontname("Helvetica");
                    pl_ffontsize (9.5);
                }
                
                pl_ftextangle(90.0);
                pl_fmove(j*box_x_foot+offset+0.7*boxsize, tree_o+2.0*tree_y);
                pl_alabel('c', 'c', "CDS");
                
            }
            //}
            
            if (bitmap == TRUE) {
                pl_fontname ("HersheySans-Bold");
            }
            else {
                pl_fontname("Helvetica");
            }
        }
        
        for(j=0; j<num_genes; j++) {
            //get_max_prob_pattern(curr_exchange, the_matrix, the_model, start+j, max_prob_pattern, taxa_track_ids[j], max_prob);
            if (IPC_call==FALSE)
              cout<<"Taxa "<<taxa<<" Site "<<j<<"\t"<<start+j<<" max pattern is "<<max_prob_pattern[j]<<" id: "
            <<taxa_track_ids[j][taxa]<<" Prob: "<<max_prob[j]<<endl;
            
	      
            for(level=0; level<the_homologs->get_dupl_level(); level++) {
                my_track=0;
                while(the_model->tracking_permutes[taxa_track_ids[j][taxa]][my_track] != level) my_track++;
                if (IPC_call==FALSE)
                    cout<<"Track for level "<<level<<" is "<<my_track<<endl;
                
                
                if ((*the_model->the_tracks)[taxa].get_gene_track(start+j, level)->my_locus != 0) {
                    
                    
                    if ((*the_model->the_tracks)[taxa].get_gene_track(start+j, level)->last != 0) {
                        find_loc_pillar=-1;
                        for(k=0; k<j; k++) {
                            for(other_level=0; other_level<the_homologs->get_dupl_level(); other_level++) {
                              if (figure_points[other_level][k] == (*the_model->the_tracks)[taxa].get_gene_track(start+j, level)->last) {
                                  find_loc_pillar=k;
                                  find_loc_track=other_level;
                              }
                            }
                        }
                        
                        
                        if (find_loc_pillar != -1) {
                                pl_fline(find_loc_pillar*box_x_foot+boxsize+offset, tracky[find_loc_track]-taxa*box_y_foot+width/2.0,
                                             j*box_x_foot+offset,tracky[my_track]-taxa*box_y_foot+width/2.0);
                        }
                        else {
                            pl_fline(offset-0.1*offset, tracky[my_track]-taxa*box_y_foot+width/2.0,
                                     j*box_x_foot+offset,tracky[my_track]-taxa*box_y_foot+width/2.0);
                        }
                    }
              
                    if (start+j == arrow_loc) {
                        if ((arrow_taxa==taxa) &&
                            (tracked_name == (*the_model->the_tracks)[taxa].get_gene_track(start+j, level)->my_locus->get_gene_obj((*the_model->the_tracks)[taxa].get_gene_track(start+j, level)->index_num)->get_name_string())) {
                            
                            if (bitmap==FALSE) {
                                pl_ffontsize (5.5);
                            }
                            else {
                                font_size=((3.5/(double)the_genomes->get_num_genomes()))*9.0;
                                
                                if (the_homologs->get_dupl_level() ==3) font_size=0.66*font_size;
                                
                                pl_ffontsize(font_size);
                            }
                            
                            text_width = pl_flabelwidth (outname.c_str());
                            
                            if (arrow_left==TRUE) {
                                //pl_fmove(j*box_x_foot+offset-(0.3*box_x_foot)-text_width, tracky[my_track]-taxa*box_y_foot+width/2.0 + width*(2.0/6.0));
                                pl_fmove(j*box_x_foot+offset-(0.5*spacer_x), tracky[my_track]-taxa*box_y_foot+1.1*width);
                                pl_ftextangle(0.0);
                                pl_alabel ('c', 'c',outname.c_str());
                                
                                pl_pencolorname("red");
                                pl_fline(j*box_x_foot+offset-(0.5*spacer_x), tracky[my_track]-taxa*box_y_foot+width/2.0, j*box_x_foot+offset-(0.5*spacer_x), tracky[my_track]-taxa*box_y_foot+width);
                                pl_fline(j*box_x_foot+offset-(0.5*spacer_x), tracky[my_track]-taxa*box_y_foot+width/2.0, j*box_x_foot+offset-(0.75*spacer_x), tracky[my_track]-taxa*box_y_foot+width/2.0+ width*(2.5/12.0));
                                pl_fline(j*box_x_foot+offset-(0.5*spacer_x), tracky[my_track]-taxa*box_y_foot+width/2.0, j*box_x_foot+offset-(0.25*spacer_x), tracky[my_track]-taxa*box_y_foot+width/2.0+ width*(2.5/12.0));
                                pl_pencolorname("black");
                                
                            }
                            else {
                                pl_fmove(j*box_x_foot+offset-(0.5*spacer_x), tracky[my_track]-taxa*box_y_foot+1.1*width);
                                pl_ftextangle(0.0);
                                pl_alabel ('c', 'c',outname.c_str());
                                
                                pl_pencolorname("red");
                                pl_fline(j*box_x_foot+boxsize+offset+(0.5*spacer_x), tracky[my_track]-taxa*box_y_foot+width/2.0, j*box_x_foot+boxsize+offset+(0.5*spacer_x), tracky[my_track]-taxa*box_y_foot+width);
                                pl_fline(j*box_x_foot+boxsize+offset+(0.5*spacer_x), tracky[my_track]-taxa*box_y_foot+width/2.0, j*box_x_foot+boxsize+offset+(0.75*spacer_x), tracky[my_track]-taxa*box_y_foot+width/2.0+ width*(2.5/12.0));
                                pl_fline(j*box_x_foot+boxsize+offset+(0.5*spacer_x), tracky[my_track]-taxa*box_y_foot+width/2.0, j*box_x_foot+boxsize+offset+(0.25*spacer_x), tracky[my_track]-taxa*box_y_foot+width/2.0+ width*(2.5/12.0));
                                pl_pencolorname("black");
                            }
                        }
                    }
                    
                    //pl_fbox (j*box_x_foot+offset, tracky[the_model->tracking_permutes[taxa_track_ids[j][taxa]][level]]-taxa*box_y_foot, j*box_x_foot+boxsize+offset, tracky[the_model->tracking_permutes[taxa_track_ids[j][taxa]][level]]-taxa*box_y_foot+width);
                    if (IPC_call==FALSE)
                        cout<<"Box "<<j<<" is "<<colors[color_ids[j]]<<endl;
                    pl_filltype(1);
                    pl_fillcolorname(colors[color_ids[j]].c_str());
                    
                    if (start+j == focus) {
                        pl_flinewidth(2.5);    /* set line thickness */
                        if (focus_name == (*the_model->the_tracks)[taxa].get_gene_track(start+j, level)->my_locus->get_gene_obj((*the_model->the_tracks)[taxa].get_gene_track(start+j, level)->index_num)->get_name_string())
                            pl_pencolorname("orchid");
                        else
                            pl_pencolorname("orchid4");
                    }
                    //pl_fillcolorname("yellow");
                    pl_fbox (j*box_x_foot+offset, tracky[my_track]-taxa*box_y_foot, j*box_x_foot+boxsize+offset, tracky[my_track]-taxa*box_y_foot+width);
                    if (start+j == focus) {
                        pl_flinewidth(1.0);    /* set line thickness */
                        pl_pencolorname("black");
                    }
                    
                    pl_filltype(0);
                    //figure_points[the_model->tracking_permutes[taxa_track_ids[j][taxa]][level]][j] = (*the_model->the_tracks)[taxa].get_gene_track(start+j, level);
                    
                    figure_points[my_track][j] = (*the_model->the_tracks)[taxa].get_gene_track(start+j, level);
                    //last_drawn[the_model->tracking_permutes[taxa_track_ids[j][taxa]][level]]=(*the_model->the_tracks)[taxa].get_gene_track(start+j, level);
                    last_drawn[my_track]=(*the_model->the_tracks)[taxa].get_gene_track(start+j, level);
                    if (IPC_call==FALSE)
                        cout<<"Extracting name from "<<start+j<<" level "<<level<<": "<<(*the_model->the_tracks)[taxa].get_gene_track(start+j, level)->my_locus->get_gene_obj((*the_model->the_tracks)[taxa].get_gene_track(start+j, level)->index_num)->get_name_string()<<" at "<<tracky[my_track]-taxa*box_y_foot<<" ie: "<<my_track<<endl;
                    fullname=(*the_model->the_tracks)[taxa].get_gene_track(start+j, level)->my_locus->get_gene_obj((*the_model->the_tracks)[taxa].get_gene_track(start+j, level)->index_num)->get_name_string();
                    
                    if ((prefix_map!=0) &&(prefix_map[taxa] != "")) {
                        name_head=fullname.substr(0, prefix_map[taxa].length());
                        cout<<"Checking header "<<name_head<<" against "<<prefix_map[taxa]<<endl;
                        if (name_head == prefix_map[taxa]) {
                            fullname=fullname.substr(prefix_map[taxa].length(), fullname.length());
                        }
                    }
                    
                    subname=extract_name(fullname);
                    
                    my_line_size=base_line_size;
                    //my_line_size=LINE_SIZE;
                    //if (the_genomes->get_num_genomes()>6) my_line_size -=2;
                    
                    //if (the_genomes->get_num_genomes()<5) my_line_size++;
                    
                    if (subname.length()==my_line_size+1) {
                        my_line_size++;
                        num_lines=1;
                    }
                    else {
                        num_lines = subname.length()/my_line_size;
                        
                        
                        
                        if (subname.length()%my_line_size !=0) num_lines++;
                        
                        if (num_lines > 3) num_lines=3;
                        if ((x_y_aspect>1.6) && num_lines==3) num_lines=2;
                    }
                    
                    if (vertical==TRUE) {
                        if (num_lines==1) center=boxsize*0.5;
                        else {
                            if (num_lines==2) center=boxsize*0.333;
                            if (num_lines==3) center=boxsize*0.15;
                        }
                    }
                    else {
                        if (num_lines==1) center=width*0.5;
                        else {
                            if (num_lines==2) center=width*0.666;
                            if (num_lines==3) center=width*0.85;
                        }
                    }
                    
                    for(name_line=0; name_line <num_lines; name_line++) {
                        split_name[name_line]=subname.substr((name_line*my_line_size), my_line_size);
                        //if (name_line != (num_lines-1)) split_name[name_line]+="-";
                    }

                    pl_ffontsize(name_font_size);
                    //if (bitmap==FALSE) {
                        //pl_ffontsize (7.5);
                    //}
                    //else {
                        //font_size=((5.0/(double)the_genomes->get_num_genomes()))*9.0;
                        
                        //if (the_homologs->get_dupl_level() ==3) font_size=0.66*font_size;
                        
                        //pl_ffontsize(font_size);
                   // }
                    
                    //text_div = (0.4*box_x_foot)/(double)(num_lines);
                    if (vertical==TRUE)
                        text_div = 0.7*name_font_size;
                    else
                        text_div = 0.79*name_font_size;
                    //text_div = font_size;
                    
                    my_offset = (0.5*box_x_foot)/ (double)num_lines;
                    
                    for(name_line =0; name_line < num_lines; name_line++) {
                        if (IPC_call==FALSE)
                            cout<<"Adding line "<<name_line<<" of "<<subname<<" ( total= "<<num_lines<<" which is "<<split_name[name_line]<<" x = "<<endl;
                        text_width = pl_flabelwidth (split_name[name_line].c_str());
                        //pl_fmove(j*box_x_foot+offset+0.55*text_width, tracky[the_model->tracking_permutes[taxa_track_ids[j][taxa]][level]]-taxa*box_y_foot+0.5*width);
                        //pl_fmove(j*box_x_foot+offset+0.5*box_x_foot, tracky[the_model->tracking_permutes[taxa_track_ids[j][taxa]][level]]-taxa*box_y_foot+0.5*width);
                        //pl_fmove(j*box_x_foot+offset+(name_line*text_div) + 0.2*box_x_foot, tracky[my_track]-taxa*box_y_foot+0.5*width);
                        
                        if (vertical==TRUE) {
                            cout<<"DIM is "<<box_dim<<" Placing line "<<name_line<<" at "<<j*box_x_foot+offset<<" + "<<center<<" + "<<name_line*text_div<<" to give : "<<split_name[name_line]<<endl;
                            //pl_fmove(j*box_x_foot+offset+(name_line*text_div) + my_offset, tracky[my_track]-taxa*box_y_foot+0.5*width);
                            pl_fmove(j*box_x_foot+offset+center+(name_line*text_div), tracky[my_track]-taxa*box_y_foot+0.5*width);
                            pl_ftextangle(90.0);
                            
                        }
                        else {
                            pl_fmove(j*box_x_foot+offset+(0.5*text_width)+(0.03*box_x_foot), tracky[my_track]-taxa*box_y_foot-(name_line*text_div)+center);
                            //pl_fmove(j*box_x_foot+offset+(name_line*text_div) + my_offset, tracky[my_track]-taxa*box_y_foot+0.5*width);
                            pl_ftextangle(0.0);
                            
                        }
                        pl_alabel ('c', 'c',split_name[name_line].c_str());
                    }
                    
                    //track_first[the_model->tracking_permutes[taxa_track_ids[j][taxa]][level]] = FALSE;
                    track_first[my_track] = FALSE;
                }
              
	      
                if ((*the_model->the_tracks)[taxa].get_gene_track(start+j, level)->my_locus != 0) {
                    if ((*the_model->the_tracks)[taxa].get_gene_track(start+j, level)->next != 0) {
                        //has_last[the_model->tracking_permutes[taxa_track_ids[j][taxa]][level]]=TRUE;
                        //lastline[the_model->tracking_permutes[taxa_track_ids[j][taxa]][level]]=j*box_x_foot+boxsize+offset;
                        has_last[my_track]=TRUE;
                        lastline[my_track]=j*box_x_foot+boxsize+offset;
                        
                    }
                    else
                        //has_last[the_model->tracking_permutes[taxa_track_ids[j][taxa]][level]]=FALSE;
                        has_last[my_track]=FALSE;
                    
                }
            }
            last_max_prob_pattern[j]=max_prob_pattern[j];
        }

        for(level=0; level<the_homologs->get_dupl_level(); level++) {
            if (last_drawn[level] !=0) {
                if (last_drawn[level]->next != 0 ) {
                    find_loc_pillar=-1;
                    for(k=0; k<num_genes; k++) {
                        for(other_level=0; other_level<the_homologs->get_dupl_level(); other_level++) {
                            if (last_drawn[level]->next == figure_points[other_level][k])
                                find_loc_pillar = k;
                        }
                    }
                    if (find_loc_pillar == -1)
                        pl_fline(lastline[level], tracky[level]-taxa*box_y_foot+width/2.0,
                                 (num_genes-1)*box_x_foot+boxsize+offset+offset*0.1,tracky[level]-taxa*box_y_foot+width/2.0);
              
                }
            }
        }
	  
    }


  
    delete[] max_prob_pattern;
    delete[] max_prob;
    delete[] last_max_prob_pattern;
    delete[] track_first;
    delete[] tracky;
    delete[] color_ids;
    for (level=0; level<the_homologs->get_dupl_level(); level++)
        delete[] figure_points[level];
    delete[] figure_points;
    delete[] lastline;
    
    for(i=0; i<num_per_page; i++) delete[] taxa_track_ids[i];
    delete[] taxa_track_ids;
	
   
    //fflush(outfile);
  if (pl_closepl () < 0)      /* close Plotter */ {
      fprintf (stderr, "Couldn't close Plotter\n");
      return 1;
    }
#if 1
    pl_selectpl (0);            /* select default Plotter */
    if (pl_deletepl (thandle) < 0) /* delete Plotter we used */ {
        fprintf (stderr, "Couldn't delete Plotter\n");
        return 1;
    }
    if (last_drawn!=0)
        delete[] last_drawn;
#endif

    if (IPC_call==TRUE) {
        std::cout<<"Extracting from buffer\n";
        fflush(outfile);
        
        //printf ("At close buf = `%s', size = %zu\n", buf_pt, size);
        mypos=(size_t)0;
        
        //std::cout<<0<<": "<<buf_pt[0]<<std::endl<<std::flush;
        //std::cout<<1<<": "<<buf_pt[1]<<std::endl<<std::flush;
        //std::cout<<2<<": "<<buf_pt[2]<<std::endl<<std::flush;
        //std::cout<<mypos<<": "<<buf_pt[mypos]<<std::endl<<std::flush;
        
        while(mypos < size) {
            plot_ss->put((char)buf_pt[mypos]);
            //cout<<"P: "<<mypos<<" = "<<(int)buf_pt[mypos]<<endl;
            ++mypos;
        }
        //plot_ss->put('\0');
        //while ((readchar = fgetc (outfile)) != EOF)
        //for (std::size_t cntm = 0; cntm < size; cntm++) {
         //   std::cout<<cntm<<": "<<buf_pt[cntm]<<std::endl<<std::flush;
        //    *plot_ss<<buf_pt[cntm];
        //}
        
        //printf ("Got %c\n", ch);
        //fclose (outfile);
        //free(buf_pt);
    }
#if 0
  pl_selectpl (0);            /* select default Plotter */
  if (pl_deletepl (thandle) < 0) /* delete Plotter we used */ {
      fprintf (stderr, "Couldn't delete Plotter\n");
      return 1;
    }
    if (last_drawn!=0)
        delete[] last_drawn;
#endif
    fclose(outfile);

 return(0);
}

#if 0
void get_max_prob_pattern(Exchange *curr_exchange, Phylo_Matrix *the_matrix, Ploidy_Like_model *the_model, int site, int &opt_track_id, int *taxa_track_ids, double &max_prob)
{
    int  track_id;
    
   
    
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
        
        the_model->set_track_states(opt_track_id, taxa_track_ids);
    }
		
}
#endif

string extract_name(string input_name)
{
    string ret_val, match, end_string;
    std::size_t found, found2;
    
    ret_val=input_name;
    
    found=ret_val.find("AT");
    
    if (found!=std::string::npos) {
        if (found == 0) {
            ret_val=ret_val.substr(0, ret_val.length()-2);
        }
    }
    
    
    found=ret_val.find("__AT");
    
    if (found!=std::string::npos) {
        found2=ret_val.find("_scaffold");
        if (found2!=std::string::npos) {
            ret_val=ret_val.substr(found+2, found2-(found+2)-2);
        }
    }
    
    
    
    found=ret_val.find("PAC:");
    
    if (found!=std::string::npos) {
        if (found == 0) {
            ret_val=ret_val.substr(4, ret_val.length()-4);
        }
    }
    
    found=ret_val.find("_scaffold_");
    
    if (found!=std::string::npos) {
        if (found >= ret_val.length()-13) {
            ret_val=ret_val.substr(0, found);
        }
    }
    
    found=ret_val.find("_ch");
    
    if (found!=std::string::npos) {
        if (found >= ret_val.length()-7) {
            ret_val=ret_val.substr(0,  found);
        }
    }
    
    
    found=ret_val.find("_Scaffold_");
    
    if (found!=std::string::npos) {
        if (found >= ret_val.length()-13) {
            ret_val=ret_val.substr(0,  found);
        }
    }


    if (ret_val.length()>25) {
        ret_val=ret_val.substr(0,15);
    }
    
    match="_gene";
    
    if (ret_val.length()>5) {
        end_string=ret_val.substr(ret_val.length()-5,ret_val.length());
        
        if (end_string == match) {
            ret_val=ret_val.substr(0, ret_val.length()-5);
        }
    }
    
    return(ret_val);
}

void make_prefixes (Clade *the_genomes, std::string *prefixes, std::string *&prefix_map)
{
    int i, taxa;
    string spp_name, prefix;
    
    prefix_map=0;
    
    if (prefixes !=0) {
        prefix_map=new string [the_genomes->get_num_genomes()];
        for(i=0; i<the_genomes->get_num_genomes(); i++) prefix_map[i]="";
        
        for(i=0; i<the_genomes->get_num_genomes(); i++) {
            //cout<<"HAve prefix entry "<<i<<" of "<<prefixes[i]<<endl;
            
            std::size_t found = prefixes[i].find("=");
            if (found!=std::string::npos) {
                spp_name = prefixes[i].substr(0, found);
                prefix = prefixes[i].substr(found+1, prefixes[i].length());
                
                cout<<"Extracted spp "<<spp_name<<" and pre="<<prefix<<" from "<<prefixes[i]<<endl;
                
                taxa=0;
                
                while ((spp_name != (*the_genomes)[taxa].get_name()) && (taxa <the_genomes->get_num_genomes() )) taxa++;
                    
                if (taxa <the_genomes->get_num_genomes()) {
                    cout<<"Stripping prefix "<<prefix<<" from "<<(*the_genomes)[taxa].get_name()<<endl;
                    prefix_map[taxa]=prefix;
                }
            }
        }
    }
}

void construct_full_genome_lookup(string genomefile, Clade *the_genomes, std::map<std::string, int> gene_hash, std::map<std::string, int> &left_matches, std::map<std::string, int> &right_matches, std::map<std::string, int> &taxa_lookup, std::map<std::string, std::string> &left_name, std::map<std::string, std::string> &right_name, std::map<std::string, std::string> &aliases)
{
    int i, taxa_num=0, site, chrom_num, chrom, found_site, last_chrom=-1, curr_left, curr_right, taxa;
    string spp_name, gene_name, line, alias, alias_uc, alias_lc;
    ifstream fin;
    std::map<int, int> chrom_lookup;
    
   
    
    fin.open(genomefile.c_str());
    
    if (fin.fail()) {
        std::cerr<<"Error: Could not open full genome order file "<<genomefile<<std::endl;
        return;
    }
    fin>>spp_name;

    while (!(spp_name == (*the_genomes)[taxa_num].get_name()) && (taxa_num<the_genomes->get_num_genomes()))
        taxa_num++;
    
    if (taxa_num == the_genomes->get_num_genomes()){
        std::cerr<<"Error: genome  "<<spp_name<<" not found in the genomes"<<std::endl;
        return;
    }
    
    while (!(fin.eof())) {
        //fin>>chrom_num>>gene_name;
        std::getline(fin, line);
        std::stringstream ss (line);
        
        ss>>chrom_num>>gene_name;
        
        aliases[gene_name]=gene_name;
        
        while(!ss.eof()) {
            alias="NA";
            ss>>alias;
            
            if (alias != "NA") {
                alias_lc=alias;
                alias_uc=alias;
                std::transform(alias_lc.begin(), alias_lc.end(),alias_lc.begin(), ::tolower);
                std::transform(alias_uc.begin(), alias_uc.end(),alias_uc.begin(), ::toupper);
                
                if ((aliases.find(alias_uc) == aliases.end()) &&(aliases.find(alias_lc) == aliases.end())) {
                    aliases[alias_lc]=gene_name;
                    aliases[alias_uc]=gene_name;
                }
                else {
                    aliases[alias_lc]="NONE";
                    aliases[alias_uc]="NONE";
                }
            }
        }
        
        if (chrom_lookup.find(chrom_num) == chrom_lookup.end()) {
            
            for(chrom=0; chrom<(*the_genomes)[taxa_num].get_num_contigs(); chrom++ ) {
                found_site=-1;
                for(site=0; site<(*the_genomes)[taxa_num][chrom].get_num_genes(); site++) {
                    if ((*the_genomes)[taxa_num][chrom][site].get_name_string() == gene_name)
                        found_site=site;
                }
                if (found_site != -1) chrom_lookup[chrom_num]=chrom;
            }
        }
    }
    
    fin.close();
    fin.clear();
    fin.open(genomefile.c_str());

    fin>>spp_name;

    while (!(fin.eof())) {
        //fin>>chrom_num>>gene_name;
        
        std::getline(fin, line);
        std::stringstream ss (line);
        
        ss>>chrom_num>>gene_name;
        
        
        while(!ss.eof())
            ss>>alias;
        
        
        
        if (chrom_lookup.find(chrom_num) == chrom_lookup.end()) {
            left_matches[gene_name]=-1;
            right_matches[gene_name]=-1;
            taxa_lookup[gene_name]=-1;
        }
        else {
            if (chrom_num != last_chrom) {
                curr_left=-1;
                curr_right=0;
                last_chrom=chrom_num;
            }
            
            if ((*the_genomes)[taxa_num][chrom_lookup[chrom_num]][curr_right].get_name_string() == gene_name) {
                curr_left=curr_right;
                curr_right++;
                
                if (curr_right == (*the_genomes)[taxa_num][chrom_lookup[chrom_num]].get_num_genes()) curr_right=-1;
            }
            else {
                taxa_lookup[gene_name]=taxa_num;
                if (curr_left!=-1) {
                    left_matches[gene_name]=gene_hash[(*the_genomes)[taxa_num][chrom_lookup[chrom_num]][curr_left].get_name_string()];
                    left_name[gene_name]=(*the_genomes)[taxa_num][chrom_lookup[chrom_num]][curr_left].get_name_string();
                }
                else
                    left_matches[gene_name]=-1;
                if (curr_right != -1) {
                    right_matches[gene_name]=gene_hash[(*the_genomes)[taxa_num][chrom_lookup[chrom_num]][curr_right].get_name_string()];
                    right_name[gene_name]=(*the_genomes)[taxa_num][chrom_lookup[chrom_num]][curr_right].get_name_string();
                }
                
                else
                    right_matches[gene_name]=-1;
            }
            
            
        }
    }
    fin.close();
    

}


int draw_model_diag(Phylo_Matrix *the_matrix, string plotfile, BOOL bitmap, BOOL IPC_call, std::stringstream *plot_ss, int matrix_num)
{
    int level, state_num, state_num2, thandle, realy, realx, max_level, param_cnt, my_level, max_state, *state_colors, color_cnt=0;
    double spacex, spacey, border=10, spacer, name_font_size=1.0, center, my_width, text_x_offset, text_y_offset,
      xcenter, ycenter, label_size, width, boxsize, text_div, font_size, old_fontsize, circle_size, level_size, slope,
        labelx, labely, text_width, x_start, x_end, y_start, y_end, angle, delta_x, delta_y, arrow_len, new_dx, new_dy;
    char write_string[100],  attrib[100], value[100], size_string[25], readchar;
    BOOL same_level, draw_line;
    string statename, state_head, colors[8], line_colors[4], param_string;
    char* buf_pt;
    FILE *outfile;
    std::map<string, int> color_lookup, state_xcoords, state_ycoords, state_order;
    std::size_t found, size, mypos;
    
    colors[0]="mistyrose";
    colors[1]="lightsteelblue1";
    colors[2]="springgreen3";
    colors[3]="bisque3";
    colors[4]="lightsalmon2";
    colors[5]="lavender";
    colors[6]="palegoldenrod";
    colors[7]="wheat1";
    
    line_colors[0]="thistle4";
    line_colors[1]="midnightblue";
    line_colors[2]="purple";
    line_colors[3]="red";
    
    state_colors = new int [the_matrix->get_num_states()];
    
    for(state_num=0; state_num<the_matrix->get_num_states(); state_num++) state_colors[state_num] =-1;
    
    the_matrix->initialize_rate_matrix(matrix_num);
    
    max_level=1;
    for (level=0; level<=the_matrix->get_num_levels(); level++) {
        if (the_matrix->get_num_full_level_states(level)>max_level) max_level=the_matrix->get_num_full_level_states(level);
        for(state_num=0; state_num<the_matrix->get_num_full_level_states(level); state_num++) {
            state_order[the_matrix->get_ith_full_state_level(level, state_num)->get_state_name()]=state_num;
            if (level == 3)
                state_head=the_matrix->get_ith_full_state_level(level, state_num)->get_state_name().substr(0, 5);
            else {
                if ((level == 2) || (level == 4))
                    state_head=the_matrix->get_ith_full_state_level(level, state_num)->get_state_name().substr(0, 4);
                else
                    state_head=the_matrix->get_ith_full_state_level(level, state_num)->get_state_name();
            }
        
        
            cout<<"Looking at level "<<level<<" sub state "<<state_num<<endl;
        
            color_lookup[the_matrix->get_ith_full_state_level(level, state_num)->get_state_name()]=7;
            if (state_head == "Tripl") color_lookup[the_matrix->get_ith_full_state_level(level, state_num)->get_state_name()]=5;
            if (state_head == "Dupl") color_lookup[the_matrix->get_ith_full_state_level(level, state_num)->get_state_name()]=0;
            if (state_head == "Copy1") color_lookup[the_matrix->get_ith_full_state_level(level, state_num)->get_state_name()]=1;
            if (state_head == "Copy2") color_lookup[the_matrix->get_ith_full_state_level(level, state_num)->get_state_name()]=2;
            if (state_head == "Copy3") color_lookup[the_matrix->get_ith_full_state_level(level, state_num)->get_state_name()]=4;
        }
    }

    
   
    if (IPC_call==FALSE)
        outfile=fopen(plotfile.c_str(), "w");
    else {
        size=0;
        //printf ("At open buf = `%s', size = %zu\n", buf_pt, size);
        outfile= open_memstream (&buf_pt, &size);
        
        //fprintf(outfile, "TESTING");
        fflush(outfile);
        //printf ("At test buf = `%s', size = %zu\n", buf_pt, size);
    }
   
    spacex=500.0;
    spacey=500.0;
   
    level_size= spacey/(1.0*the_matrix->get_num_levels());
    
    //offset=offset_y=0;
    
    if ((the_matrix->get_num_levels())>max_level)
        //circle_size = (0.8*(spacey/(1.0*(the_matrix->get_num_levels()-1))))/2.0;
        circle_size=0.65*level_size;
    else
        circle_size=(0.65*(spacex/(1.0*(max_level))))/2.0;

    arrow_len=0.3*circle_size;
    
    cout<<"Circle Radius is "<<circle_size<<endl;
    if (bitmap == FALSE) {
        strcpy(attrib, "PAGESIZE");
        strcpy(value, "letter,xsize=5in,ysize=5in");
        pl_parampl (attrib, value);
        /* create a Postscript Plotter that writes to standard output */
        
        if ((thandle = pl_newpl ("ps", stdin, outfile, stderr)) < 0) {
              fprintf (stderr, "Couldn't create Plotter\n");
              return 1;
        }
       
        pl_selectpl (thandle);       /* select the Plotter for use */
    }
    else {
        
        strcpy(size_string, "1000x1000");
        realx=1000;
        realy=1000;
                
        /* set a Plotter parameter */
        pl_parampl ("BITMAPSIZE", size_string);
       
        if ((thandle = pl_newpl ("png", stdin, outfile, stderr)) < 0) {
            fprintf (stderr, "Couldn't create Plotter\n");
            return 1;
        }
       
        pl_selectpl (thandle);       /* select the Plotter for use */
    }

    
    
    if (pl_openpl () < 0)       /* open Plotter */ {
      fprintf (stderr, "Couldn't open Plotter\n");
      return 1;
    }
    
    pl_fspace (0.0, 0.0, spacex, spacey); /* specify user coor system */
    pl_flinewidth (1.0);       /* line thickness in user coordinates */
    pl_pencolorname ("black");    /* path will be drawn in red */
    pl_erase ();                /* erase Plotter's graphics display */
    /*pl_fmove (600.0, 300.0);*/    /* position the graphics cursor */
    if (bitmap == TRUE) {
        pl_fontname ("HersheySans-BoldOblique"); /* choose a Postscript font */
        pl_ftextangle (0);         /* text inclination angle (degrees) */
        pl_ffontsize (17);
        font_size=17.0;
    }
    else {
        pl_fontname("Helvetica-Oblique");
        pl_ftextangle (0);         /* text inclination angle (degrees) */
        pl_ffontsize (17);
        font_size=17;
    }
    
    cout<<"Model has "<<the_matrix->get_num_levels()<<" levels"<<endl;
    for (level=0; level<=the_matrix->get_num_levels(); level++) {
       // ycenter = spacey*((1.0*(level-1)/(1.0*(the_matrix->get_num_levels()))) - 0.5*(1.0/(1.0*(the_matrix->get_num_levels()))) );
        ycenter=level_size*level - 0.5*level_size;
        cout<<"At level "<<level<<". Substates: "<<the_matrix->get_num_full_level_states(level)<<endl;
        for(state_num=0; state_num<the_matrix->get_num_full_level_states(level); state_num++) {
            xcenter = spacex* ((1.0*(state_num+1))/(1.0*the_matrix->get_num_full_level_states(level)) - 0.5*(1.0/(1.0*the_matrix->get_num_full_level_states(level))) );
            
            state_xcoords[the_matrix->get_ith_full_state_level(level, state_num)->get_state_name()]=xcenter;
            state_ycoords[the_matrix->get_ith_full_state_level(level, state_num)->get_state_name()]=ycenter;
            
            cout<<"State "<<the_matrix->get_ith_full_state_level(level, state_num)->get_state_name()<<" x: "<<xcenter<<" y: "<<ycenter<<" Color: "
            <<colors[color_lookup[the_matrix->get_ith_full_state_level(level, state_num)->get_state_name()]]<<endl;
            
            pl_filltype(1);
            pl_fillcolorname(colors[color_lookup[the_matrix->get_ith_full_state_level(level, state_num)->get_state_name()]].c_str());
            pl_fcircle (xcenter, ycenter, circle_size);
            
            pl_fmove(xcenter, ycenter);
            pl_alabel('c', 'c', the_matrix->get_ith_full_state_level(level, state_num)->get_state_name().c_str());
        }
    }
    
    pl_ffontsize (10.75);
    pl_filltype(0);
    for(state_num=0; state_num<the_matrix->get_num_states(); state_num++) {
        for(state_num2=state_num+1; state_num2<the_matrix->get_num_states(); state_num2++) {
            
            if (the_matrix->rate_matrix[state_num][state_num2] != 0) {
                my_width=4.0*the_matrix->rate_matrix[state_num][state_num2];
                pl_flinewidth (my_width);
                
                cout<<the_matrix->get_nth_state(state_num)->get_state_name()<<" to "<<the_matrix->get_nth_state(state_num2)->get_state_name()<<" Width is "<<my_width<<endl;
                
                if (the_matrix->get_nth_state(state_num)->get_state_level()>
                    the_matrix->get_nth_state(state_num2)->get_state_level()) {
                    max_state=state_num;
                }
                else{
                    if (the_matrix->get_nth_state(state_num)->get_state_level() ==
                        the_matrix->get_nth_state(state_num2)->get_state_level()) {
                        
                        if (state_num > state_num2) max_state=state_num2;
                        else max_state=state_num;
                    }
                    else
                    max_state=state_num2;
                }
                
                if (state_colors[max_state] == -1) {
                    state_colors[max_state]=color_cnt;
                    color_cnt++;
                }
                
                draw_line=TRUE;
                
                if (the_matrix->get_nth_state(state_num)->get_state_level() <  the_matrix->get_nth_state(state_num2)->get_state_level()) my_level=the_matrix->get_nth_state(state_num)->get_state_level();
                else
                    my_level=the_matrix->get_nth_state(state_num2)->get_state_level();
                
                if (the_matrix->get_nth_state(state_num)->get_state_level() !=  the_matrix->get_nth_state(state_num2)->get_state_level()) {
                    same_level=FALSE;
                    if (state_xcoords[the_matrix->get_nth_state(state_num)->get_state_name()]>
                        state_xcoords[the_matrix->get_nth_state(state_num2)->get_state_name()]) {
                        //labelx = state_xcoords[the_matrix->get_nth_state(state_num)->get_state_name()] -state_xcoords[the_matrix->get_nth_state(state_num2)->get_state_name()];
                        //labelx = labelx/2.0;
                        //labelx +=state_xcoords[the_matrix->get_nth_state(state_num2)->get_state_name()];
                    }
                    else {
                       // labelx = state_xcoords[the_matrix->get_nth_state(state_num2)->get_state_name()] -state_xcoords[the_matrix->get_nth_state(state_num)->get_state_name()];
                        //labelx = labelx/2.0;
                        //labelx +=state_xcoords[the_matrix->get_nth_state(state_num)->get_state_name()];
                    }
                    if (state_ycoords[the_matrix->get_nth_state(state_num)->get_state_name()]>
                        state_ycoords[the_matrix->get_nth_state(state_num2)->get_state_name()]) {
                        
                        x_start=state_xcoords[the_matrix->get_nth_state(state_num)->get_state_name()];
                        y_start=state_ycoords[the_matrix->get_nth_state(state_num)->get_state_name()]-circle_size;
                        
                        x_end=state_xcoords[the_matrix->get_nth_state(state_num2)->get_state_name()];
                        y_end=state_ycoords[the_matrix->get_nth_state(state_num2)->get_state_name()]+circle_size;
                        
                        
                        //labely=state_ycoords[the_matrix->get_nth_state(state_num)->get_state_name()]-
                        //state_ycoords[the_matrix->get_nth_state(state_num2)->get_state_name()];
                        //labely = labely/2.0;
                        //labely+=state_ycoords[the_matrix->get_nth_state(state_num2)->get_state_name()];
                        
                      
                    }
                    else {
                        x_start=state_xcoords[the_matrix->get_nth_state(state_num2)->get_state_name()];
                        y_start=state_ycoords[the_matrix->get_nth_state(state_num2)->get_state_name()]+circle_size;
                        
                        x_end=state_xcoords[the_matrix->get_nth_state(state_num)->get_state_name()];
                        y_end=state_ycoords[the_matrix->get_nth_state(state_num)->get_state_name()]-circle_size;
                        
                        
                        //labely=state_ycoords[the_matrix->get_nth_state(state_num2)->get_state_name()]-
                        state_ycoords[the_matrix->get_nth_state(state_num)->get_state_name()];
                        //labely = labely/2.0;
                        //labely+=state_ycoords[the_matrix->get_nth_state(state_num)->get_state_name()];
                        
                        
                    }
                    
                    
                }
                else {
                    same_level=TRUE;
                    if (abs(state_order[the_matrix->get_nth_state(state_num2)->get_state_name()]-
                            state_order[the_matrix->get_nth_state(state_num)->get_state_name()])==1) {
                        labely=state_ycoords[the_matrix->get_nth_state(state_num2)->get_state_name()]+22;
                        
                        
                        
                        if (state_xcoords[the_matrix->get_nth_state(state_num)->get_state_name()]>
                            state_xcoords[the_matrix->get_nth_state(state_num2)->get_state_name()]) {
                            labelx = ((state_xcoords[the_matrix->get_nth_state(state_num)->get_state_name()] -state_xcoords[the_matrix->get_nth_state(state_num2)->get_state_name()])/2.0)+
                            state_xcoords[the_matrix->get_nth_state(state_num2)->get_state_name()];
                            
                            
                            x_start=state_xcoords[the_matrix->get_nth_state(state_num2)->get_state_name()]+circle_size;
                            y_start=state_ycoords[the_matrix->get_nth_state(state_num2)->get_state_name()];
                                                  
                            x_end=state_xcoords[the_matrix->get_nth_state(state_num)->get_state_name()]-circle_size;
                            y_end=state_ycoords[the_matrix->get_nth_state(state_num)->get_state_name()];
                                           
                            //pl_fline(state_xcoords[the_matrix->get_nth_state(state_num)->get_state_name()]-circle_size, state_ycoords[the_matrix->get_nth_state(state_num)->get_state_name()], state_xcoords[the_matrix->get_nth_state(state_num2)->get_state_name()]+circle_size, state_ycoords[the_matrix->get_nth_state(state_num2)->get_state_name()]);
                            
                        }
                        else {
                            labelx = ((state_xcoords[the_matrix->get_nth_state(state_num2)->get_state_name()] -state_xcoords[the_matrix->get_nth_state(state_num)->get_state_name()])/2.0)+
                            state_xcoords[the_matrix->get_nth_state(state_num)->get_state_name()];
                        
                            x_start=state_xcoords[the_matrix->get_nth_state(state_num)->get_state_name()]+circle_size;
                            y_start=state_ycoords[the_matrix->get_nth_state(state_num)->get_state_name()];
                                                  
                            x_end=state_xcoords[the_matrix->get_nth_state(state_num2)->get_state_name()]-circle_size;
                            y_end=state_ycoords[the_matrix->get_nth_state(state_num2)->get_state_name()];
                            
                            
                            //pl_fline(state_xcoords[the_matrix->get_nth_state(state_num)->get_state_name()]+circle_size, state_ycoords[the_matrix->get_nth_state(state_num)->get_state_name()], state_xcoords[the_matrix->get_nth_state(state_num2)->get_state_name()]-circle_size, state_ycoords[the_matrix->get_nth_state(state_num2)->get_state_name()]);
                        }
                    }
                    else {
                        //labely=state_ycoords[the_matrix->get_nth_state(state_num2)->get_state_name()]+2.1*circle_size;
                        draw_line=FALSE;
                        delta_x=abs(state_xcoords[the_matrix->get_nth_state(state_num)->get_state_name()]-
                                    state_xcoords[the_matrix->get_nth_state(state_num2)->get_state_name()]);
                        if (state_xcoords[the_matrix->get_nth_state(state_num)->get_state_name()]>
                            state_xcoords[the_matrix->get_nth_state(state_num2)->get_state_name()]) {
                            x_start=state_xcoords[the_matrix->get_nth_state(state_num2)->get_state_name()];
                            x_end =state_xcoords[the_matrix->get_nth_state(state_num)->get_state_name()];
                            //labelx = 1.5*((state_xcoords[the_matrix->get_nth_state(state_num)->get_state_name()] -state_xcoords[the_matrix->get_nth_state(state_num2)->get_state_name()])/2.0)+
                            //state_xcoords[the_matrix->get_nth_state(state_num2)->get_state_name()];
                            
                        }
                        else {
                            x_start=state_xcoords[the_matrix->get_nth_state(state_num)->get_state_name()];
                            x_end =state_xcoords[the_matrix->get_nth_state(state_num2)->get_state_name()];
                            //labelx = 1.5*((state_xcoords[the_matrix->get_nth_state(state_num2)->get_state_name()] -state_xcoords[the_matrix->get_nth_state(state_num)->get_state_name()])/2.0)+
                            //state_xcoords[the_matrix->get_nth_state(state_num)->get_state_name()];
                        }
                        
                        if (state_num2%2 ==0) {
                            delta_y=1.5*circle_size;
                            labelx = 0.75*((x_end-x_start)/2.0)+x_start;
                        }
                        else {
                            labelx = 1.5*((x_end-x_start)/2.0)+x_start;
                            delta_y=circle_size;
                        }
                        
                        labely=state_ycoords[the_matrix->get_nth_state(state_num2)->get_state_name()]+0.1*delta_x+delta_y;
                        
                        pl_pencolorname (line_colors[state_colors[max_state]].c_str());
                        pl_fbezier2 (state_xcoords[the_matrix->get_nth_state(state_num)->get_state_name()],
                                     state_ycoords[the_matrix->get_nth_state(state_num)->get_state_name()]+circle_size,
                                     //labelx,
                                     0.5*delta_x+x_start,
                                     //state_ycoords[the_matrix->get_nth_state(state_num2)->get_state_name()]+(2.0*circle_size),
                                     state_ycoords[the_matrix->get_nth_state(state_num2)->get_state_name()]+0.1*delta_x+delta_y,
                                     state_xcoords[the_matrix->get_nth_state(state_num2)->get_state_name()], state_ycoords[the_matrix->get_nth_state(state_num2)->get_state_name()]+circle_size);
                        x_end=state_xcoords[the_matrix->get_nth_state(state_num2)->get_state_name()];
                        y_end=state_ycoords[the_matrix->get_nth_state(state_num2)->get_state_name()]+circle_size;
                        //pl_farc(labelx,state_ycoords[the_matrix->get_nth_state(state_num)->get_state_name()], state_xcoords[the_matrix->get_nth_state(state_num)->get_state_name()], state_ycoords[the_matrix->get_nth_state(state_num)->get_state_name()]+circle_size, state_xcoords[the_matrix->get_nth_state(state_num2)->get_state_name()], state_ycoords[the_matrix->get_nth_state(state_num2)->get_state_name()]+circle_size);
                        
                    }
                }
            }
            
            pl_pencolorname (line_colors[state_colors[max_state]].c_str());
            if (draw_line == TRUE) {
                
                if (my_width<1e-3) {
                    pl_linemod ("dotted");
                }
                else
                    pl_linemod ("solid");
                
                pl_fline(x_start, y_start, x_end, y_end);
                delta_x=x_end-x_start;
                delta_y=y_end-y_start;
                slope = delta_y/delta_x;
                angle = atan(delta_y/delta_x);
                
                if (same_level==FALSE) {
                    if (delta_x !=0) {
                        if (slope > 0) {
                            cout<<"High?\n";
                            labelx = 0.57*delta_x+x_start;
                            labely = slope*0.57*delta_x+y_start;
                            
                        }
                        if (slope < 0) {
                            cout<<"Low?\n";
                            labelx = 0.25*delta_x+x_start;
                            labely = slope*0.25*delta_x+y_start;
                            
                        }
                    }
                    else {
                        if (the_matrix->get_nth_state(state_num) == the_matrix->get_ith_full_state_level(the_matrix->get_nth_state(state_num)->get_state_level(), 0))
                        {
                            labelx = x_start;
                            labely = 0.57*delta_y+y_start;
                        }
                        else {
                            labelx = x_start;
                            labely = 0.25*delta_y+y_start;
                        }
                    }
                    cout<<"Labeling at "<<labelx<<", "<<labely<<" (slope="<<slope<<", delta_y="<<delta_y<<")"<<endl;
                    
                    
                }
                
            }
            else
                angle=2.75;
            
            pl_linemod ("solid");
            pl_flinewidth (2.0);
            
           
            new_dy=arrow_len*sin(angle+0.4);
            new_dx=arrow_len*cos(angle+0.4);
            if (draw_line == TRUE) {
                if (delta_x >=0)
                    pl_fline(x_end, y_end, x_end-new_dx, y_end-new_dy);
                else
                    pl_fline(x_end, y_end, x_end+new_dx, y_end+new_dy);
            }
            else {
                pl_fline(x_end, y_end, x_end+new_dx, y_end+new_dy);
            }
            
            cout<<"Line has "<<x_start<<", "<<y_start<<" to "<<x_end<<", "<<y_end<<" with DL "<<draw_line<<" ad " <<delta_x<<", "<<delta_y<<" angle: "<<angle<<" Arrow dx dy: "<<new_dx<<", "<<new_dy<<endl;
            
            new_dy=arrow_len*sin(angle-0.4);
            new_dx=arrow_len*cos(angle-0.4);
            if (draw_line == TRUE) {
                if (delta_x >=0)
                    pl_fline(x_end, y_end, x_end-new_dx, y_end-new_dy);
                else
                    pl_fline(x_end, y_end, x_end+new_dx, y_end+new_dy);
            }
            else {
                pl_fline(x_end, y_end, x_end+new_dx, y_end+new_dy);
            }
            
            pl_pencolorname ("black");
            param_string="";
            
            if (the_matrix->get_num_levels() == 3) {
                if (delta_x >0) labely+=10;
                else labely-=10;
            }
            
            
            
            if (the_matrix->get_trans_type(state_num,state_num2)==PARAM) {
                for(param_cnt=0; param_cnt< the_matrix->num_trans_params(state_num, state_num2); param_cnt++) {
                    param_string += the_matrix->get_nth_trans_param(state_num, state_num2, param_cnt)->get_global_param_name();
                    param_string += "*";
                }
                param_string = param_string.substr(0,param_string.length()-1);
                param_string += "=";
                text_width = pl_flabelwidth (param_string.c_str());
                
                cout<<"Printing "<<param_string<<" at "<<labelx<<" + "<<text_width<<endl;
                
                
                text_x_offset=0.0;
                text_y_offset=8;
                if (same_level==FALSE) {
                    text_y_offset=0;
                    //pl_fmove(labelx+0.8*text_width, labely);
                    if (delta_x ==0) {
                        if (the_matrix->get_nth_state(state_num2) == the_matrix->get_ith_full_state_level(the_matrix->get_nth_state(state_num2)->get_state_level(), 0))
                            text_x_offset=-0.57*text_width;
                        else
                            text_x_offset=0.57*text_width;
                    }
                    else {
                        if (slope <0)
                            text_x_offset=((0.0015*delta_x)*(1.0/abs(slope)))*text_width;
                        else
                            text_x_offset=0.25*text_width;
                    }
                }
              
                pl_fmove(labelx+text_x_offset, labely+text_y_offset);
                
                pl_alabel('c', 'c', param_string.c_str());
                
                param_string="";
                for(param_cnt=0; param_cnt< the_matrix->num_trans_params(state_num, state_num2); param_cnt++) {
                    std::ostringstream ss;
                    ss<< the_matrix->get_nth_trans_param(state_num, state_num2, param_cnt)->get_parameter(matrix_num);
                    param_string += ss.str();
                    param_string += "*";
                }
                param_string = param_string.substr(0,param_string.length()-1);
                
                
                //text_width = pl_flabelwidth (param_string.c_str());
                text_x_offset=0.0;
                text_y_offset=-5;
                if (same_level ==FALSE) {
                    text_y_offset=-12;
                        if (delta_x ==0) {
                            if (the_matrix->get_nth_state(state_num2) == the_matrix->get_ith_full_state_level(the_matrix->get_nth_state(state_num2)->get_state_level(), 0))
                                text_x_offset=-0.57*text_width;
                            else
                                text_x_offset=0.57*text_width;
                        }
                        else {
                            if (slope >0)
                                text_x_offset=0.2*text_width;
                            else
                                text_x_offset=((0.002*delta_x)*(1.0/abs(slope)))*text_width;
                        }
                }
               
                    
                pl_fmove(labelx+text_x_offset, labely+text_y_offset);
                pl_alabel('c', 'c', param_string.c_str());
            
            }
            
            
            
        }
    }
    
    delete[] state_colors;
    
    if (pl_closepl () < 0)      /* close Plotter */ {
        fprintf (stderr, "Couldn't close Plotter\n");
        return 1;
      }

      pl_selectpl (0);            /* select default Plotter */
      if (pl_deletepl (thandle) < 0) /* delete Plotter we used */ {
          fprintf (stderr, "Couldn't delete Plotter\n");
          return 1;
      }
      
      if (IPC_call==TRUE) {
          std::cout<<"Extracting from buffer\n";
          fflush(outfile);
          
          //printf ("At close buf = `%s', size = %zu\n", buf_pt, size);
          mypos=(size_t)0;
          
         
          while(mypos < size) {
              plot_ss->put((char)buf_pt[mypos]);
              //cout<<"P: "<<mypos<<" = "<<(int)buf_pt[mypos]<<endl;
              ++mypos;
          }
          
      }
 
      fclose(outfile);

   return(0);
}
