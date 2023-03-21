#include <sys/types.h>
#include <sys/socket.h>
#include <sys/un.h>
#include <unistd.h>
#include <stdio.h>
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
#include <random>

#define BUFFER_SIZE 2048
#define MAX_FRAME 150;
#define NUM_FIELDS 5

using namespace::std;

enum MSG_TYPE {MESSAGE, DATABLOCK, SIZE, PING, ENDCOMM, PHYLOTREE, BATCHBLOCK, MODELDIAG, TREEDIAG, GENOMELIST, INFOBLOCK};
enum IMAGE_SIZE {SMALL, MEDIUM, BIGGER, LARGE, MASSIVE};

extern int draw_region (Exchange *curr_exchange,  Tree *the_tree, Phylo_Matrix *the_matrix, Ploidy_Like_model *the_model,  Clade *the_genomes, WGX_Data *the_homologs, int start, int num_per_page, int focus, string gene_name, string filename, BOOL bitmap, BOOL IPC_call, BOOL rescale, std::stringstream *plot_ss, IMAGE_SIZE mysize, int *coords, int arrow_loc, BOOL arrow_left, int arrow_taxa, string tracked_name, string outname, string *prefix_map);

extern string make_gene_tree(Exchange *curr_exchange,  Tree *the_tree, Phylo_Matrix *the_matrix, TREE_TYPE my_type, WGX_Data *the_homologs, Ploidy_Like_model *the_model,  Clade *the_genomes, int pillar_num);

extern void construct_full_genome_lookup(string genomefile, Clade *the_genomes, std::map<std::string, int> gene_hash, std::map<std::string, int> &left_matches, std::map<std::string, int> &right_matches, std::map<std::string, int> &taxa_lookup, std::map<std::string, std::string> &left_name, std::map<std::string, std::string> &right_name, std::map<std::string, std::string> &aliases, BOOL have_loc_data);

extern void make_prefixes (Clade *the_genomes, std::string *prefixes, std::string *&prefix_map);

extern void get_max_prob_pattern(Exchange *curr_exchange, Phylo_Matrix *the_matrix, Ploidy_Like_model *the_model, int site, int &opt_track_id, int *taxa_track_ids, double &prob);

extern int draw_model_diag(Phylo_Matrix *the_matrix, string plotfile, BOOL bitmap, BOOL IPC_call, std::stringstream *plot_ss, int matrix_num);

extern int draw_rooted_tree(Exchange *curr_exchange, Tree *the_tree, string plotfile, double lnL, BOOL bitmap, BOOL IPC_call, std::stringstream *plot_ss);

MSG_TYPE get_message_type(char *buffer);

int communicate_plots(Exchange *curr_exchange, Tree *the_tree, Phylo_Matrix *the_matrix, Ploidy_Like_model *the_model,  Clade *the_genomes, WGX_Data *the_homologs, std::string socketname, std::string *&full_genome_files, std::string *&prefixes, BOOL have_loc_data)
{
    
    int i, frame_size, pillar, new_pillar, dupl_level, taxa, start, end, max, sock, msgsock, focus, my_track, gene_cutoff,
    rval, len, data_size, coords[7], tree_pillar, tree_size, string_pos, arrow_loc, arrow_taxa, min_genes,
    num_valid_pillars, *pillar_list, max_prob_pattern, *taxa_track_ids, num_genes, total_genes, **ortho_genes, *count_orthos, num_one, num_zero;
    char size_char, tree_type_char, rescale_char, ortho_char, subgenome_char, source_sub, rand_char;
    double percent_cutoff, max_prob, max_rate[3];
    struct sockaddr_un server;
    BOOL size_valid, is_bmp, name_error=FALSE, have_full=FALSE, arrow_left, full_match_no_syn=FALSE, alias_non_unique=FALSE, rescale=FALSE,
        orthologs_only, *use_pillars, all_one, all_zero, get_rand=FALSE;
    string gene_name, quit, size, filename, number_str, type_string, tree_string,tracked_name, outname, *prefix_map, alias_lc, alias;
    char buffer[BUFFER_SIZE], send_buffer[BUFFER_SIZE], dummy[1];
    MSG_TYPE my_msg;
    IMAGE_SIZE mysize=LARGE;
    Gene *my_gene;
    
    std::map<char, int> subgenome_lookup;
    std::map<std::string, int> gene_hash, full_left_matches, full_right_matches, full_taxa_ids, state_map;
    std::map<std::string, std::string> left_name, right_name, aliases, reverse_aliases;
    std::fstream fout;
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> distrib(0, the_homologs->get_num_homologs()-1);
    
    taxa_track_ids=new int[the_genomes->get_num_genomes()];
    use_pillars=new BOOL [the_homologs->get_num_homologs()];
    ortho_genes=new int* [the_homologs->get_dupl_level()];
    count_orthos = new int [the_homologs->get_dupl_level()];
    
    for (dupl_level=0; dupl_level<the_homologs->get_dupl_level(); dupl_level++) ortho_genes[dupl_level]=new int [the_genomes->get_num_genomes()];
        
    for (i=0; i<the_matrix->get_num_states(); i++)
        state_map[the_matrix->get_nth_state(i)->get_state_name()]=i;
    
    if (the_homologs->get_dupl_level()==2) {
        if (the_matrix->rate_matrix[state_map["Dupl"]][state_map["Copy1"]] > the_matrix->rate_matrix[state_map["Dupl"]][state_map["Copy2"]]) {
            subgenome_lookup['L']=0;
            subgenome_lookup['M']=1;
        }
        else {
            subgenome_lookup['L']=1;
            subgenome_lookup['M']=0;
        }
    }
    if (the_homologs->get_dupl_level()==3) {
        max_rate[0]=the_matrix->rate_matrix[state_map["Dupl"]][state_map["Copy1"]];
        max_rate[1]=the_matrix->rate_matrix[state_map["Dupl"]][state_map["Copy2"]];
        max_rate[2]=the_matrix->rate_matrix[state_map["Dupl"]][state_map["Copy3"]];
        
        if ((max_rate[0] > max_rate[1]) && (max_rate[0] > max_rate[2])) {
            subgenome_lookup['L']=0;
            if ((max_rate[1] >max_rate[2])) {
                subgenome_lookup['I']=1;
                subgenome_lookup['M']=2;
            }
            else {
                subgenome_lookup['I']=2;
                subgenome_lookup['M']=1;
            }
        }
        if ((max_rate[1] > max_rate[0]) && (max_rate[1] > max_rate[2])) {
            subgenome_lookup['L']=1;
            if ((max_rate[0] >max_rate[2])) {
                subgenome_lookup['I']=0;
                subgenome_lookup['M']=2;
            }
            else {
                subgenome_lookup['I']=2;
                subgenome_lookup['M']=0;
            }
        }
        if ((max_rate[2] > max_rate[0]) && (max_rate[2] > max_rate[1])) {
            subgenome_lookup['L']=2;
            if ((max_rate[1] >max_rate[0])) {
                subgenome_lookup['I']=1;
                subgenome_lookup['M']=0;
            }
            else {
                subgenome_lookup['I']=0;
                subgenome_lookup['M']=1;
            }
        }
    }
    
    //NOTE THAT THE BROWSER DOES NOT CURRENTLY SUPPORT WGQs!
    total_genes=0;
    for(pillar=0; pillar<the_homologs->get_num_homologs(); pillar++) {
        for(taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) {
            for(dupl_level=0; dupl_level< the_homologs->get_dupl_level(); dupl_level++) {
                if ((*the_model->the_tracks)[taxa].get_gene_track(pillar, dupl_level)->my_locus != 0) {
                    gene_hash[(*the_model->the_tracks)[taxa].get_gene_track(pillar, dupl_level)->my_locus->get_gene_obj((*the_model->the_tracks)[taxa].get_gene_track(pillar, dupl_level)->index_num)->get_name_string()]=pillar;
                    total_genes++;
                }
            }
        }
    }
    
    make_prefixes (the_genomes, prefixes, prefix_map);
    
    if (full_genome_files!=0) {
        for(taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) {
            if (full_genome_files[taxa] != "NONE") {
                have_full=TRUE;
                construct_full_genome_lookup(full_genome_files[taxa], the_genomes, gene_hash, full_left_matches, full_right_matches, full_taxa_ids, left_name, right_name, aliases, have_loc_data);
                
                
                for (std::map<string,string>::iterator it=aliases.begin(); it!=aliases.end(); ++it){
                    if (it->second != it->first) {
                        if(reverse_aliases.find(it->second) == reverse_aliases.end()) {
                            reverse_aliases[it->second]=it->first;
                        }
                    }
                }
            }
        }
    }
    
    the_model->get_gene_conditional_probs();
    
    ifstream f(socketname.c_str());
    if (f.good()) unlink(socketname.c_str());
    
    sock = socket(AF_UNIX, SOCK_STREAM, 0);
    if (sock < 0) {
        perror("opening stream socket");
        return(-1);
    }
    server.sun_family = AF_UNIX;
    strcpy(server.sun_path, socketname.c_str());
    if (bind(sock, (struct sockaddr *) &server, sizeof(struct sockaddr_un))) {
        perror("binding stream socket");
        exit(-1);
    }
    cout<<"Printing Socket has name"<<server.sun_path<<endl;
    listen(sock, 5);
    
    for (;;) {
        msgsock = accept(sock, 0, 0);
        if (msgsock == -1)
            perror("accept");
        else {
            do {
                bzero(buffer, sizeof(buffer));
                if ((rval = read(msgsock, buffer, BUFFER_SIZE)) < 0)
                    perror("reading stream message");
                else if (rval == 0)
                    printf("Ending connection\n");
                else {
                    
                    std::cout<<"Read buffer of "<<buffer<<std::endl<<std::flush;
                    my_msg=get_message_type(buffer);
                    std::stringstream ss (buffer);
                    
                    switch (my_msg){
                        case PING:
                            bzero(send_buffer, sizeof(send_buffer));
                            sprintf(send_buffer, "P\tPINGBACK\t%d\t%d\t%d\t%.2f", the_genomes->get_num_genomes(), the_homologs->get_num_homologs(), total_genes, POInT_version);
                            //strcpy(send_buffer, "P");
                            //strcat(send_buffer, "PINGBACK");
                            //strcat(send_buffer, "\t");
                            
                           // if (the_matrix->model_is_symmetric()==TRUE)
                            //    strcat(send_buffer, "S");
                           // else
                           //     strcat(send_buffer, "N");
                            
                            len=strlen(send_buffer);
                            
                            for (i=len; i<BUFFER_SIZE; i++) send_buffer[i]='\0';
                            
                            if (write(msgsock, send_buffer, sizeof(send_buffer)) < 0)
                                perror("writing on stream socket");
                            
                            bzero(send_buffer, sizeof(send_buffer));
                            strcpy(send_buffer, "E");
                            strcat(send_buffer, "ENDDATA");
                            
                            len=strlen(send_buffer);
                            
                            for (i=len; i<BUFFER_SIZE; i++) send_buffer[i]='\0';
                            
                            if (write(msgsock, send_buffer, sizeof(send_buffer)) < 0)
                                perror("writing on stream socket");
                            cout<<"Wrote ping back to query process"<<endl<<flush;
                            break;
                        case DATABLOCK:
                            ss>>dummy[0]>>gene_name>>size>>type_string>>new_pillar>>size_char>>rescale_char>>tree_pillar>>tree_type_char>>rand_char;
                            
                            get_rand=FALSE;
                            name_error=FALSE;
                            
                            switch (size_char) {
                                case 'L':
                                    mysize=LARGE;
                                    break;
                                case 'B':
                                    mysize=BIGGER;
                                    break;
                                case 'M':
                                    mysize=MEDIUM;
                                    break;
                                case 'S':
                                    mysize=SMALL;
                                    break;
                                case 'X':
                                    mysize=MASSIVE;
                                    break;
                                default:
                                    mysize=LARGE;
                                    break;
                            }
                            
                            if ((rescale_char=='r') || (rescale_char=='R')) rescale=TRUE;
                            
                            if ((rand_char == 'y') || (rand_char=='Y')) get_rand=TRUE;
                            
                            if (tree_pillar != -1) {
                                if ((tree_type_char=='p') || (tree_type_char=='P'))
                                    tree_string=make_gene_tree(curr_exchange,  the_tree, the_matrix, PHYLIP_TREE,
                                                               the_homologs, the_model,  the_genomes, tree_pillar);
                                else
                                    tree_string=make_gene_tree(curr_exchange,  the_tree, the_matrix, NEXUS_TREE,
                                                           the_homologs, the_model,  the_genomes, tree_pillar);
                                
                                cout<<"Tree is "<<tree_string<<endl;
                            }
                            
                            if (get_rand==TRUE)
                                gene_name="NONE";
                            
                            
                            cout<<"Rand char is "<<rand_char<<" Val is "<<get_rand<<" name is "<<gene_name<<endl;
                           
                            if (gene_name != "NONE") {
                                    name_error=TRUE;
                                alias_non_unique=FALSE;
                                tracked_name="NONE";
                                outname="NONE";
                                
                                if(gene_hash.find(gene_name) == gene_hash.end()) {
                                    alias_lc=gene_name;
                                    std::transform(alias_lc.begin(), alias_lc.end(),alias_lc.begin(), ::tolower);
                                    
                                    if (aliases.find(alias_lc) != aliases.end()) {
                                        if (aliases[alias_lc]=="NONE") alias_non_unique=TRUE;
                                        else
                                            gene_name=aliases[alias_lc];
                                    }
                                }
                                
                                if(gene_hash.find(gene_name) == gene_hash.end()) {
                                    if (have_full==TRUE) {
                                        full_match_no_syn=FALSE;
                                        if (!(full_left_matches.find(gene_name) == full_left_matches.end())) {
                                            if (full_left_matches[gene_name] != -1) {
                                                pillar=full_left_matches[gene_name];
                                                tracked_name=left_name[gene_name];
                                                outname=gene_name;
                                                focus=-1;
                                                arrow_left=TRUE;
                                                arrow_taxa=full_taxa_ids[gene_name];
                                                arrow_loc=pillar;
                                                name_error=FALSE;
                                            }
                                            else {
                                                if (full_right_matches[gene_name] != -1) {
                                                    pillar=full_right_matches[gene_name];
                                                    tracked_name=right_name[gene_name];
                                                    outname=gene_name;
                                                    focus=-1;
                                                    arrow_left=FALSE;
                                                    arrow_loc=pillar;
                                                    arrow_taxa=full_taxa_ids[gene_name];
                                                    name_error=FALSE;
                                                }
                                            }
                                            if ((full_right_matches[gene_name] ==-1) && full_left_matches[gene_name] == -1) full_match_no_syn=TRUE;
                                        }
                                    }
                                    
                                    
                                }
                                else {
                                    pillar = gene_hash[gene_name];
                                    focus=pillar;
                                    arrow_loc=-1;
                                    name_error=FALSE;
                                }
                                if (name_error==TRUE) {
                                    bzero(send_buffer, sizeof(send_buffer));
                                    if (alias_non_unique==TRUE) {
                                        strcpy(send_buffer, "M");
                                        strcat(send_buffer, "ERROR: Gene name ");
                                        strcat(send_buffer, gene_name.c_str());
                                        strcat(send_buffer, " is not a unquie identifier in this context.");
                                    }
                                    else {
                                        if (full_match_no_syn==TRUE) {
                                            strcpy(send_buffer, "M");
                                            strcat(send_buffer, "ERROR: Gene ");
                                            strcat(send_buffer, gene_name.c_str());
                                            strcat(send_buffer, " is not found within a valid synteny block.");
                                        }
                                        else {
                                            strcpy(send_buffer, "M");
                                            strcat(send_buffer, "ERROR: Gene ");
                                            strcat(send_buffer, gene_name.c_str());
                                            strcat(send_buffer, " was not found in these data.");
                                        }
                                    }
                                    len=strlen(send_buffer);
                                    
                                    for (i=len; i<BUFFER_SIZE; i++) send_buffer[i]='\0';
                                    
                                }
                                else {
                                
                                    cout<<"Gene is "<<gene_name<<" size is "<<size<<" type is "<<type_string<<endl;
                                    
                                    //pillar = gene_hash[gene_name];
                                    
                                    frame_size=std::stoi(size,nullptr);
                                    frame_size=abs(frame_size);
                                    
                                    max=MAX_FRAME;
                                    if (frame_size > max) {
                                        frame_size=MAX_FRAME;
                                        cout<<"Note: maximum frame size exceeded. Reseting to "<<frame_size<<"\n";
                                    }
                                    
                                    std::transform(type_string.begin(), type_string.end(), type_string.begin(),::tolower);
                                    
                                    if (type_string.front() == 'y') is_bmp=TRUE;
                                    else is_bmp=FALSE;
                                    
                                    start = pillar - (frame_size /2);
                                    
                                    if (start < 0) start=0;
                                    
                                    end=start+frame_size;
                                    
                                    if (end >= the_homologs->get_num_homologs()) {
                                        end =the_homologs->get_num_homologs()-1;
                                        start = end - frame_size+1;
                                    }
                                    
                                    stringstream ss2;
                                    ss2 << frame_size;
                                    number_str = ss2.str();
                                    filename = gene_name;
                                    filename += "tracking_";
                                    filename += number_str;
                                    if (is_bmp == TRUE)
                                    filename += ".png";
                                    else
                                    filename += ".eps";
                                }
                            }
                            else {
                                frame_size=std::stoi(size,nullptr);
                                frame_size=abs(frame_size);
                                
                                max=MAX_FRAME;
                                if (frame_size > max) {
                                    frame_size=MAX_FRAME;
                                    cout<<"Note: maximum frame size exceeded. Reseting to "<<frame_size<<"\n";
                                }
                                
                                if (get_rand==TRUE) {
                                    pillar = distrib(gen);
                                    if (pillar < (frame_size/2)) pillar = frame_size/2;
                                    
                                    if (pillar > ((the_homologs->get_num_homologs()-1)-(frame_size/2))) pillar = (the_homologs->get_num_homologs()-1)-(frame_size/2);
                                    
                                    focus=pillar;
                                    cout<<"Chosen random pillar is "<<pillar<<endl;
                                }
                                else {
                                    if ((new_pillar < 0) || (new_pillar>=the_homologs->get_num_homologs())) {
                                        name_error=TRUE;
                                        
                                        strcpy(send_buffer, "M");
                                        strcat(send_buffer, "ERROR: Requested pillar ");
                                        stringstream ss2;
                                        ss2 << new_pillar;
                                        number_str = ss2.str();
                                        strcat(send_buffer, number_str.c_str());
                                        strcat(send_buffer, " out of bounds.");
                                        len=strlen(send_buffer);
                                        
                                        for (i=len; i<BUFFER_SIZE; i++) send_buffer[i]='\0';
                                        
                                    }
                                    else {
                                        pillar=new_pillar;
                                        focus=pillar;
                                        
                                    }
                                }
                                
                                if (name_error==FALSE) {
                                        std::transform(type_string.begin(), type_string.end(), type_string.begin(),::tolower);
                                        
                                        if (type_string.front() == 'y') is_bmp=TRUE;
                                        else is_bmp=FALSE;
                                        
                                        start = pillar - (frame_size /2);
                                        
                                        if (start < 0) start=0;
                                        
                                        end=start+frame_size;
                                        
                                        if (end >= the_homologs->get_num_homologs()) {
                                            end =the_homologs->get_num_homologs()-1;
                                            start = end - frame_size+1;
                                        }
                                        
                                        stringstream ss2;
                                        ss2 << pillar << frame_size;
                                        number_str = ss2.str();
                                        filename = "Pillar";
                                        filename += number_str;
                                        filename += "_tracking";
                                        
                                        if (is_bmp == TRUE)
                                            filename += ".png";
                                        else
                                            filename += ".eps";
                                        
                                }
                            }
                            
                            
                            if (name_error == TRUE) {
                                if (write(msgsock, send_buffer, sizeof(send_buffer)) < 0)
                                    perror("writing on stream socket");
                                
                                bzero(send_buffer, sizeof(send_buffer));
                                strcpy(send_buffer, "E");
                                strcat(send_buffer, "ENDDATA");
                                
                                len=strlen(send_buffer);
                                
                                for (i=len; i<BUFFER_SIZE; i++) send_buffer[i]='\0';
                                
                                if (write(msgsock, send_buffer, sizeof(send_buffer)) < 0)
                                    perror("writing on stream socket");
                            }
                            else {
                                
                                
                                std::stringstream plot_ss;
                                //auto old_buf = std::cout.rdbuf(plot_ss.rdbuf());
                                draw_region (curr_exchange, the_tree, the_matrix, the_model, the_genomes, the_homologs, start, frame_size, focus, gene_name, filename, is_bmp,TRUE, rescale, &plot_ss, mysize, coords, arrow_loc, arrow_left, arrow_taxa, tracked_name, outname, prefix_map);

                                std::stringstream genedata_ss;
                                
                                
                                data_size=0;
                                while(!(plot_ss.eof())) {
                                    dummy[0]=(char)plot_ss.get();
                                //while (dummy[0]=(char)plot_ss.get()) {
                                //while (plot_ss>>dummy[0]) {
                                    data_size++;
                                }
                                data_size--;
                                cout<<"Computed total size of "<<data_size<<endl<<flush;
                                bzero(send_buffer, sizeof(send_buffer));
                                if (tree_pillar!=-1) {
                                    tree_size=tree_string.length();
                                }
                                else tree_size=0;
                                
                                //send_buffer[0]='D';
                                //strcat(send_buffer, "\t");
                                sprintf(send_buffer, "S\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d", data_size, pillar, start, end, the_homologs->get_num_homologs(), the_homologs->get_dupl_level(), the_genomes->get_num_genomes(),
                                        coords[0],coords[1],coords[2],coords[3],coords[4],coords[5],coords[6], tree_size);
                                
                                cout<<"Sending coords: "<<coords[0]<<" "<<coords[1]<<" "<<coords[2]<<" "<<coords[3]<<" "<<coords[4]<<" "<<coords[5]<<" "<<coords[6]<<endl;
                                
                                //for(i=0; i<7; i++)
                                //    sprintf(send_buffer, "\t%d", coords[i]);
                                
                                cout<<"Sending size data: "<<send_buffer<<endl<<flush;
                                
                                if (write(msgsock, send_buffer, sizeof(send_buffer)) < 0)
                                    perror("writing on stream socket");
                                
                                plot_ss.clear();
                                plot_ss.seekg(0, ios::beg);
                                
                                while(!(plot_ss.eof())) {
                                    bzero(send_buffer, sizeof(send_buffer));
                                    send_buffer[0]='D';
                                    i=1;
                                    
                                    while((i<BUFFER_SIZE) &&(!(plot_ss.eof()))) {
                                        //send_buffer[i]=(char)plot_ss.get();
                                        dummy[0]=(char)plot_ss.get();
                                        if (!(plot_ss.eof()))
                                            send_buffer[i]=dummy[0];
                                        //std::cout<<"i ="<<i<<" = "<<send_buffer[i]<<std::endl<<std::flush;
                                        i++;
                                    }
                                    
                                    //std::cout<<"BUF is "<<send_buffer<<endl;
                                    
                                    if(i<BUFFER_SIZE) std::cout<<"Final msg i = "<<i<<std::endl<<std::flush;
                                    
                                    while(i<BUFFER_SIZE) {
                                        send_buffer[i]='\0';
                                        i++;
                                    }
                                    
                                    std::cout<<"Sending image message of size "<<sizeof(send_buffer)<<std::endl<<std::flush;
                                    
                                    if (write(msgsock, send_buffer, sizeof(send_buffer)) < 0)
                                        perror("writing on stream socket");
                                }
                                
                                if (tree_pillar !=-1) {
                                    string_pos=0;
                                    
                                    while(string_pos<tree_string.length()) {
                                        bzero(send_buffer, sizeof(send_buffer));
                                        strcpy(send_buffer, "T");
                                        i=1;
                                        
                                        while((i<BUFFER_SIZE) && (string_pos < tree_string.length())) {
                                            send_buffer[i]=tree_string[string_pos];
                                            i++;
                                            string_pos++;
                                        }
                                        
                                        while (i<BUFFER_SIZE) {
                                            send_buffer[i]='\0';
                                            i++;
                                        }
                                        cout<<"Tree msg "<<send_buffer<<endl;
                                        if (write(msgsock, send_buffer, sizeof(send_buffer)) < 0)
                                            perror("writing on stream socket");
                                    }
                                    
                                }
                                
                                if (have_loc_data ==TRUE) {
                                    std::stringstream loc_ss;
                                    for(taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) {
                                        loc_ss<<(*the_genomes)[taxa].get_name_string()<<"\t"<<(*the_genomes)[taxa].get_web_link()<<"\t";
                                        cout<<"Sending "<<(*the_genomes)[taxa].get_web_link()<<" for "<<(*the_genomes)[taxa].get_name_string()<<endl;
                                    }
                                    
                                    for(i=0; i<frame_size; i++) {
                                        get_max_prob_pattern(curr_exchange, the_matrix, the_model, i+start, max_prob_pattern, taxa_track_ids, max_prob);
                                        for (dupl_level=0; dupl_level<the_homologs->get_dupl_level(); dupl_level++) {
                                            for(taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) {
                                                cout<<"Position is "<<i+start<<" d: "<<dupl_level<<" Taxa "<<taxa
                                                <<" TP: "<<taxa_track_ids[taxa]<<endl<<flush;
                                                if ((*the_model->the_tracks)[taxa].get_gene_track(i+start, the_model->tracking_permutes[taxa_track_ids[taxa]][dupl_level])->my_locus != 0) {
                                                    my_gene=(*the_model->the_tracks)[taxa].get_gene_track(i+start, the_model->tracking_permutes[taxa_track_ids[taxa]][dupl_level])->my_locus->get_gene_obj((*the_model->the_tracks)[taxa].get_gene_track(i+start, the_model->tracking_permutes[taxa_track_ids[taxa]][dupl_level])->index_num);
                                                    cout<<"Looking at "<<my_gene->get_name_string()<<endl;
                                                    loc_ss<<my_gene->get_name_string()<<"\t";
                                                    
                                                    //alias="NONE";
                                                    //for (std::map<string,string>::iterator it=aliases.begin(); it!=aliases.end(); ++it){
                                                     //   if (it->second == my_gene->get_name_string()) {
                                                      //      if (( it->first != it->second ) &&(alias == "NONE")) alias=it->first;
                                                      //  }
                                                   // }
                                                    if(reverse_aliases.find(my_gene->get_name_string()) != reverse_aliases.end()) {
                                                        alias=reverse_aliases[my_gene->get_name_string()];
                                                    }
                                                    else alias="NONE";
                                                    
                                                    cout<<"Found alias of "<<alias<<endl;
                                                    loc_ss<<alias<<"\t";
                                                    loc_ss<<my_gene->get_chrom_name()<<"\t";
                                                    loc_ss<<my_gene->get_start_pos()<<"\t";
                                                    loc_ss<<my_gene->get_end_pos()<<"\t";
                                                }
                                                else {
                                                    loc_ss<<"NONE\tNONE\tNONE\tNONE\tNONE\t";
                                                }
                                            }
                                        }
                                    }
                                    
                                    data_size=0;
                                    while(!(loc_ss.eof())) {
                                        dummy[0]=(char)loc_ss.get();
                                    //while (dummy[0]=(char)plot_ss.get()) {
                                    //while (plot_ss>>dummy[0]) {
                                        data_size++;
                                    }
                                    data_size--;
                                    
                                    loc_ss.clear();
                                    loc_ss.seekg(0, ios::beg);
                                    
                                    while(!(loc_ss.eof())) {
                                        bzero(send_buffer, sizeof(send_buffer));
                                        send_buffer[0]='G';
                                        i=1;
                                        
                                        while((i<BUFFER_SIZE) &&(!(loc_ss.eof()))) {
                                            //send_buffer[i]=(char)plot_ss.get();
                                            dummy[0]=(char)loc_ss.get();
                                            if (!(loc_ss.eof()))
                                                send_buffer[i]=dummy[0];
                                            //std::cout<<"i ="<<i<<" = "<<send_buffer[i]<<std::endl<<std::flush;
                                            i++;
                                        }
                                        
                                        std::cout<<"GENE BUF is "<<send_buffer<<endl;
                                        
                                        if(i<BUFFER_SIZE) std::cout<<"Final msg i = "<<i<<std::endl<<std::flush;
                                        
                                        while(i<BUFFER_SIZE) {
                                            send_buffer[i]='\0';
                                            i++;
                                        }
                                        
                                        std::cout<<"Sending message of size "<<sizeof(send_buffer)<<std::endl<<std::flush;
                                        
                                        if (write(msgsock, send_buffer, sizeof(send_buffer)) < 0)
                                            perror("writing on stream socket");
                                    }
                                    
                                }
                                
                                bzero(send_buffer, sizeof(send_buffer));
                                strcpy(send_buffer, "E");
                                strcat(send_buffer, "ENDDATA");
                                len=strlen(send_buffer);
                                
                                for (i=len; i<BUFFER_SIZE; i++) send_buffer[i]='\0';
                                std::cout<<"Sending end of data stream"<<std::endl<<std::flush;
                                if (write(msgsock, send_buffer, sizeof(send_buffer)) < 0)
                                    perror("writing on stream socket");
                            }
                            break;
                        case BATCHBLOCK:
                            ss>>dummy[0]>>percent_cutoff>>min_genes>>ortho_char>>tree_type_char>>subgenome_char;
                            
                            cout<<"PC: "<<percent_cutoff<<"| MG: |"<<min_genes<<" OC: "<<ortho_char<<" TC: "<<tree_type_char<<" SGC: "<<subgenome_char<<endl;
                            
                            if ((ortho_char=='y') || (ortho_char=='Y')) {
                                orthologs_only=TRUE;
                                cout<<"Downloading only single-copy orthologs. Ortho char is "<<ortho_char<<" Subgenome L:"<<subgenome_lookup['L']<<" Subgenome M:"<<subgenome_lookup['M']<<"\n";
                            }
                            else {
                                orthologs_only=FALSE;
                                gene_cutoff = (the_homologs->get_dupl_level()*the_genomes->get_num_genomes())-min_genes;
                            }
                            num_valid_pillars=0;
                            
                            if (orthologs_only==TRUE) {
                                for(pillar=0; pillar<the_homologs->get_num_homologs(); pillar++) {
                                    use_pillars[pillar]=FALSE;
                                    get_max_prob_pattern(curr_exchange, the_matrix, the_model, pillar, max_prob_pattern, taxa_track_ids, max_prob);
                                  
                                    if (max_prob>=percent_cutoff) {
                                        for(taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) {
                                            for (dupl_level=0; dupl_level<the_homologs->get_dupl_level(); dupl_level++)
                                                ortho_genes[dupl_level][taxa]=0;
                                        }
                                        for (dupl_level=0; dupl_level<the_homologs->get_dupl_level(); dupl_level++) count_orthos[dupl_level]=0;
                                        
                                        for(taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) {
                                            for (dupl_level=0; dupl_level<the_homologs->get_dupl_level(); dupl_level++) {
                                                my_track=0;
                                                while(the_model->tracking_permutes[taxa_track_ids[taxa]][my_track] != dupl_level) my_track++;
                                                if ((*the_model->the_tracks)[taxa].get_gene_track(pillar, dupl_level)->my_locus != 0) ortho_genes[my_track][taxa]=1;
                                            }
                                        }
                                        
                                        for (dupl_level=0; dupl_level<the_homologs->get_dupl_level(); dupl_level++) {
                                            for(taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) count_orthos[dupl_level]+= ortho_genes[dupl_level][taxa];
                                        }
                                        
                                        num_one=0;
                                        num_zero=0;
                                        
                                        for (dupl_level=0; dupl_level<the_homologs->get_dupl_level(); dupl_level++) {
                                            if (count_orthos[dupl_level]==the_genomes->get_num_genomes()) num_one++;
                                            if (count_orthos[dupl_level]==0) num_zero++;
                                        }
                                        
                                        cout<<"Pillar "<<pillar<<" N1: "<<num_one<<" NZ: "<<num_zero<<endl;
                                        
                                        if ((num_one ==1) && (num_zero == (the_homologs->get_dupl_level()-1))) {
                                            if (subgenome_char=='A') {
                                                use_pillars[pillar]=TRUE;
                                                num_valid_pillars++;
                                            }
                                            else {
                                                for (dupl_level=0; dupl_level<the_homologs->get_dupl_level(); dupl_level++)
                                                    if (count_orthos[dupl_level]==the_genomes->get_num_genomes()) source_sub=dupl_level;
                                                
                                                if (subgenome_lookup[subgenome_char] == source_sub) {
                                                    use_pillars[pillar]=TRUE;
                                                    num_valid_pillars++;
                                                }
                                            }
                                        }
                                        
                                        if (use_pillars[pillar]==TRUE) cout<<"Pillar "<<pillar<<" is used\n";
                                        
                                    }
                                }
                                    
                            }
                            else {
                                for(pillar=0; pillar<the_homologs->get_num_homologs(); pillar++) {
                                    use_pillars[pillar]=FALSE;
                                    get_max_prob_pattern(curr_exchange, the_matrix, the_model, pillar, max_prob_pattern, taxa_track_ids, max_prob);
                                    
                                    num_genes=0;
                                    if (max_prob>=percent_cutoff) {
                                        for(taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) {
                                            for (dupl_level=0; dupl_level<the_homologs->get_dupl_level(); dupl_level++) {
                                                if ((*the_model->the_tracks)[taxa].get_gene_track(pillar, dupl_level)->my_locus != 0) num_genes++;
                                            }
                                        }
                                        
                                        if (num_genes>=gene_cutoff) {
                                            use_pillars[pillar]=TRUE;
                                            num_valid_pillars++;
                                        }
                                    }
                                }
                            }
                                
                            if (num_valid_pillars==0) {
                                bzero(send_buffer, sizeof(send_buffer));
                                strcpy(send_buffer, "M");
                                strcat(send_buffer, "ERROR: No pillars matching these criteria were found.");
                                len=strlen(send_buffer);
                                cout<<"Sending "<<send_buffer<<endl;
                                for (i=len; i<BUFFER_SIZE; i++) send_buffer[i]='\0';
                                if (write(msgsock, send_buffer, sizeof(send_buffer)) < 0)
                                    perror("writing on stream socket");
                                
                                bzero(send_buffer, sizeof(send_buffer));
                                strcpy(send_buffer, "E");
                                strcat(send_buffer, "ENDDATA");
                                
                                len=strlen(send_buffer);
                                
                                for (i=len; i<BUFFER_SIZE; i++) send_buffer[i]='\0';
                                
                                if (write(msgsock, send_buffer, sizeof(send_buffer)) < 0)
                                    perror("writing on stream socket");
                            }
                            else {
                                cout<<"Found "<<num_valid_pillars<<" to send\n";
                                std::stringstream tree_ss;
                                //tree_ss<<num_valid_pillars;
                                
                                for(pillar=0; pillar<the_homologs->get_num_homologs(); pillar++) {
                                    if (use_pillars[pillar]==TRUE) {
                                        tree_string.clear();
                                        if ((tree_type_char=='p') || (tree_type_char=='P'))
                                            tree_string=make_gene_tree(curr_exchange,  the_tree, the_matrix, PHYLIP_TREE,
                                                                       the_homologs, the_model,  the_genomes, pillar);
                                        else
                                            tree_string=make_gene_tree(curr_exchange,  the_tree, the_matrix, NEXUS_TREE,
                                                                       the_homologs, the_model,  the_genomes, pillar);
                                        
                                        tree_ss<<pillar<<'\0'<<tree_string<<'\0';
                                    }
                                }
                                
                                data_size=0;
                                while(!(tree_ss.eof())) {
                                    dummy[0]=(char)tree_ss.get();
                                    data_size++;
                                }
                                data_size--;
                                cout<<"Computed total size of "<<data_size<<endl<<flush;
                                bzero(send_buffer, sizeof(send_buffer));
                                
                                sprintf(send_buffer, "S\t%d\t%d", data_size, num_valid_pillars);
                                
                                //for(i=0; i<7; i++)
                                //    sprintf(send_buffer, "\t%d", coords[i]);
                                
                                cout<<"Sending size data: "<<send_buffer<<endl<<flush;
                                
                                if (write(msgsock, send_buffer, sizeof(send_buffer)) < 0)
                                    perror("writing on stream socket");
                                
                                tree_ss.clear();
                                tree_ss.seekg(0, ios::beg);
                                
                                while(!(tree_ss.eof())) {
                                    bzero(send_buffer, sizeof(send_buffer));
                                    send_buffer[0]='B';
                                    i=1;
                                    
                                    while((i<BUFFER_SIZE) &&(!(tree_ss.eof()))) {
                                        //send_buffer[i]=(char)plot_ss.get();
                                        dummy[0]=(char)tree_ss.get();
                                        if (!(tree_ss.eof()))
                                            send_buffer[i]=dummy[0];
                                        //std::cout<<"i ="<<i<<" = "<<send_buffer[i]<<std::endl<<std::flush;
                                        i++;
                                    }
                                    
                                    //std::cout<<"BUF is "<<send_buffer<<endl;
                                    
                                    if(i<BUFFER_SIZE) std::cout<<"Final msg i = "<<i<<std::endl<<std::flush;
                                    
                                    while(i<BUFFER_SIZE) {
                                        send_buffer[i]='\0';
                                        i++;
                                    }
                                    
                                    std::cout<<"Sending message of size "<<sizeof(send_buffer)<<std::endl<<std::flush;
                                    
                                    if (write(msgsock, send_buffer, sizeof(send_buffer)) < 0)
                                        perror("writing on stream socket");
                                }
                                    
                                bzero(send_buffer, sizeof(send_buffer));
                                strcpy(send_buffer, "E");
                                strcat(send_buffer, "ENDDATA");
                                len=strlen(send_buffer);
                                
                                for (i=len; i<BUFFER_SIZE; i++) send_buffer[i]='\0';
                                std::cout<<"Sending end of data stream"<<std::endl<<std::flush;
                                if (write(msgsock, send_buffer, sizeof(send_buffer)) < 0)
                                    perror("writing on stream socket");
                                    
                                break;
                                
                            }
                            case MODELDIAG:
                                {
                                    filename="ModelDiag_Dummy.png";
                                    std::stringstream model_ss;
                                    draw_model_diag(the_matrix, filename, TRUE, TRUE, &model_ss, 0);
                                
                                    data_size=0;
                                    while(!(model_ss.eof())) {
                                        dummy[0]=(char)model_ss.get();
                                        data_size++;
                                    }
                                    data_size--;
                                    cout<<"Computed total size of "<<data_size<<endl<<flush;
                                    bzero(send_buffer, sizeof(send_buffer));
                                
                                    sprintf(send_buffer, "S\t%d", data_size);
                                
                                
                                    cout<<"Sending size data: "<<send_buffer<<endl<<flush;
                                
                                    if (write(msgsock, send_buffer, sizeof(send_buffer)) < 0)
                                        perror("writing on stream socket");
                                
                                    model_ss.clear();
                                    model_ss.seekg(0, ios::beg);
                                
                                    while(!(model_ss.eof())) {
                                        bzero(send_buffer, sizeof(send_buffer));
                                        send_buffer[0]='I';
                                        i=1;
                                    
                                        while((i<BUFFER_SIZE) &&(!(model_ss.eof()))) {
                                            //send_buffer[i]=(char)plot_ss.get();
                                            dummy[0]=(char)model_ss.get();
                                            if (!(model_ss.eof()))
                                                send_buffer[i]=dummy[0];
                                            //std::cout<<"i ="<<i<<" = "<<send_buffer[i]<<std::endl<<std::flush;
                                            i++;
                                        }
                                    
                                        //std::cout<<"BUF is "<<send_buffer<<endl;
                                    
                                        if(i<BUFFER_SIZE) std::cout<<"Final msg i = "<<i<<std::endl<<std::flush;
                                    
                                        while(i<BUFFER_SIZE) {
                                            send_buffer[i]='\0';
                                            i++;
                                        }
                                    
                                        std::cout<<"Sending message of size "<<sizeof(send_buffer)<<std::endl<<std::flush;
                                    
                                        if (write(msgsock, send_buffer, sizeof(send_buffer)) < 0)
                                            perror("writing on stream socket");
                                    }
                                
                                    bzero(send_buffer, sizeof(send_buffer));
                                    strcpy(send_buffer, "E");
                                    strcat(send_buffer, "ENDDATA");
                                    len=strlen(send_buffer);
                                
                                    for (i=len; i<BUFFER_SIZE; i++) send_buffer[i]='\0';
                                    std::cout<<"Sending end of data stream"<<std::endl<<std::flush;
                                    if (write(msgsock, send_buffer, sizeof(send_buffer)) < 0)
                                        perror("writing on stream socket");
                                }
                            break;
                            
                            case TREEDIAG:
                                {
                                    filename="TreeDiag_Dummy.png";
                                    std::stringstream treediag_ss;
                                    draw_rooted_tree(curr_exchange, the_tree, filename, curr_exchange->get_saved_lnL(), TRUE, TRUE, &treediag_ss);
                            
                                    data_size=0;
                                    while(!(treediag_ss.eof())) {
                                        dummy[0]=(char)treediag_ss.get();
                                        data_size++;
                                    }
                                    data_size--;
                                    cout<<"Computed total size of "<<data_size<<endl<<flush;
                                    bzero(send_buffer, sizeof(send_buffer));
                            
                                    sprintf(send_buffer, "S\t%d", data_size);
                            
                            
                                    cout<<"Sending size data: "<<send_buffer<<endl<<flush;
                            
                                    if (write(msgsock, send_buffer, sizeof(send_buffer)) < 0)
                                        perror("writing on stream socket");
                            
                                    treediag_ss.clear();
                                    treediag_ss.seekg(0, ios::beg);
                            
                                    while(!(treediag_ss.eof())) {
                                        bzero(send_buffer, sizeof(send_buffer));
                                        send_buffer[0]='V';
                                        i=1;
                                
                                        while((i<BUFFER_SIZE) &&(!(treediag_ss.eof()))) {
                                        //send_buffer[i]=(char)plot_ss.get();
                                            dummy[0]=(char)treediag_ss.get();
                                            if (!(treediag_ss.eof()))
                                                send_buffer[i]=dummy[0];
                                            //std::cout<<"i ="<<i<<" = "<<send_buffer[i]<<std::endl<<std::flush;
                                            i++;
                                        }
                                
                                    //std::cout<<"BUF is "<<send_buffer<<endl;
                                
                                        if(i<BUFFER_SIZE) std::cout<<"Final msg i = "<<i<<std::endl<<std::flush;
                                
                                        while(i<BUFFER_SIZE) {
                                            send_buffer[i]='\0';
                                            i++;
                                        }
                                
                                        std::cout<<"Sending message of size "<<sizeof(send_buffer)<<std::endl<<std::flush;
                                
                                        if (write(msgsock, send_buffer, sizeof(send_buffer)) < 0)
                                            perror("writing on stream socket");
                                    }
                            
                                    bzero(send_buffer, sizeof(send_buffer));
                                    strcpy(send_buffer, "E");
                                    strcat(send_buffer, "ENDDATA");
                                    len=strlen(send_buffer);
                            
                                    for (i=len; i<BUFFER_SIZE; i++) send_buffer[i]='\0';
                                    std::cout<<"Sending end of data stream"<<std::endl<<std::flush;
                                    if (write(msgsock, send_buffer, sizeof(send_buffer)) < 0)
                                        perror("writing on stream socket");
                                }
                            
                        break;
                            
                        case GENOMELIST:
                            {
                                stringstream name_ss;
                                bzero(send_buffer, sizeof(send_buffer));
                                name_ss<<"G\t"<<the_genomes->get_num_genomes();
                                
                                //sprintf(send_buffer, "G\t%d", the_genomes->get_num_genomes());
                                for(taxa=0; taxa<the_genomes->get_num_genomes(); taxa++)
                                    name_ss<<"\t"<<(*the_genomes)[taxa].get_name_string();
                                    //sprintf(send_buffer, "\t%s", (*the_genomes)[i].get_name_string().c_str());
                                name_ss<<"\0";
                                i=0;
                                name_ss.clear();
                                name_ss.seekg(0, ios::beg);
                                while(!(name_ss.eof())) {
                                //send_buffer[i]=(char)plot_ss.get();
                                    dummy[0]=(char)name_ss.get();
                                    if (!(name_ss.eof()))
                                    send_buffer[i]=dummy[0];
                                    i++;
                                }
                                //send_buffer[i]='\0';
                                cout<<"Sending buffer "<<send_buffer<<endl;
                              
                                if (write(msgsock, send_buffer, sizeof(send_buffer)) < 0)
                                    perror("writing on stream socket");
                        
                               
                        
                                bzero(send_buffer, sizeof(send_buffer));
                                strcpy(send_buffer, "E");
                                strcat(send_buffer, "ENDDATA");
                                len=strlen(send_buffer);
                        
                                for (i=len; i<BUFFER_SIZE; i++) send_buffer[i]='\0';
                                std::cout<<"Sending end of data stream"<<std::endl<<std::flush;
                                if (write(msgsock, send_buffer, sizeof(send_buffer)) < 0)
                                    perror("writing on stream socket");
                            }
                        
                    break;
                          
                            
                    }
                    

                  
                }
            } while (rval > 0);
        }
            close(msgsock);
    }
    close(sock);
    unlink(socketname.c_str());
    delete[] taxa_track_ids;
    delete[] use_pillars;
    delete[] count_orthos;
    
    for (dupl_level=0; dupl_level<the_homologs->get_dupl_level(); dupl_level++) delete[] ortho_genes[dupl_level];
    delete[] ortho_genes;
    
}

MSG_TYPE get_message_type(char *buffer)
{
    MSG_TYPE msg_type=ENDCOMM;
    
    if(buffer[0] == 'M') msg_type=MESSAGE;
    if(buffer[0] == 'D') msg_type=DATABLOCK;
    if(buffer[0] == 'T') msg_type=PHYLOTREE;
    if(buffer[0] == 'P') msg_type=PING;
    if(buffer[0] == 'B') msg_type=BATCHBLOCK;
    if(buffer[0] == 'I') msg_type=MODELDIAG;
    if(buffer[0] == 'V') msg_type=TREEDIAG;
    if(buffer[0] == 'G') msg_type=GENOMELIST;
    return(msg_type);
}
