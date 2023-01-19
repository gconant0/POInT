#include "anneal_template.cpp"
#include "track_anneal_recombo.h"
#include <queue>

#ifdef _OPEN_MP_VERSION_
#include "omp.h"
#endif


//#include <sys/time.h>
//#include <sys/resource.h>

//#define __SHOW_STEPS__

//unsigned long RunTime();  
using namespace::std;

Track_point::Track_point(Exchange *cexchange)
{

  int i, taxa, dupl_level, *uniform_order, choose_point, choose_place;
  
  curr_exchange=new Exchange();
  
    
  //Copies the exchange object into the new exchange
  //(Prevents conflicts between the exchanges of different walkers)
  (*curr_exchange)=(*cexchange);
    max_breaks=curr_exchange->get_max_breaks();
    
    cout<<"Allowing "<<max_breaks<<" without inducing block break\n";
    
  the_genomes=global_genomes;
  the_homologs=global_homologs;

    //the_tracks=0;

    the_tracks = new WGX_Tracks(global_homologs, global_genomes);
    homolog_order=new int [the_homologs->get_num_homologs()];
 
    for(i=0; i<the_homologs->get_num_homologs(); i++)
        homolog_order[i]=i;
    
    if (cexchange->use_guessed_order() == TRUE)
        guess_order();
    
    in_valid_break=new BOOL* [the_genomes->get_num_genomes()];
    for (taxa=0; taxa<the_genomes->get_num_genomes(); taxa++)
        in_valid_break[taxa]=new BOOL[the_homologs->get_dupl_level()];
    
    location_density=new int [the_homologs->get_num_homologs()*the_genomes->get_num_genomes()*the_homologs->get_dupl_level()];
    
    has_ends=new BOOL*[the_genomes->get_num_genomes()];
    for (taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) {
        has_ends[taxa]=new BOOL [the_homologs->get_dupl_level()];
    }
    
    block_set = new Gene_block[the_homologs->get_num_homologs()];
    //set_blocks();
    set_blocks_v2();
    //Make the density array as big as it can ever get
    
    
    
    score_valid=TRUE;
  
    cout<<"Initial block count: "<<num_blocks<<endl;
    cout<<"Initial score: "<<get_score()<<endl;
}




Track_point& Track_point::operator=(Track_point & assign_from)
{
    int i;
    
    for(i=0; i<the_homologs->get_num_homologs(); i++)
       homolog_order[i]=assign_from.homolog_order[i];
    
    //set_blocks();
    score_valid=FALSE;
    return(*this);
}



int Track_point::get_block_start(int bl)
{
    return(block_set[bl].start);
}

int Track_point::get_block_end(int bl)
{
    return(block_set[bl].end);
}


void Track_point::set_blocks()
{
    int curr_pillar, next_pillar, taxa, curr_block_start, min_end, dupl_index, num_links, break_count, num_complete=0;
    char link;
    Gene_Track_List_DX *my_pillar_site, *next_pillar_site;
    BOOL has_all_neighbors,  end_before;
    
    for (taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) {
        for(dupl_index=0; dupl_index<the_homologs->get_dupl_level(); dupl_index++) in_valid_break[taxa][dupl_index]=FALSE;
    }
    check_order();
    
    the_tracks->change_order(homolog_order);
    the_tracks->update_tracking();
   

    curr_pillar=0;
    curr_block_start=0;
    
    num_blocks=0;
    
    while(curr_pillar < (the_homologs->get_num_homologs()-1)) {
        next_pillar=curr_pillar+1;
        
        //has_all_neighbors=TRUE;
        //taxa=0;
        //end_before=FALSE;
        
        break_count=0;
        for(taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) {
            for(dupl_index=0; dupl_index<the_homologs->get_dupl_level(); dupl_index++) {
                if ((*the_tracks)[taxa].has_back_link(next_pillar, dupl_index) ==FALSE) {break_count++; link ='B';}
                else link = 'L';
                cout<<taxa<<": "<<dupl_index<<" = "
                <<(*the_tracks)[taxa].get_gene_track(next_pillar, dupl_index)->my_locus<<" has "<<link<<endl;
            }
        }
        //while((taxa<the_genomes->get_num_genomes()) && (has_all_neighbors == TRUE)) {
        //    for(dupl_index=0; dupl_index<the_homologs->get_dupl_level(); dupl_index++) {

          //   v has_all_neighbors=FALSE;

          //  }
           // taxa++;
        //}
        
        if (break_count ==0 ) num_complete++;
        if (break_count > max_breaks) {
        //if (has_all_neighbors==FALSE) {
            
                block_set[num_blocks].start=curr_block_start;
                block_set[num_blocks].end=curr_pillar;
                num_blocks++;
                curr_block_start=next_pillar;
        }
        cout<<"Pillar "<<curr_pillar<<" at block "<<num_blocks<<" has "<<break_count<<" breaks "<<" Complete: "<<num_complete<<endl;
        curr_pillar++;
    }
    
    //Account for the final block
    block_set[num_blocks].start=curr_block_start;
    block_set[num_blocks].end=curr_pillar+1;
    num_blocks++;

    set_density();
    score_valid=TRUE;
    
    //cout<<" Setting "<<num_blocks<<" blocks and num_complete locations : "<<num_complete<<"\n";
   // cout<<"Setting "<<num_blocks<<" blocks. End is "<<block_set[num_blocks-1].end<<" Full Breaks: "<<the_tracks->get_num_full_breaks()<<" Total breaks: "<<the_tracks->get_num_breaks()<<endl;
    
    
}

void Track_point::set_blocks_v2()
{
    int i,curr_pillar, next_pillar, taxa, curr_block_start, min_end, dupl_index, num_links, break_count, num_complete=0;
    char link;
    Gene_Track_List_DX *my_pillar_site, *next_pillar_site;
    BOOL has_all_neighbors,  end_before;
    
    
    for (taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) {
        for(dupl_index=0; dupl_index<the_homologs->get_dupl_level(); dupl_index++) has_ends[taxa][dupl_index]=FALSE;
    }
    
    
    for (taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) {
        for(dupl_index=0; dupl_index<the_homologs->get_dupl_level(); dupl_index++) in_valid_break[taxa][dupl_index]=FALSE;
    }
    check_order();
    
    the_tracks->change_order(homolog_order);
    the_tracks->update_tracking();
    
    
    curr_pillar=0;
    curr_block_start=0;
    
    num_blocks=0;
    break_count=0;
    while(curr_pillar < (the_homologs->get_num_homologs()-1)) {
        next_pillar=curr_pillar+1;
        
        
        for(taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) {
            for(dupl_index=0; dupl_index<the_homologs->get_dupl_level(); dupl_index++) {
                if (((*the_tracks)[taxa].get_gene_track(curr_pillar, dupl_index)->my_locus !=0) &&
                    has_ends[taxa][dupl_index] ==TRUE)
                    break_count++;
                else {
                    if (((*the_tracks)[taxa].get_gene_track(curr_pillar, dupl_index)->my_locus !=0) &&
                        (*the_tracks)[taxa].has_forward_link(curr_pillar, dupl_index) ==FALSE)
                        has_ends[taxa][dupl_index]=TRUE;
                }
            }
        }
        
       // break_count=0;
       // for(taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) {
       //     for(dupl_index=0; dupl_index<the_homologs->get_dupl_level(); dupl_index++) {
       //         if (has_ends[taxa][dupl_index]== TRUE) break_count++;
       //     }
       // }
        
        //has_all_neighbors=TRUE;
        //taxa=0;
        //end_before=FALSE;
        
      
        if (break_count > max_breaks) {
            //if (has_all_neighbors==FALSE) {
            
            block_set[num_blocks].start=curr_block_start;
            //block_set[num_blocks].end=curr_pillar;
            block_set[num_blocks].end=curr_pillar-1;
            num_blocks++;
            curr_block_start=curr_pillar;
            break_count=0;
            
            for (taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) {
                for(dupl_index=0; dupl_index<the_homologs->get_dupl_level(); dupl_index++) has_ends[taxa][dupl_index]=FALSE;
            }
        }
        //cout<<"Pillar "<<curr_pillar<<" at block "<<num_blocks<<" has "<<break_count<<" breaks "<<" Complete: "<<num_complete<<endl;
        else
            curr_pillar++;
    }
    
    //Account for the final block
    block_set[num_blocks].start=curr_block_start;
    block_set[num_blocks].end=curr_pillar+1;
    num_blocks++;
    
    
    set_density();
    score_valid=TRUE;
    
    //cout<<" Setting "<<num_blocks<<" blocks and num_complete locations : "<<num_complete<<"\n";
    // cout<<"Setting "<<num_blocks<<" blocks. End is "<<block_set[num_blocks-1].end<<" Full Breaks: "<<the_tracks->get_num_full_breaks()<<" Total breaks: "<<the_tracks->get_num_breaks()<<endl;
   
    //for(i=0; i<num_blocks; i++) cout<<" BLOCK "<<i<<" S: "<<block_set[i].start<<" E: "<<block_set[i].end<<endl;
    
}

int Track_point::get_score()
{
    int part_score, taxa_index;
    
    if (score_valid == FALSE)
        //set_blocks();
        set_blocks_v2();
    
    if (curr_exchange->count_break_pos() ==FALSE) {
        if (curr_exchange->use_all_taxa() == TRUE)
            return(the_tracks->get_num_breaks());
        else {
            part_score=0;
            for(taxa_index=0; taxa_index<curr_exchange->get_num_taxa(); taxa_index++) {
                if (curr_exchange->is_taxa_used(taxa_index) == TRUE) {
                    part_score+=(*the_tracks)[taxa_index].count_num_breaks();
                    //cout<<"Score for "<<taxa_index<<" is "<<(*the_tracks)[taxa_index].count_num_breaks()<<endl;
                }
            }
            return(part_score);
        }
    }
    else
        //return(the_tracks->get_num_positions_w_breaks());
        return(num_blocks);
}

BOOL Track_point::check_order()
{
    int i, *new_order;
    BOOL retval=TRUE;
    
    new_order=new int [the_homologs->get_num_homologs()];
    
    for(i=0; i<the_homologs->get_num_homologs(); i++) new_order[i]=0;
    for(i=0; i<the_homologs->get_num_homologs(); i++) {
        new_order[homolog_order[i]]++;
        
        if (new_order[homolog_order[i]] >1) {
            retval=FALSE;
            cerr<<"ERROR: Found two entries for pillar "<<i<<" Current is "<<homolog_order[i]<<endl;
        }
    }
    
    for(i=0; i<the_homologs->get_num_homologs(); i++) {
        if (new_order[i]==0)  {
            cerr<<"ERROR: Never found an entry for pillar "<<i<<endl;
            retval=FALSE;
        }
    }
    
    delete[] new_order;
    return(retval);
    
}

void Track_point::set_density()
{
    int i, taxa, dupl_level;
    
    density_size=0;
    
    for(i=0; i<num_blocks; i++) {
        for(taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) {
            for(dupl_level=0; dupl_level<the_homologs->get_dupl_level(); dupl_level++) {
                if ((*the_tracks)[taxa].has_back_link(block_set[i].start, dupl_level) ==FALSE) {
                    location_density[density_size]=i;
                    density_size++;
                }
            }
        }
    }
}

Track_point::~Track_point ()
{
    int taxa;
    delete the_tracks;
    delete[] homolog_order;
    delete curr_exchange;
    delete[] block_set;
    
    for(taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) delete[] in_valid_break[taxa];
    delete[] in_valid_break;
}


void Track_point::initialize_pillar_nums()
{
    int taxa, i, dupl_level;
    
    the_tracks->change_order(homolog_order);
    the_tracks->update_tracking();
    
    for(i=0; i<the_homologs->get_num_homologs(); i++) {
        for(taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) {
            for(dupl_level=0; dupl_level<the_homologs->get_dupl_level(); dupl_level++) {
                if ((*the_tracks)[taxa].get_gene_track(i,dupl_level)->my_locus!=0)
                    (*the_tracks)[taxa].get_gene_track(i,dupl_level)->my_locus->get_gene_obj((*the_tracks)[taxa].get_gene_track(i, dupl_level)->index_num)->my_pillar=i;
                }
            }
        }
    
}

void Track_point::guess_order()
{
    int i, j, taxa, dupl_level, cnt, num_dupl_pillars=0, *dupl_pillar_ids, new_pillar1, new_pillar2, my_pillar, start_right1, start_right2, start_left1, start_left2, curr_loc=0, left_loc, right_loc, *local_order, offset, used=0, num_singl_pillars=0, *singl_pillar_ids, num_not_full;
    BOOL *used_dupl_pillars, *used_singl_pillars, *mark_pillars, is_singl, full, done=FALSE, mydir;
    queue <int> add_pillars;
    queue <BOOL> pillar_dir;
    
    initialize_pillar_nums();
    
    for(i=0; i<the_homologs->get_num_homologs(); i++) {
        full=TRUE;
        num_not_full=0;
        for(taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) {
            for(dupl_level=0; dupl_level<the_homologs->get_dupl_level(); dupl_level++) {
                if ((*the_tracks)[taxa].get_gene_track(i,dupl_level)->my_locus==0) {
                    full=FALSE;
                    num_not_full++;
                }
            }
        }
        if (full ==TRUE) num_dupl_pillars++;
        //if (num_not_full<2) num_dupl_pillars++;
    }
    
    search_array_size=2*the_genomes->get_num_genomes() * the_homologs->get_dupl_level();
    
    cutoff=the_genomes->get_num_genomes()-3;
    pillar_scores=new int [search_array_size];
    pillar_lookup=new int [search_array_size];
    pillar_neighbor_dirs=new int * [search_array_size];
    dupl_pillar_ids=new int [num_dupl_pillars];
    used_dupl_pillars=new BOOL [num_dupl_pillars];
    used_pillars=new BOOL[the_homologs->get_num_homologs()];
    local_order=new int[the_homologs->get_num_homologs()];
    mark_pillars=new BOOL [the_homologs->get_num_homologs()];
    
    
    for(i=0; i<search_array_size; i++)
        pillar_neighbor_dirs[i]=new int [the_genomes->get_num_genomes()];
    
    
    for(i=0; i<the_homologs->get_num_homologs(); i++)
        used_pillars[i]=FALSE;
    
    
    cnt=0;
    for(i=0; i<the_homologs->get_num_homologs(); i++) {
        full=TRUE;
        num_not_full=0;
        for(taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) {
            for(dupl_level=0; dupl_level<the_homologs->get_dupl_level(); dupl_level++) {
                if ((*the_tracks)[taxa].get_gene_track(i,dupl_level)->my_locus==0) {
                    full=FALSE;
                    num_not_full++;
                }
            }
        }
        //if (num_not_full<2) {
        if (full ==TRUE){
            dupl_pillar_ids[cnt]=i;
            used_dupl_pillars[cnt]=FALSE;
            cnt++;
        }
    }
    while (done==FALSE) {
        for(i=0; i<num_dupl_pillars; i++) {
            if (used_pillars[dupl_pillar_ids[i]]==TRUE)
                used_dupl_pillars[i]=TRUE;
        }
        
        i=0;
        
        while ((i<num_dupl_pillars) && (used_dupl_pillars[i]==TRUE)) i++;
        
        if (i==num_dupl_pillars) done=TRUE;
        else {
            for(j=0; j<the_homologs->get_num_homologs(); j++) mark_pillars[j]=FALSE;
            
            
            my_pillar=dupl_pillar_ids[i];
            mark_pillars[dupl_pillar_ids[i]] =TRUE;
            used_pillars[dupl_pillar_ids[i]]=TRUE;
            used_dupl_pillars[i]=TRUE;
            used++;
            cout<<"Initialing block from pillar "<<dupl_pillar_ids[i]<<endl;
            
           
            score_neighbors(my_pillar);
            start_right1=get_best_pillar();
            
           
            
            right_loc=left_loc=curr_loc;
            
            local_order[my_pillar]=curr_loc;
            
            if ((start_right1!= -1) && (pillar_scores[start_right1] >= cutoff)) {
                right_loc=curr_loc+1;
                
               // cout<<"R1 is good: "<<pillar_lookup[start_right1]<<endl;
                
                add_pillars.push(pillar_lookup[start_right1]);
                pillar_dir.push(TRUE);
                local_order[pillar_lookup[start_right1]]=right_loc;
                mark_pillars[pillar_lookup[start_right1]]=TRUE;
                used_pillars[pillar_lookup[start_right1]]=TRUE;
                used++;
                
                start_left1=get_compl_pillar(start_right1);
                if ((start_left1 != -1) && (pillar_scores[start_left1] >= cutoff) && (used_pillars[pillar_lookup[start_left1]]==FALSE)) {
                 //   cout<<"L1 (cR1) is good: "<<pillar_lookup[start_left1]<<endl;
                    
                    add_pillars.push(pillar_lookup[start_left1]);
                    pillar_dir.push(FALSE);
                    left_loc=curr_loc-1;
                    local_order[pillar_lookup[start_left1]]=left_loc;
                    mark_pillars[pillar_lookup[start_left1]]=TRUE;
                    used_pillars[pillar_lookup[start_left1]]=TRUE;
                    used++;
                }
                
                start_right2=get_best_pillar();
                if ((start_right2!= -1) && (pillar_scores[start_right2] >= cutoff)) {
                    
                   // cout<<"R2 is good: "<<pillar_lookup[start_right2]<<endl;
                    
                    
                    right_loc=right_loc+1;
                    add_pillars.push(pillar_lookup[start_right2]);
                    pillar_dir.push(TRUE);
                    local_order[pillar_lookup[start_right2]]=right_loc;
                    mark_pillars[pillar_lookup[start_right2]]=TRUE;
                    used_pillars[pillar_lookup[start_right2]]=TRUE;
                    used++;
                    
                    start_left2=get_compl_pillar(start_right2);
                    if ((start_left2 != -1) && (pillar_scores[start_left2] >= cutoff) && (used_pillars[pillar_lookup[start_left2]]==FALSE)) {
                        left_loc=left_loc-1;
                     //   cout<<"L2 (cR2) is good: "<<pillar_lookup[start_left2]<<endl;
                        add_pillars.push(pillar_lookup[start_left2]);
                        pillar_dir.push(FALSE);
                        local_order[pillar_lookup[start_left2]]=left_loc;
                        mark_pillars[pillar_lookup[start_left2]]=TRUE;
                        used_pillars[pillar_lookup[start_left2]]=TRUE;
                        used++;
                    }
                }
            }
            
            if ((start_left1==-1) && (start_left2==-1)) {
                start_left1=get_best_pillar();
                
                if ((start_left1 != -1) && (pillar_scores[start_left1] >= cutoff)) {
                    left_loc=curr_loc-1;
                    //cout<<"L1 (raw) is good: "<<pillar_lookup[start_left1]<<endl;
                    add_pillars.push(pillar_lookup[start_left1]);
                    pillar_dir.push(FALSE);
                    
                    local_order[pillar_lookup[start_left1]]=left_loc;
                    mark_pillars[pillar_lookup[start_left1]]=TRUE;
                    used_pillars[pillar_lookup[start_left1]]=TRUE;
                    used++;
                    
                    start_left2=get_best_pillar();
                    if ((start_left2 != -1) && (pillar_scores[start_left2] >= cutoff)) {
                        left_loc=left_loc-1;
                        
                      //  cout<<"L2 (raw) is good: "<<pillar_lookup[start_left2]<<endl;
                        add_pillars.push(pillar_lookup[start_left2]);
                        pillar_dir.push(FALSE);
                        
                        local_order[pillar_lookup[start_left2]]=left_loc;
                        mark_pillars[pillar_lookup[start_left2]]=TRUE;
                        used_pillars[pillar_lookup[start_left2]]=TRUE;
                        used++;
                    }
                }
                
            }
            
           
            while(!add_pillars.empty()) {
                my_pillar=add_pillars.front();
                add_pillars.pop();
                mydir=pillar_dir.front();
                pillar_dir.pop();
                
                
                //cout<<"Processing "<<my_pillar<<" Right: "<<mydir<<endl;
                
                score_neighbors(my_pillar);
                new_pillar1=get_best_pillar();
                
                if ((new_pillar1!= -1) && (pillar_scores[new_pillar1] >= cutoff)) {
                    add_pillars.push(pillar_lookup[new_pillar1]);
                    mark_pillars[pillar_lookup[new_pillar1]]=TRUE;
                    used_pillars[pillar_lookup[new_pillar1]]=TRUE;
                    used++;
                    if (mydir ==TRUE) {
                        pillar_dir.push(TRUE);
                        right_loc=right_loc+1;
                        local_order[pillar_lookup[new_pillar1]]=right_loc;
                    }
                    else {
                        pillar_dir.push(FALSE);
                        left_loc=left_loc-1;
                        local_order[pillar_lookup[new_pillar1]]=left_loc;
                    }
                    
                    new_pillar2=get_best_pillar();
                    
                    if ((new_pillar2!= -1) && (pillar_scores[new_pillar2] >= cutoff)) {
                        add_pillars.push(pillar_lookup[new_pillar2]);
                        mark_pillars[pillar_lookup[new_pillar2]]=TRUE;
                        used_pillars[pillar_lookup[new_pillar2]]=TRUE;
                        used++;
                        if (mydir == TRUE) {
                            pillar_dir.push(TRUE);
                            right_loc=right_loc+1;
                            local_order[pillar_lookup[new_pillar2]]=right_loc;
                            
                        }
                        else {
                            pillar_dir.push(FALSE);
                            left_loc=left_loc-1;
                            local_order[pillar_lookup[new_pillar2]]=left_loc;
                        }
                        
                    }
                }
            }
            
            cout<<"Added: "<<right_loc-curr_loc<<" pillars to the right. Used: "<<used<<"\n";
            
           
            cout<<"Added: "<<curr_loc-left_loc<<" pillars to the left\n";
            offset=curr_loc-left_loc;
            cout<<"reseting indices by "<<offset<<": ";
            for(j=0; j<the_homologs->get_num_homologs(); j++) {
                //local_order[j]+=offset;
                if (mark_pillars[j] == TRUE) {
                    homolog_order[local_order[j]+offset]=j;
                    cout<<j<<"\t";
                  //  cout<<"Assigning "<<j<<" to position "<<local_order[j]+offset<<endl;
                    //used_pillars[j]=TRUE;
                }
            }
            
            cout<<"\nMade block of size "<<right_loc-left_loc+1<<endl;
            curr_loc+=right_loc-left_loc+1;
            cout<<"Now at location "<<curr_loc<<" used: "<<used<<endl;
            done=TRUE;
            
            for(i=0; i<num_dupl_pillars; i++){
                if (used_dupl_pillars[i]==FALSE) done=FALSE;
            }
        }
    }
    
    
#if 1
    num_dupl_pillars=0;
    for(i=0; i<the_homologs->get_num_homologs(); i++) {
        num_not_full=0;
        for(taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) {
            for(dupl_level=0; dupl_level<the_homologs->get_dupl_level(); dupl_level++) {
                if ((*the_tracks)[taxa].get_gene_track(i,dupl_level)->my_locus==0)
                    num_not_full++;
            }
        }
        if ((num_not_full<2) && (used_pillars[i] ==FALSE)) num_dupl_pillars++;
    }
    
    delete[] dupl_pillar_ids;
    delete[] used_dupl_pillars;
    
    dupl_pillar_ids=new int[num_dupl_pillars];
    used_dupl_pillars=new BOOL[num_dupl_pillars];
   
    cout<<"Initializing 2 from "<<num_dupl_pillars<<" mostly duplicated pillars\n";
    cnt=0;
    for(i=0; i<the_homologs->get_num_homologs(); i++) {
        num_not_full=0;
        for(taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) {
            for(dupl_level=0; dupl_level<the_homologs->get_dupl_level(); dupl_level++) {
                if ((*the_tracks)[taxa].get_gene_track(i,dupl_level)->my_locus==0)
                    num_not_full++;
            }
        }
        if ((num_not_full<2) && (used_pillars[i] ==FALSE)){
            dupl_pillar_ids[cnt]=i;
            used_dupl_pillars[cnt]=FALSE;
            cnt++;
        }
    }
    
    done =FALSE;
    
    while (done==FALSE) {
        for(i=0; i<num_dupl_pillars; i++) {
            if (used_pillars[dupl_pillar_ids[i]]==TRUE)
                used_dupl_pillars[i]=TRUE;
        }
        
        i=0;
        
        while ((i<num_dupl_pillars) && (used_dupl_pillars[i]==TRUE)) i++;
        
        if (i==num_dupl_pillars) done=TRUE;
        else {
            for(j=0; j<the_homologs->get_num_homologs(); j++) mark_pillars[j]=FALSE;
            
            
            my_pillar=dupl_pillar_ids[i];
            mark_pillars[dupl_pillar_ids[i]] =TRUE;
            used_pillars[dupl_pillar_ids[i]]=TRUE;
            used_dupl_pillars[i]=TRUE;
            used++;
            cout<<"Initialing -1 block from pillar "<<dupl_pillar_ids[i]<<endl;
            
            
            score_neighbors(my_pillar);
            start_right1=get_best_pillar();
            
            
            
            right_loc=left_loc=curr_loc;
            
            local_order[my_pillar]=curr_loc;
            
            if ((start_right1!= -1) && (pillar_scores[start_right1] >= cutoff)) {
                right_loc=curr_loc+1;
                
                // cout<<"R1 is good: "<<pillar_lookup[start_right1]<<endl;
                
                add_pillars.push(pillar_lookup[start_right1]);
                pillar_dir.push(TRUE);
                local_order[pillar_lookup[start_right1]]=right_loc;
                mark_pillars[pillar_lookup[start_right1]]=TRUE;
                used_pillars[pillar_lookup[start_right1]]=TRUE;
                used++;
                
                start_left1=get_compl_pillar(start_right1);
                if ((start_left1 != -1) && (pillar_scores[start_left1] >= cutoff) && (used_pillars[pillar_lookup[start_left1]]==FALSE)) {
                    //   cout<<"L1 (cR1) is good: "<<pillar_lookup[start_left1]<<endl;
                    
                    add_pillars.push(pillar_lookup[start_left1]);
                    pillar_dir.push(FALSE);
                    left_loc=curr_loc-1;
                    local_order[pillar_lookup[start_left1]]=left_loc;
                    mark_pillars[pillar_lookup[start_left1]]=TRUE;
                    used_pillars[pillar_lookup[start_left1]]=TRUE;
                    used++;
                }
                
                start_right2=get_best_pillar();
                if ((start_right2!= -1) && (pillar_scores[start_right2] >= cutoff)) {
                    
                    // cout<<"R2 is good: "<<pillar_lookup[start_right2]<<endl;
                    
                    
                    right_loc=right_loc+1;
                    add_pillars.push(pillar_lookup[start_right2]);
                    pillar_dir.push(TRUE);
                    local_order[pillar_lookup[start_right2]]=right_loc;
                    mark_pillars[pillar_lookup[start_right2]]=TRUE;
                    used_pillars[pillar_lookup[start_right2]]=TRUE;
                    used++;
                    
                    start_left2=get_compl_pillar(start_right2);
                    if ((start_left2 != -1) && (pillar_scores[start_left2] >= cutoff) && (used_pillars[pillar_lookup[start_left2]]==FALSE)) {
                        left_loc=left_loc-1;
                        //   cout<<"L2 (cR2) is good: "<<pillar_lookup[start_left2]<<endl;
                        add_pillars.push(pillar_lookup[start_left2]);
                        pillar_dir.push(FALSE);
                        local_order[pillar_lookup[start_left2]]=left_loc;
                        mark_pillars[pillar_lookup[start_left2]]=TRUE;
                        used_pillars[pillar_lookup[start_left2]]=TRUE;
                        used++;
                    }
                }
            }
            
            if ((start_left1==-1) && (start_left2==-1)) {
                start_left1=get_best_pillar();
                
                if ((start_left1 != -1) && (pillar_scores[start_left1] >= cutoff)) {
                    left_loc=curr_loc-1;
                    //cout<<"L1 (raw) is good: "<<pillar_lookup[start_left1]<<endl;
                    add_pillars.push(pillar_lookup[start_left1]);
                    pillar_dir.push(FALSE);
                    
                    local_order[pillar_lookup[start_left1]]=left_loc;
                    mark_pillars[pillar_lookup[start_left1]]=TRUE;
                    used_pillars[pillar_lookup[start_left1]]=TRUE;
                    used++;
                    
                    start_left2=get_best_pillar();
                    if ((start_left2 != -1) && (pillar_scores[start_left2] >= cutoff)) {
                        left_loc=left_loc-1;
                        
                        //  cout<<"L2 (raw) is good: "<<pillar_lookup[start_left2]<<endl;
                        add_pillars.push(pillar_lookup[start_left2]);
                        pillar_dir.push(FALSE);
                        
                        local_order[pillar_lookup[start_left2]]=left_loc;
                        mark_pillars[pillar_lookup[start_left2]]=TRUE;
                        used_pillars[pillar_lookup[start_left2]]=TRUE;
                        used++;
                    }
                }
                
            }
            
            
            while(!add_pillars.empty()) {
                my_pillar=add_pillars.front();
                add_pillars.pop();
                mydir=pillar_dir.front();
                pillar_dir.pop();
                
                
                //cout<<"Processing "<<my_pillar<<" Right: "<<mydir<<endl;
                
                score_neighbors(my_pillar);
                new_pillar1=get_best_pillar();
                
                if ((new_pillar1!= -1) && (pillar_scores[new_pillar1] >= cutoff)) {
                    add_pillars.push(pillar_lookup[new_pillar1]);
                    mark_pillars[pillar_lookup[new_pillar1]]=TRUE;
                    used_pillars[pillar_lookup[new_pillar1]]=TRUE;
                    used++;
                    if (mydir ==TRUE) {
                        pillar_dir.push(TRUE);
                        right_loc=right_loc+1;
                        local_order[pillar_lookup[new_pillar1]]=right_loc;
                    }
                    else {
                        pillar_dir.push(FALSE);
                        left_loc=left_loc-1;
                        local_order[pillar_lookup[new_pillar1]]=left_loc;
                    }
                    
                    new_pillar2=get_best_pillar();
                    
                    if ((new_pillar2!= -1) && (pillar_scores[new_pillar2] >= cutoff)) {
                        add_pillars.push(pillar_lookup[new_pillar2]);
                        mark_pillars[pillar_lookup[new_pillar2]]=TRUE;
                        used_pillars[pillar_lookup[new_pillar2]]=TRUE;
                        used++;
                        if (mydir == TRUE) {
                            pillar_dir.push(TRUE);
                            right_loc=right_loc+1;
                            local_order[pillar_lookup[new_pillar2]]=right_loc;
                            
                        }
                        else {
                            pillar_dir.push(FALSE);
                            left_loc=left_loc-1;
                            local_order[pillar_lookup[new_pillar2]]=left_loc;
                        }
                        
                    }
                }
            }
            
            cout<<"Added: "<<right_loc-curr_loc<<" pillars to the right. Used: "<<used<<"\n";
            
            
            cout<<"Added: "<<curr_loc-left_loc<<" pillars to the left\n";
            offset=curr_loc-left_loc;
            cout<<"reseting indices by "<<offset<<": ";
            for(j=0; j<the_homologs->get_num_homologs(); j++) {
                //local_order[j]+=offset;
                if (mark_pillars[j] == TRUE) {
                    homolog_order[local_order[j]+offset]=j;
                    cout<<j<<"\t";
                    //  cout<<"Assigning "<<j<<" to position "<<local_order[j]+offset<<endl;
                    //used_pillars[j]=TRUE;
                }
            }
            
            cout<<"\nMade block of size "<<right_loc-left_loc+1<<endl;
            curr_loc+=right_loc-left_loc+1;
            cout<<"Now at location "<<curr_loc<<" used: "<<used<<endl;
            done=TRUE;
            
            for(i=0; i<num_dupl_pillars; i++){
                if (used_dupl_pillars[i]==FALSE) done=FALSE;
            }
        }
    }
    
    
#endif
    
    for(i=0; i<the_homologs->get_num_homologs(); i++) {
        is_singl=TRUE;
        for(taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) {
                if (((*the_tracks)[taxa].get_gene_track(i,0)->my_locus!=0)  &&
                    ((*the_tracks)[taxa].get_gene_track(i,1)->my_locus!=0)) is_singl=FALSE;
        }
        if ((is_singl ==TRUE) && (used_pillars[i]==FALSE)) num_singl_pillars++;
    }
            
            
#if 1
    cout<<"Initializing from "<<num_singl_pillars<<" singl pillars"<<" and "<<curr_loc<<endl;
    
    used_singl_pillars=new BOOL [num_singl_pillars];
    singl_pillar_ids=new int [num_singl_pillars];
    
    for(i=0; i<num_singl_pillars; i++) used_singl_pillars[i]=FALSE;
    cnt=0;
    
    for(i=0; i<the_homologs->get_num_homologs(); i++) {
        is_singl=TRUE;
        for(taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) {
            if (((*the_tracks)[taxa].get_gene_track(i,0)->my_locus!=0)  &&
                ((*the_tracks)[taxa].get_gene_track(i,1)->my_locus!=0)) is_singl=FALSE;
        }
        if ((is_singl ==TRUE) && (used_pillars[i]==FALSE)) {
            singl_pillar_ids[cnt]=i;
            cnt++;
        }
    }
    

    
    done=FALSE;
    
    while(done == FALSE) {
        for(i=0; i<num_singl_pillars; i++) {
            if (used_pillars[singl_pillar_ids[i]]==TRUE)
                used_singl_pillars[i]=TRUE;
        }
        
        i=0;
        
        while ((i<num_singl_pillars) && (used_singl_pillars[i]==TRUE)) i++;
        
        if (i==num_singl_pillars) done=TRUE;
        else {
            for(j=0; j<the_homologs->get_num_homologs(); j++) mark_pillars[j]=FALSE;
            
            my_pillar=singl_pillar_ids[i];
            mark_pillars[singl_pillar_ids[i]] =TRUE;
            used_pillars[singl_pillar_ids[i]]=TRUE;
            used_singl_pillars[i]=TRUE;
            used++;
            cout<<"Initialing single block from pillar "<<singl_pillar_ids[i]<<endl;
            
            
            score_neighbors(my_pillar);
            start_right1=get_best_pillar();
            
            start_left1 = -1;
            right_loc=curr_loc;
            left_loc=curr_loc;
            
            local_order[my_pillar]=curr_loc;
            
            if ((start_right1!= -1) && (pillar_scores[start_right1] >= cutoff)) {
                right_loc=curr_loc+1;
                
                // cout<<"R1 is good: "<<pillar_lookup[start_right1]<<endl;
                
                add_pillars.push(pillar_lookup[start_right1]);
                pillar_dir.push(TRUE);
                local_order[pillar_lookup[start_right1]]=right_loc;
                mark_pillars[pillar_lookup[start_right1]]=TRUE;
                used_pillars[pillar_lookup[start_right1]]=TRUE;
                used++;
                
                start_left1=get_compl_pillar(start_right1);
                if ((start_left1 != -1) && (pillar_scores[start_left1] >= cutoff) && (used_pillars[pillar_lookup[start_left1]]==FALSE)) {
                    //   cout<<"L1 (cR1) is good: "<<pillar_lookup[start_left1]<<endl;
                    
                    add_pillars.push(pillar_lookup[start_left1]);
                    pillar_dir.push(FALSE);
                    left_loc=curr_loc-1;
                    local_order[pillar_lookup[start_left1]]=left_loc;
                    mark_pillars[pillar_lookup[start_left1]]=TRUE;
                    used_pillars[pillar_lookup[start_left1]]=TRUE;
                    used++;
                }
                
            }
            
            if (start_left1==-1) {
                start_left1=get_best_pillar();
                
                if ((start_left1 != -1) && (pillar_scores[start_left1] >= cutoff)) {
                    left_loc=curr_loc-1;
                    //cout<<"L1 (raw) is good: "<<pillar_lookup[start_left1]<<endl;
                    add_pillars.push(pillar_lookup[start_left1]);
                    pillar_dir.push(FALSE);
                    
                    local_order[pillar_lookup[start_left1]]=left_loc;
                    mark_pillars[pillar_lookup[start_left1]]=TRUE;
                    used_pillars[pillar_lookup[start_left1]]=TRUE;
                    used++;
                }
            }
            
            
            while(!add_pillars.empty()) {
                my_pillar=add_pillars.front();
                add_pillars.pop();
                mydir=pillar_dir.front();
                pillar_dir.pop();
                
                score_neighbors(my_pillar);
                new_pillar1=get_best_pillar();
                
                if ((new_pillar1!= -1) && (pillar_scores[new_pillar1] >= cutoff)) {
                    add_pillars.push(pillar_lookup[new_pillar1]);
                    mark_pillars[pillar_lookup[new_pillar1]]=TRUE;
                    used_pillars[pillar_lookup[new_pillar1]]=TRUE;
                    used++;
                    if (mydir ==TRUE) {
                        pillar_dir.push(TRUE);
                        right_loc=right_loc+1;
                        local_order[pillar_lookup[new_pillar1]]=right_loc;
                    }
                    else {
                        pillar_dir.push(FALSE);
                        left_loc=left_loc-1;
                        local_order[pillar_lookup[new_pillar1]]=left_loc;
                    }
                    
                    new_pillar2=get_best_pillar();
                    
                    if ((new_pillar2!= -1) && (pillar_scores[new_pillar2] >= cutoff)) {
                        add_pillars.push(pillar_lookup[new_pillar2]);
                        mark_pillars[pillar_lookup[new_pillar2]]=TRUE;
                        used_pillars[pillar_lookup[new_pillar2]]=TRUE;
                        used++;
                        if (mydir == TRUE) {
                            pillar_dir.push(TRUE);
                            right_loc=right_loc+1;
                            local_order[pillar_lookup[new_pillar2]]=right_loc;
                        }
                        else {
                            pillar_dir.push(FALSE);
                            left_loc=left_loc-1;
                            local_order[pillar_lookup[new_pillar2]]=left_loc;
                        }
                    }
                }
            }
            
            //cout<<"Added: "<<right_loc-curr_loc<<" pillars to the right. Used: "<<used<<"\n";
            
            cout<<"Made singl block of size "<<right_loc-left_loc+1<<": ";
            //cout<<"Added: "<<curr_loc-left_loc<<" pillars to the left\n";
            offset=curr_loc-left_loc;
            //cout<<"reseting indices by "<<offset<<endl;
            for(j=0; j<the_homologs->get_num_homologs(); j++) {
                //local_order[j]+=offset;
                if (mark_pillars[j] == TRUE) {
                    homolog_order[local_order[j]+offset]=j;
                    cout<<j<<"\t";
                    //  cout<<"Assigning "<<j<<" to position "<<local_order[j]+offset<<endl;
                    //used_pillars[j]=TRUE;
                }
                
            }
            cout<<endl;
           // cout<<"Made singl block of size "<<right_loc-left_loc+1<<endl;
            curr_loc+=right_loc-left_loc+1;
            cout<<"Now at location "<<curr_loc<<" used: "<<used<<endl;
            done=TRUE;
            
            for(i=0; i<num_singl_pillars; i++){
                if (used_singl_pillars[i]==FALSE) done=FALSE;
            }
        }
        
    }
#endif
    for(i=0; i<the_homologs->get_num_homologs(); i++)
    {
        if (used_pillars[i] == FALSE) {
            homolog_order[curr_loc]=i;
            curr_loc++;
        }
    }
    cout<<"Added old order to "<<curr_loc<<endl;
    
    
    
    //cout<<"Guessed order:\n";
    //for(i=0; i<the_homologs->get_num_homologs(); i++)
      //  cout<<i<<"\t"<<homolog_order[i]<<endl;
    //cout<<"Running check"<<endl<<flush;
    check_order();
    
    delete[] pillar_scores;
    delete[] pillar_lookup;
    delete[] dupl_pillar_ids;
    delete[] used_dupl_pillars;
    delete[] used_pillars;
    delete[] local_order;
    delete[] mark_pillars;
        delete[] singl_pillar_ids;
        delete[] used_singl_pillars;
    
    for(i=0; i<search_array_size; i++)
        delete [] pillar_neighbor_dirs[i];
    delete[] pillar_neighbor_dirs;
}



void Track_point::score_neighbors(int current_pillar)
{
    int i, taxa, dupl_level, myloc;
    Gene *my_neighbor;
    
    for(i=0; i<search_array_size; i++ ) {
        pillar_lookup[i]=-1;
        pillar_scores[i]=0;
    }
    
    for(taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) {
        for(dupl_level=0; dupl_level<the_homologs->get_dupl_level(); dupl_level++) {
            if ((*the_tracks)[taxa].get_gene_track(current_pillar,dupl_level)->my_locus !=0) {
                my_neighbor=(*the_tracks)[taxa].get_gene_track(current_pillar,dupl_level)->my_locus->get_gene_obj((*the_tracks)[taxa].get_gene_track(current_pillar, dupl_level)->index_num)->get_neighbor(0);
                if (my_neighbor !=0) {
                    myloc = get_pillar_search_num(my_neighbor->my_pillar);
                    pillar_scores[myloc]++;
                    pillar_neighbor_dirs[myloc][taxa]=0;
                }
                
                my_neighbor=(*the_tracks)[taxa].get_gene_track(current_pillar,dupl_level)->my_locus->get_gene_obj((*the_tracks)[taxa].get_gene_track(current_pillar, dupl_level)->index_num)->get_neighbor(1);
                if (my_neighbor !=0) {
                    myloc = get_pillar_search_num(my_neighbor->my_pillar);
                    pillar_scores[myloc]++;
                    pillar_neighbor_dirs[myloc][taxa]=1;
                }
                
            }
        }
    }
    
}



int Track_point::get_pillar_search_num(int pillar_num)
{
    int ret_val=0;
    
    while ((pillar_lookup[ret_val] != -1) && (pillar_lookup[ret_val] != pillar_num)) ret_val++;

    if (pillar_lookup[ret_val] == -1) pillar_lookup[ret_val]=pillar_num;
    return(ret_val);
}


int Track_point::get_best_pillar()
{
    int myloc, max=-1, max_loc=-1;
    
    
    for(myloc=0; myloc<search_array_size; myloc++)
    {
        
        //cout<<"BEST SEARCH: "<<myloc<<": "<<pillar_lookup[myloc]<<" and used "<<used_pillars[pillar_lookup[myloc]]<<" S= "<<pillar_scores[myloc]<<endl;
        if ((pillar_lookup[myloc]!=-1) && (used_pillars[pillar_lookup[myloc]]==FALSE)) {
            if (max <pillar_scores[myloc])
            {
                max=pillar_scores[myloc];
                max_loc=myloc;
            }
        }
    }
   // cout<<"Best score "<<max<<" at "<<max_loc<<endl;
    if (max_loc == -1) return(-1);
    else
        return(max_loc);
    
}


int Track_point::get_compl_pillar(int pillar_lookup_id)
{
    int myloc=0, taxa;
    BOOL complement, found=FALSE;
    
    while((found == FALSE) && ( myloc<search_array_size))
    {
        complement=TRUE;
        
        if(pillar_lookup[myloc]!=-1) {
            for(taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) {
                if (pillar_neighbor_dirs[myloc][taxa] != (1 ^ pillar_neighbor_dirs[pillar_lookup_id][taxa]))
                    complement=FALSE;
            }
        }
        if (complement ==TRUE) found=TRUE;
        else
            myloc++;
    }
    
    if (found == FALSE) return(-1);
    else return(myloc);
    
}

Track_anneal::Track_anneal(Exchange *cexchange, int nwalkers) : Anneal<Track_point>(cexchange, nwalkers)
{
    int i, j, k, l;
    string ***orthos;
    Genome **new_genomes;
    Clade *extra_genomes;
#ifdef _OPEN_MP_VERSION_
    cout<<"Allocating "<<cexchange->get_num_open_mp_threads()<<" extra points for parallel version\n";
    pre_scores =new double [cexchange->get_num_open_mp_threads()];
    test_points=new Track_point * [cexchange->get_num_open_mp_threads()];
    
    
    #pragma omp parallel for private (i,j, k,l, new_genomes, extra_genomes, orthos)
    for(i=0; i<cexchange->get_num_open_mp_threads(); i++) {
        test_points[i]= new Track_point(cexchange);
        new_genomes = new Genome*[global_genomes->get_num_genomes()];
        
        for (j=0; j<global_genomes->get_num_genomes(); j++) {
            new_genomes[j] = new Genome;
            (*new_genomes[j])=(*global_genomes)[j];
        }
        extra_genomes=new Clade(global_genomes->get_num_genomes(), new_genomes);
        test_points[i]->the_genomes=extra_genomes;
        
        orthos=new string ** [global_homologs->get_num_homologs()];
        
        for(j=0; j<global_homologs->get_num_homologs(); j++) {
            orthos[j]=new string * [test_points[i]->the_genomes->get_num_genomes()];
            for(k=0; k<test_points[i]->the_genomes->get_num_genomes(); k++) {
                orthos[j][k]=new string [global_homologs->get_dupl_level()];
            }
        }

        for(j=0; j<global_homologs->get_num_homologs(); j++) {
            for(k=0; k<test_points[i]->the_genomes->get_num_genomes(); k++) {
                for(l=0; l<global_homologs->get_dupl_level(); l++) {
                    if ((*global_homologs)[j][k].get_gene_obj(l) != 0)
                        orthos[j][k][l]=(*global_homologs)[j][k].get_gene_obj(l)->get_name_string();
                    else
                        orthos[j][k][l]="NONE";
                }
            }
        }
        
        test_points[i]->the_homologs=new WGX_Data(global_homologs->get_num_homologs(), global_homologs->get_dupl_level(), orthos, test_points[i]->the_genomes);
        
        test_points[i]->the_tracks = new WGX_Tracks(test_points[i]->the_homologs, test_points[i]->the_genomes);
        
        for(j=0; j<global_homologs->get_num_homologs(); j++) {
            for(k=0; k<test_points[i]->the_genomes->get_num_genomes(); k++)
                delete[] orthos[j][k];
            delete[] orthos[j];
        }
        delete[] orthos;
        
    }

#endif

}


void Track_anneal::move(int walker)
{
    int mythread, i, j,  dscore, old,  block_start, block_length, new_block_loc, start, length, new_loc, *new_order, first, second, max, max_loc;
    double choice, val, lambda;
    BOOL prev, do_reverse=FALSE, move_far;
    Track_point *the_point;



    *pre_move=*(*current_walkers)[walker];
#ifdef _OPEN_MP_VERSION_
#pragma omp parallel for private (i)
    for(i=0; i<omp_get_num_threads(); i++)
        (*test_points[i])=(*(*current_walkers)[walker]->point);
    
#else
    the_point=(*current_walkers)[walker]->point;
#endif

    //the_point->the_tracks->print_all_tracks(FALSE);

    //if (the_point->the_tracks!=0) delete the_point->the_tracks;

#ifdef _OPEN_MP_VERSION_
#pragma omp parallel for private (i,j,mythread, block_start, choice, move_far, block_length, start, length, new_block_loc, new_loc, do_reverse, new_order, the_point)
for(mythread=0; mythread<omp_get_num_threads(); mythread++) {
    the_point=test_points[mythread];
#endif
    block_start =choose_move_start(walker);
  
#pragma omp critical
    choice=ranf();
    move_far=FALSE;
    
    if (choice <0.15) {
#pragma omp critical
        block_length = (int)gennor(40, 10);
        move_far=TRUE;
    }
    else {
        if (choice> 0.5)
#pragma omp critical
            block_length = (int)gennor(6, 3);
        else
            block_length=1;
    }
    
    while ((block_length < 1) || ((block_start+block_length) >= the_point->get_num_blocks())) {
        #pragma omp critical
        block_start =ignuin(0, the_point->get_num_blocks()-1);
        #pragma omp critical
        block_length = (int)gennor(5, 2);
    }
    
    start=the_point->get_block_start(block_start);
    length = the_point->get_block_end(block_start+(block_length-1))-the_point->get_block_start(block_start)+1;
    
    if ((ranf()<0.35) || (move_far == TRUE)) {
        //cout<<"Moving far\n";
        
        new_block_loc=choose_far_loc(walker, block_start, block_length);
        
        new_loc=the_point->get_block_start(new_block_loc);
        
        //if (new_loc >= start) new_loc+=length;
    
    }
    else {
        //cout<<"Moving near: "<<block_start<<"-->"<<block_length<<" ";
        #pragma omp critical
        new_block_loc=gennor(0,3)+block_start;
        while(((new_block_loc >=block_start) && (new_block_loc < (block_start+block_length))) || (new_block_loc <0) || (new_block_loc >= the_point->get_num_blocks()))
            #pragma omp critical
            new_block_loc=gennor(0,3)+block_start;
        
        new_loc=the_point->get_block_start(new_block_loc);
        
    }
    
    
    #pragma omp critical
    if (ranf() > 0.3) do_reverse=TRUE;
    
   //if (do_reverse==FALSE) cout<<"Moving forward "<<start<<" for "<<length<<" to "<<new_loc<"\t";
   //else cout<<"Moving reverse "<<start<<" for "<<length<<" to "<<new_loc<<"\t";
    
    new_order=new int[the_point->the_homologs->get_num_homologs()];
    
    if (start > new_loc) {
        for(i=0; i<new_loc; i++) new_order[i]=the_point->homolog_order[i];
        if (do_reverse==FALSE) {
            for(i=new_loc; i<new_loc+length; i++) new_order[i] = the_point->homolog_order[start+(i-new_loc)];
        }
        else {
            for(i=new_loc; i<new_loc+length; i++) new_order[i] = the_point->homolog_order[(start+length-1)-(i-new_loc)];
        }
        for(i=new_loc+length; i<start+length; i++) new_order[i] = the_point->homolog_order[i-length];
        for(i=start+length; i<the_point->the_homologs->get_num_homologs(); i++)
            new_order[i]=the_point->homolog_order[i];
    }
    else {
        for(i=0; i<start; i++) new_order[i]=the_point->homolog_order[i];
        for(i=start+length; i<new_loc; i++) new_order[i-length]=the_point->homolog_order[i];
        if (do_reverse==FALSE) {
            for(i=new_loc; i<new_loc+length; i++) new_order[i-length]=the_point->homolog_order[start+(i-new_loc)];
        }
        else {
            for(i=new_loc; i<new_loc+length; i++) new_order[i-length]=the_point->homolog_order[(start+length-1)-(i-new_loc)];
        }
        for(i=new_loc; i<the_point->the_homologs->get_num_homologs(); i++)
                new_order[i]=the_point->homolog_order[i];
        
    }
    
    //for(i=0; i<the_point->the_homologs->get_num_homologs(); i++) cout<<new_order[i]<<" ";
    //cout<<endl;
    
    for(i=0; i<the_point->the_homologs->get_num_homologs(); i++) the_point->homolog_order[i]=new_order[i];
    delete[] new_order;
    
    
#ifdef _OPEN_MP_VERSION_
    pre_scores[mythread]=the_point->get_score();
}  //End parallel for
max=pre_scores[0];
max_loc=0;
for(mythread=1; mythread<omp_get_num_threads(); mythread++) {
    if (max > pre_scores[mythread]) {
        max=pre_scores[mythread];
        max_loc=mythread;
    }
}
(*current_walkers)[walker]->new_score=pre_scores[max_loc];
(*(*current_walkers)[walker]->point)=(*test_points[max_loc]);
//(*current_walkers)[walker]->point->set_blocks();
(*current_walkers)[walker]->point->set_blocks_v2();

#else
    //the_point->set_blocks();
    the_point->set_blocks_v2();
    (*current_walkers)[walker]->new_score= the_point->get_score();
#endif
    //(*current_walkers)[walker]->new_score= the_point->the_tracks->get_num_full_breaks();

	//cout<<"Walker: "<<walker<<endl;
	//cout<<" Old score: "<<(*current_walkers)[walker]->score<<"\t";
	//cout<<" New score: "<<(*current_walkers)[walker]->new_score<<"\n";

}  //End Groups_anneal::move







void Track_anneal::initialize_params(Space_point<Track_point> *curr_condition)
{
    //curr_condition->point->the_tracks = new WGX_Tracks(global_homologs, global_genomes);
    //curr_condition->point->the_tracks->update_tracking();

	cout<<"Initial score: "<<curr_condition->point->get_score()<<endl;
	
    //cout<<"Initial score: "<<curr_condition->point->the_tracks->get_num_full_breaks()<<endl;
    
	curr_condition->new_score=curr_condition->point->get_score();
    //curr_condition->new_score=curr_condition->point->the_tracks->get_num_full_breaks();
}


void Track_anneal::after_move(int walker_num)
{
	

}

int Track_anneal::choose_move_start(int walker)
{
    int retval;
    #pragma omp critical
    retval=ignuin(0, (*current_walkers)[walker]->point->get_num_blocks()-1);
   
    return(retval);
}


int Track_anneal::choose_far_loc(int walker, int block_start, int block_length)
{
    int new_block_loc;
    #pragma omp critical
    new_block_loc=ignuin(0,(*current_walkers)[walker]->point->get_num_blocks()-(1+block_length));
    
    if (new_block_loc >= block_start) new_block_loc+=block_length;
    
    return(new_block_loc);
}

void Track_anneal::output ()
{
  int i; 

 
  cout<<"Final Score "<<best->score<<endl;
  cout<<"No other output written yet"<<endl;
}

Track_anneal_density::Track_anneal_density() :Track_anneal()
{
}


Track_anneal_density::Track_anneal_density(Exchange *cexchange, int nwalkers) : Track_anneal(cexchange, nwalkers)
{
}

Track_anneal_density::~Track_anneal_density()
{
}



int Track_anneal_density::choose_move_start(int walker)
{
    int retval;
    #pragma omp critical
    retval=(*current_walkers)[walker]->point->location_density[ignuin(0, (*current_walkers)[walker]->point->density_size-1)];
    
    return(retval);
}


int Track_anneal_density::choose_far_loc(int walker, int block_start, int block_length)
{
    int new_block_loc;
    #pragma omp critical
    do {
        
        new_block_loc=(*current_walkers)[walker]->point->location_density[ignuin(0, (*current_walkers)[walker]->point->density_size-1)];
    } while ((new_block_loc >= block_start) && (new_block_loc < (block_start+block_length)));
    
    return(new_block_loc);
}
