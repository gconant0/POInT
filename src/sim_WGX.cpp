#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "exchange.h"
#include "sim_data.h"
#include "read_seq.h"
#include "maxlike.h"
#include "codon_like.h"
#include "tree.h"
#include "random.h"
#include "read_tree.h"
#include "write_tree.h"
#include "genome_list.h"
#include "genome_ploidy_like.h"
#include "genome_tripl_list.h"
#include "phylo_model_matrix.h"
#include "genome_ploidy_like.h"
#include <sstream>
#include <map>

#ifndef POInT_version
#define POInT_version "v1.6"
#endif

void parse_args(int argc, char** argv, Exchange *curr_exchange, int &num_genomes, string *&genome_files,
				string *&output_genome_files, string &ortho_file, string &outorder, string &model_file, string &outtree);


int main (int argc, char **argv)
{
    int i, j, k, l, nchars, num_genomes, num_homologs, **num_contigs, **contig_num, dupl_pos,
	  last_start, num_genes, **gene_num, contig_cnt, *walk_masks, sum, *curr_ortho_pattern, *new_ortho_pattern;
    //long seed, seed2;
    double select_track_prob, cum_prob;
    //ifstream seedin, tree_scan;
    //ofstream seedout, genomeout, homologsout;
    ifstream tree_scan;
    ofstream genomeout, homologsout;
    string  outorder, *output_genome_files, *genome_files, ortho_file, model_file, ***gene_names, *contig_names, outtree, name_string, prog_name;
    char line[400], temp[20],  genome_name[50];
    BOOL has_valid_genes;
    Simulate_data *simulator;
    Ploidy_Like_model *the_model;
    Sequence_dataset *dummy_sequences, *writable_dataset;              //Required for use of Like_model classes
    Read_PAUP_Tree_Ex get_tree;
    Write_Tree_Arb_Model writeout_tree;
    Tree *current_tree, *new_tree;
    Exchange current_exchange;
    Contig ****new_contigs;
    Genome **list_of_genomes;
    Read_Genome *read_genome;
    Clade *the_genomes;
    WGX_Data *the_homologs;
    Read_WGX_Data read_homologs;
    Phylo_Matrix *curr_model_matrix;
    Random_Gen *mygen;
    map<string, string> name_map;
 
    if(argc<5) {
        cerr<<"Usage: sim_WGX -d# -g:<Genome1> -g:<Genome2>... -o:<ortholog_file>  -t:<Treefile> -m:<ModelFile> -x:<Output_genomes> -w:<Output_orthologs> -z:<OutputtreeFile>\n";
        return(-1);
    }
    else {
        parse_args(argc, argv, &current_exchange, num_genomes, genome_files, output_genome_files,
		  ortho_file, outorder, model_file, outtree);
        current_exchange.set_have_data(FALSE);
        current_exchange.set_num_rates(1);
      

        //Make and log random number generator
        //seedin.open("./seed.txt");
        //seedin>>seed>>seed2;
        //seedin.close();
        //setall(seed, seed2);
        mygen=new Random_Gen;
     
        current_exchange.set_num_taxa(num_genomes);
        list_of_genomes=new Genome * [num_genomes];
    
        for(i=0; i<num_genomes; i++) {
			read_genome=new Read_Genome();
			list_of_genomes[i]=read_genome->get_genome(genome_files[i]);
			delete read_genome;
		}
		the_genomes=new Clade(num_genomes, list_of_genomes);
		current_exchange.set_num_taxa(num_genomes);
		
        curr_model_matrix=new Phylo_Matrix(&current_exchange, model_file);
        
		//Clear the list we made, as the Clade object has its own copy
		for(i=0; i<num_genomes; i++)
			delete list_of_genomes[i];
		delete[] list_of_genomes;

		the_homologs=read_homologs.get_data(ortho_file, curr_model_matrix->get_num_levels(), num_homologs, the_genomes);


        dummy_sequences=new Sequence_dataset(the_genomes->get_num_genomes(), 1, ARBITRARY);


        for(i=0; i<the_genomes->get_num_genomes(); i++) {
			(*dummy_sequences)[i].Assign_name((*the_genomes)[i].get_name());
            cout<<"Assigning "<<(*the_genomes)[i].get_name()<<" to seq dummy "<<i<<endl;
        }

		
		current_exchange.set_dataset(dummy_sequences);
		current_exchange.set_num_taxa(the_genomes->get_num_genomes());
	    current_exchange.set_num_sites(the_homologs->get_num_homologs());
        current_exchange.set_model(DUPL_ARBITRARY, curr_model_matrix->get_num_states());
        
		current_tree=get_tree.create_tree_from_file(&current_exchange, current_exchange.get_dataset(), curr_model_matrix, TRUE);

	for(i=0; i<current_exchange.get_num_branches(); i++) {
		cout<<i<<" len: "<<(*current_tree)[i]->get_brnlen()<<endl;
        (*current_tree)[i]->set_brnlen((*current_tree)[i]->expect_subs_site());
		cout<<i<<" len: "<<(*current_tree)[i]->get_brnlen()<<endl;
	}
        the_model=new Ploidy_Like_model(&current_exchange, current_tree, the_genomes, the_homologs, curr_model_matrix);
//	the_model->recalculate_transprobs();
        the_model->the_tracks->update_tracking();
        
		gene_names=new string **[num_genomes];
		for(i=0; i<num_genomes; i++) {
			gene_names[i]=new string *[the_homologs->get_dupl_level()];
			for(j=0; j<the_homologs->get_dupl_level(); j++) {
				gene_names[i][j]=new string [num_homologs];
				for(k=0; k<num_homologs; k++) {
                    ostringstream oss;
                        oss << "SP_" << i <<"_g" <<k <<char(65+j);
                    
                     gene_names[i][j][k] = oss.str ();
				}
			}
		}


		num_contigs=new int*[num_genomes];
		contig_num=new int*[num_genomes];
		new_contigs=new Contig***[num_genomes];
		
		for(i=0; i<num_genomes; i++) {
			num_contigs[i]=new int [the_homologs->get_dupl_level()];
			contig_num[i]=new int[the_homologs->get_dupl_level()];
            for(j=0; j<the_homologs->get_dupl_level(); j++)
                num_contigs[i][j]=1;
            
			for(j=0; j<num_homologs-1; j++) {
                for(dupl_pos=0; dupl_pos<the_homologs->get_dupl_level(); dupl_pos++) {
                    if (((*the_model->the_tracks)[i].has_back_link(j+1,dupl_pos) == FALSE) &&
                        ((*the_model->the_tracks)[i].get_gene_track(j,dupl_pos)->my_locus !=0) )
                    {
                        num_contigs[i][dupl_pos]++;
                    }
                    else {
                        if (((*the_model->the_tracks)[i].has_back_link(j, dupl_pos) == FALSE) &&
                            ((*the_model->the_tracks)[i].get_gene_track(j, dupl_pos)->my_locus != 0) && (j != 0) &&
                            ((*the_model->the_tracks)[i].get_gene_track(j-1, dupl_pos)->my_locus ==0)) {
                            num_contigs[i][dupl_pos]++;
                        }
                    }
                }
			}
		}

		for(i=0; i<num_genomes; i++) {
			new_contigs[i]=new Contig**[the_homologs->get_dupl_level()];
			for(j=0; j<the_homologs->get_dupl_level(); j++) {
                cout<<"Creating "<<num_contigs[i][j]<<" Contigs for genome "<<i<<" strand "<<j<<endl;
				new_contigs[i][j]=new Contig*[num_contigs[i][j]];
				contig_num[i][j]=0;
		
				k=0;
				last_start=0;
				sum=0;
		
				while(k<num_homologs) {
					if (k != num_homologs-1) {
						if (((*the_model->the_tracks)[i].get_gene_track(k,j)->my_locus != 0) &&
							((*the_model->the_tracks)[i].has_back_link(k+1, j) == FALSE)) {
							num_genes=k-last_start+1;
							contig_names=new string [num_genes];
						
							for (l=0; l<num_genes; l++)	contig_names[l]=gene_names[i][j][last_start+l];

							
							new_contigs[i][j][contig_num[i][j]]=new Contig(num_genes, contig_names);
							
							delete[] contig_names;
							contig_num[i][j]++;
							last_start=k+1;
							sum+=num_genes;
							if (sum < k )
								cerr<<"ERROR: missed genes\n";
                        }
						else if (((*the_model->the_tracks)[i].get_gene_track(k,j)->my_locus != 0) &&
							((*the_model->the_tracks)[i].has_back_link(k, j) == FALSE)  && (k != 0) &&
							((*the_model->the_tracks)[i].get_gene_track(k-1, j)->my_locus ==0)) {
							num_genes=k-1-last_start+1;
							contig_names=new  string [num_genes];
							for (l=0; l<num_genes; l++) contig_names[l]= gene_names[i][j][last_start+l];

							new_contigs[i][j][contig_num[i][j]]=new Contig(num_genes, contig_names);
							
							delete[] contig_names;
							contig_num[i][j]++;
							last_start=k;
							sum+=num_genes;
							if (sum < k)
								cerr<<"ERROR: missed genes\n";
						}

					}
					else {
					
						num_genes=k-last_start+1;
						if (num_genes != 0) {
							contig_names=new string [num_genes];
							for (l=0; l<num_genes; l++) contig_names[l]= gene_names[i][j][last_start+l];
                        
							new_contigs[i][j][contig_num[i][j]]=new Contig(num_genes, contig_names);
							
							delete[] contig_names;
							contig_num[i][j]++;
							sum+=num_genes;
						}
						
					}
					k++;
				}
			}

		}


		for(i=0; i<num_genomes; i++) {
			for(j=0; j<the_homologs->get_dupl_level(); j++) {
				last_start=0;
				for(k=0; k<num_contigs[i][j]; k++)
					last_start+=new_contigs[i][j][k]->get_num_genes();
				cout<<"Total genes: "<<i<<"\t"<<j<<"\t"<<last_start<<endl;
			}
		}

		cout<<"Contigs: "<<num_contigs[0][0]<<" counted: "<<contig_num[0][0]<<endl;
	
        for(i=0; i<num_genomes; i++) {
            stringstream ss;
            ss << "Genome" << i;
            name_string=ss.str();
            name_map[(*the_genomes)[i].get_name()]=name_string;
        }
        new_tree=new Tree(&current_exchange, TRUE);
        (*new_tree)=(*current_tree);
    
        for(i=0; i<current_exchange.get_num_branches(); i++) {
            if ((*new_tree)[i]->is_tip() == TRUE) {
                name_string=name_map[(*new_tree)[i]->get_name()];
                (*new_tree)[i]->set_name(name_string.c_str());
            }
        }
     
        
        
        cout<<"Simulating "<<num_genomes<<" taxa and "<<num_homologs<<" sites\n";
for(i=0; i<current_exchange.get_num_branches(); i++) {
                cout<<i<<" len: "<<(*current_tree)[i]->get_brnlen()<<endl;     
}
	the_model->recalculate_transprobs();
        simulator=new Simulate_from_PhyloMatrix(current_tree, &current_exchange, the_model, curr_model_matrix);
for(i=0; i<current_exchange.get_num_branches(); i++) {
                cout<<i<<" len: "<<(*current_tree)[i]->get_brnlen()<<endl;   
}
        simulator->make_dataset();
	for(i=0; i<current_exchange.get_num_branches(); i++) {
                cout<<i<<" len: "<<(*current_tree)[i]->get_brnlen()<<endl;
}
        curr_ortho_pattern= new int [num_genomes];
        new_ortho_pattern=new int[num_genomes];
		

	  gene_num=new int*[num_genomes];
	  for(i=0; i<num_genomes; i++) {
		   gene_num[i]=new int[the_homologs->get_dupl_level()];
		  
          for(dupl_pos=0; dupl_pos< the_homologs->get_dupl_level(); dupl_pos++) {
              contig_num[i][dupl_pos]=0;
              gene_num[i][dupl_pos]=0;
		}
          
          curr_ortho_pattern[i]=0;
          new_ortho_pattern[i]=0;
	  }
	
        for(i=0; i<num_homologs; i++) {
            the_model->set_site_base_states(i);
            if (i >0)
                the_model->build_transition_vector(i);
            
            
            for(j=0; j<num_genomes; j++) {
                if (i>0){
                    select_track_prob=ranf();
                    new_ortho_pattern[j]=0;
                    cum_prob=0;
                    
                    cout<<"taxa "<<j<<": Site: "<<i<<" State: "<<curr_model_matrix->get_nth_state((*simulator->dataset)[j][i])->get_state_id()<<" called "<<curr_model_matrix->get_nth_state((*simulator->dataset)[j][i])->get_state_name()<<": "<<curr_ortho_pattern[j]<<" Drew: "<<select_track_prob<<" ...";
                    
                    while ((cum_prob < select_track_prob)  && (new_ortho_pattern[j]<(the_model->get_num_track_patterns()-1))) {
                        cum_prob += the_model->get_track_transprob(j, curr_ortho_pattern[j], new_ortho_pattern[j]);
                        cout<<"\t"<<new_ortho_pattern[j]<<": "<<cum_prob<<" From: "<<the_model->get_track_transprob(j, curr_ortho_pattern[j], new_ortho_pattern[j]);
                        if(cum_prob < select_track_prob)
                            new_ortho_pattern[j]++;
                    }
                   cout<<" Final track state: "<<new_ortho_pattern[j]<<endl;
                }
                
                
                for(dupl_pos=0; dupl_pos< the_homologs->get_dupl_level(); dupl_pos++) {
                    //CEHCK!!!!
			cout<<"\tFor this track, dupl level "<<dupl_pos<<" is "<<curr_model_matrix->get_nth_state((*simulator->dataset)[j][i])->position_present(the_model->get_track_position_state(new_ortho_pattern[j], dupl_pos))<<endl;
                    if (curr_model_matrix->get_nth_state((*simulator->dataset)[j][i])->position_present(the_model->get_track_position_state(new_ortho_pattern[j], dupl_pos)) ==FALSE)
                   { //if (curr_model_matrix->get_nth_level_ith_state(curr_model_matrix->get_nth_state((*simulator->dataset)[j][i])->get_state_level(),
                    //                                               the_model->get_track_position_state(new_ortho_pattern[j], dupl_pos)) == 0) {
                        (*new_contigs[j][dupl_pos][contig_num[j][dupl_pos]])[gene_num[j][dupl_pos]].discard_gene();
                	cout<<"Setting contig "<<contig_num[j][dupl_pos]<<" gene "<<gene_num[j][dupl_pos]<<" with name "<<(*new_contigs[j][dupl_pos][contig_num[j][dupl_pos]])[gene_num[j][dupl_pos]].get_name_string()
				<<" keep to "<<(*new_contigs[j][dupl_pos][contig_num[j][dupl_pos]])[gene_num[j][dupl_pos]].keep_gene()<<endl;    
		}
                    
                    if (gene_num[j][dupl_pos] == new_contigs[j][dupl_pos][contig_num[j][dupl_pos]]->get_num_genes()-1) {
                        gene_num[j][dupl_pos]=0;
                        contig_num[j][dupl_pos]++;
                    }
                    else gene_num[j][dupl_pos]++;
                }
                
                       
            }
            for(j=0; j<num_genomes; j++) curr_ortho_pattern[j]=new_ortho_pattern[j];
          // cout<<i<<": "<<curr_ortho_pattern[0]<<endl;
        }
	 
	
	 

	  for(i=0; i<num_genomes; i++) {
		contig_cnt=1;
		genomeout.open(output_genome_files[i].c_str());
		genomeout<<"Genome"<<i<<"\n";
		for(j=0; j<the_homologs->get_dupl_level(); j++) {
			for(k=0; k<num_contigs[i][j]; k++) 
			{	
				has_valid_genes=FALSE;
				for(l=0; l<new_contigs[i][j][k]->get_num_genes(); l++) {
					if ((*new_contigs[i][j][k])[l].keep_gene() == TRUE)
						has_valid_genes=TRUE;
				}

				if (has_valid_genes == TRUE) {
					for(l=0; l<new_contigs[i][j][k]->get_num_genes(); l++) {
						if ((*new_contigs[i][j][k])[l].keep_gene() == TRUE)
							genomeout<<contig_cnt<<"\t"<<(*new_contigs[i][j][k])[l].get_name()<<endl;
					}
					contig_cnt++;
				}
			}
		}
		genomeout.close();
		genomeout.clear();
	  }

	  homologsout.open(outorder.c_str());

	  for(i=0; i<num_genomes; i++) {
		homologsout<<"Genome"<<i;
          for(dupl_pos=0; dupl_pos<the_homologs->get_dupl_level(); dupl_pos++) {
              contig_num[i][dupl_pos]=0;
              gene_num[i][dupl_pos]=0;
          }
		if (i != num_genomes -1)
			homologsout<<"\t";
		else
			homologsout<<"\n";
	  }

	  for(i=0; i<num_homologs; i++) {
          for(dupl_pos=0; dupl_pos<the_homologs->get_dupl_level(); dupl_pos++) {
                for(j=0; j<num_genomes; j++) {
//                      cout<<"Pos "<<i<<" taxa "<<j<<" DP: "<<dupl_pos<<" Keep is "<<(*new_contigs[j][dupl_pos][contig_num[j][dupl_pos]])[gene_num[j][dupl_pos]].keep_gene()<<endl;
			if ((*new_contigs[j][dupl_pos][contig_num[j][dupl_pos]])[gene_num[j][dupl_pos]].keep_gene() == TRUE)
                          homologsout<<(*new_contigs[j][dupl_pos][contig_num[j][dupl_pos]])[gene_num[j][dupl_pos]].get_name();
                      else
                          homologsout<<"NONE";
                      
                      if ((j == num_genomes-1) &&( dupl_pos == the_homologs->get_dupl_level()-1))
                        homologsout<<"\n";
                      else
                        homologsout<<"\t";
                  }
          }
			  

          for(k=0; k<the_homologs->get_dupl_level(); k++) {
              for(j=0; j<num_genomes; j++) {
                 if (gene_num[j][k] == new_contigs[j][k][contig_num[j][k]]->get_num_genes()-1) {
                    gene_num[j][k]=0;
                    contig_num[j][k]++;
                  }
                  else gene_num[j][k]++;
              }
          }
			  
			  
	  }

	  homologsout.close();
        stringstream ss2;
        ss2 << "Tree topology for POInT_simulate "<< POInT_version;
        prog_name =ss2.str();
       
        writeout_tree.write_tree(outtree, prog_name, curr_model_matrix, 0, new_tree, &current_exchange);
        
      //Output seed for next run
      //seedout.open("./seed.txt");
      //getsd(&seed, &seed2);
      //seedout<<seed<<"\t"<<seed2<<endl;
      //seedout.close();
        delete mygen;
		
	  for(i=0; i<num_genomes; i++) {
		  for (j=0; j<the_homologs->get_dupl_level(); j++) {
				  delete[] gene_names[i][j];
		  }
          delete[] gene_names[i];
	  }
        delete[] gene_names;

	  delete[] genome_files;
	  delete[] output_genome_files;
      delete simulator;
      //delete model;
      delete current_tree;
      return(0);

    }
}//end main


void parse_args(int argc, char** argv, Exchange *curr_exchange, int &num_genomes,  string *&genome_files,
				string *&output_genome_files,  string &ortho_file, string &outorder, string &model_file, string &outtree)
{
	int i, j, cnt_genome=0;
	string treefile,  outgenome_stub;

	num_genomes=0;

	for (i=0; i<argc; i++) {
		if ((argv[i][1] == 'g') || (argv[i][1] == 'G')) num_genomes++;
		if ((argv[i][1] == 'x') || (argv[i][1] == 'X')) {
            outgenome_stub=argv[i];
            outgenome_stub=outgenome_stub.substr(3,outgenome_stub.length()-3);
		}
	}

	genome_files = new string [num_genomes];
	output_genome_files =new string [num_genomes];
	for(i=0; i<num_genomes; i++) {
        ostringstream oss;
        
        oss << outgenome_stub << "_" << i <<".txt";
        output_genome_files[i] = oss.str ();
	}



	for(i=0; i<argc; i++) {
		switch(argv[i][1]) {
            case 'g':
            case 'G':
            genome_files[cnt_genome]=argv[i];
            genome_files[cnt_genome]=genome_files[cnt_genome].substr(3, genome_files[cnt_genome].length()-3);
            cnt_genome++;
            break;
		case 'w':
		case 'W':
            outorder=argv[i];
            outorder=outorder.substr(3,outorder.length()-3);
			break;
        case 'z':
        case 'Z':
                outtree=argv[i];
                outtree=outtree.substr(3,outtree.length()-3);
                break;
        case 't':
        case 'T':
            treefile=argv[i];
            treefile=treefile.substr(3, treefile.length()-3);
            curr_exchange->set_treefile(treefile.c_str());
			break;
            
        case 'o':
        case 'O':
            ortho_file=argv[i];
            ortho_file=ortho_file.substr(3, ortho_file.length()-3);
            break;
            
        case 'm':
        case 'M':
            model_file=argv[i];
            model_file=model_file.substr(3, model_file.length()-3);
            break;
		}

	}

}








