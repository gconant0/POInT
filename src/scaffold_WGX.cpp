#define SEQUENCE_EXCHANGE
#include <iostream>
#include <fstream>
#include <math.h>
#include "exchange.h"
#include "gen_dna_funcs.h"
#include "random.h"
#include <string.h>
#include "libranlib.h"
#include "WGX_scaffold_anneal.cpp"


using namespace::std;

Homolog_Set *global_homologs;


void parse_args(int argc, char** argv, Exchange *curr_exchange, int &WGD_depth, string &ances_genome_file,
				string &dupl_genome_file, string &ances_tandem_file, string &dupl_tandem_file, string &homolog_file, double &homo_cutoff, double &relax_time, double &boltz, int &count_cutoff, string &outfile, string &order_file);


int main(int argc, char *argv[])
{
    int i=0, j,k, ntaxa, nchars, avg_cnt, num_genomes, num_homologs, d_level, count_cutoff;
    //long seed, seed2;
    double homo_cutoff, relax_time, boltz;
    string ances_genome_file, dupl_genome_file, ances_tandem_file, dupl_tandem_file, homolog_file, outfile, order_file, tandem_file;
    //ifstream seedin;
    //ofstream seedout;
    Random_Gen *mygen;
 
    //Classes
    Exchange current_exchange;
    Scaffold_Genome *orig_ances_genome, *orig_dupl_genome, *ances_genome, *dupl_genome;
    Scaffold_anneal *curr_annealer;
 
    if (argc>3) {
        parse_args(argc, argv, &current_exchange, d_level, ances_genome_file,
                   dupl_genome_file, ances_tandem_file, dupl_tandem_file, homolog_file, homo_cutoff, relax_time, boltz, count_cutoff, outfile, order_file);
 
        //cout<<"Parsed args\n"<<flush;
        current_exchange.set_WGX_depth(d_level);
  
        orig_ances_genome=read_genome(ances_genome_file);
        orig_dupl_genome=read_genome(dupl_genome_file);
        
        cout<<"Initial ancestral genome has "<<orig_ances_genome->get_num_genes()<<" genes\n";
        cout<<"Initial duplicated genome has "<<orig_dupl_genome->get_num_genes()<<" genes\n";
        
        ances_genome=collapse_tandems (ances_tandem_file, orig_ances_genome);
        dupl_genome=collapse_tandems(dupl_tandem_file, orig_dupl_genome);
        
        cout<<"After tandem collapse, ancestral genome has "<<ances_genome->get_num_genes()<<" genes\n";
        cout<<"After tandem collapse, duplicated genome has "<<dupl_genome->get_num_genes()<<" genes\n";
        
        
        tandem_file=ances_genome->get_genome_name() + "_" + dupl_genome->get_genome_name() + "_tandem_concat.txt";
        ances_genome->print_tandems(tandem_file);
        
        tandem_file=dupl_genome->get_genome_name() + "_tandem_concat.txt";
        dupl_genome->print_tandems(tandem_file);
        
        
        global_homologs=read_homolog_set (homolog_file, ances_genome, dupl_genome, homo_cutoff, count_cutoff);
        
        cout<<"After homolog pruning, there are "<<global_homologs->get_ancestral_genome()->get_num_genes()<<" genes in the ancestral genome\n";
        cout<<"After homolog pruning, there are "<<global_homologs->get_dupl_genome()->get_num_genes()<<" genes in the duplicated genome\n";
        cout<<"Order file is "<<order_file<<endl;
        if (order_file != "NONE") global_homologs->get_new_order(order_file);
        
        current_exchange.set_num_taxa(2);
		
	
        //Make and log random number generator
        //seedin.open("/Users/gavinconant/software/gavins/scaffold_WGX_proj/seed.txt");
        //seedin>>seed>>seed2;
        //seedin.close();
        //setall(seed, seed2);
        mygen=new Random_Gen;
        
 

        //cout<<"Begining global optimization\n"<<flush;
        curr_annealer = new Scaffold_anneal(&current_exchange, 1);
        curr_annealer->set_outfile(outfile);
        curr_annealer->set_relax(relax_time);
        curr_annealer->set_boltz(boltz);
        curr_annealer->set_temp_fact(2.0);
        curr_annealer->set_save_last(0);

        curr_annealer->begin_sim();
    
  
        //curr_annealer->best->score=curr_annealer->best->point->curr_model->find_appropriate_ln_like();
        cout<<"Final score: "<<curr_annealer->best->score<<endl;
	  
	   
   
        //Output seed for next run
        delete mygen;
        //seedout.open("seed.txt");
        //getsd(&seed, &seed2);
        //seedout<<seed<<"\t"<<seed2<<endl;
        //seedout.close();
      
        delete global_homologs;
        return(0);
  } 
  else
    {
      cerr<<"Usage: scaffold_WGD -d:WGD_depth -g:<Dupl. genome file> -a:<Ances. genome file> -t:<Ances. tandem file> -u:<Dupl. tandem file> -h:<Homolog file> -c:Homolog_Cut -r:RelaxTime -b:Boltz_const -o:<output file> (-i:initial order file)\n";
      return(-1);
    }
}  //End main




void parse_args(int argc, char** argv, Exchange *curr_exchange, int &WGD_depth, string &ances_genome_file,
                string &dupl_genome_file, string &ances_tandem_file, string &dupl_tandem_file, string &homolog_file, double &homo_cutoff, double &relax_time, double &boltz, int &count_cutoff, string &outfile, string &order_file)
{
	int i, j;
    string depth, cutoff, relax_t, boltz_s;
	
    homo_cutoff=0.5;
    WGD_depth=2;
    relax_time=200.0;
    boltz=0.002;
    order_file="NONE";
    count_cutoff=-1;
    
	for(i=1; i<argc; i++) {
        cout<<"Parsing: "<<argv[i]<<endl<<flush;
		switch(argv[i][1]) {
            case 'd':
            case 'D':
                depth=argv[i];
                depth=depth.substr(3, depth.length()-3);
                WGD_depth=string_to_int(depth.c_str());
                break;
            case 'r':
            case 'R':
                relax_t=argv[i];
                relax_t=relax_t.substr(3, relax_t.length()-3);
                relax_time=string_to_float(relax_t.c_str());
                break;
            case 'g':
            case 'G':
                dupl_genome_file=argv[i];
                dupl_genome_file=dupl_genome_file.substr(3, dupl_genome_file.length()-3);
                break;
            case 'a':
            case 'A':
                ances_genome_file=argv[i];
                ances_genome_file=ances_genome_file.substr(3, ances_genome_file.length()-3);
                break;
            case 't':
            case 'T':
                ances_tandem_file=argv[i];
                ances_tandem_file=ances_tandem_file.substr(3, ances_tandem_file.length()-3);
                break;
            case 'u':
            case 'U':
                dupl_tandem_file=argv[i];
                dupl_tandem_file=dupl_tandem_file.substr(3, dupl_tandem_file.length()-3);
                break;
            case 'h':
            case 'H':
                homolog_file=argv[i];
                homolog_file=homolog_file.substr(3, homolog_file.length()-3);
                break;
            case 'o':
            case 'O':
                outfile=argv[i];
                outfile=outfile.substr(3, outfile.length()-3);
                break;
            case 'c':
            case 'C':
                cutoff=argv[i];
                cutoff=cutoff.substr(3, cutoff.length()-3);
                homo_cutoff=string_to_float(cutoff.c_str());
                break;
            case 'b':
            case 'B':
                boltz_s=argv[i];
                boltz_s=boltz_s.substr(3, boltz_s.length()-3);
                boltz=string_to_float(boltz_s.c_str());
                break;
            case 'i':
            case 'I':
                order_file=argv[i];
                order_file=order_file.substr(3, order_file.length()-3);
                break;
            case 'n':
            case 'N':
                cutoff=argv[i];
                cutoff=cutoff.substr(3, cutoff.length()-3);
                count_cutoff=string_to_int(cutoff.c_str());
                break;
		}
	}
  
}



