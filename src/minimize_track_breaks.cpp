#define SEQUENCE_EXCHANGE
#include <iostream>
#include <fstream>
#include <math.h>
#include "exchange.h"
#include "gen_dna_funcs.h"
#include "random.h"
#include <string.h>
#include <map>
#include "libranlib.h"
#include "track_anneal_recombo.cpp"


using namespace::std;

Clade *global_genomes;
WGX_Data *global_homologs;

void parse_args(int argc, char** argv, Exchange *curr_exchange, int &num_genomes, string *&genome_files,
				string &ortho_file, int &d_level, int &max_breaks, string &outfile, double &boltz, double &relax, BOOL &use_basic_anneal, BOOL &guess_order, BOOL &use_pos_score, BOOL *&use_taxa);


int main(int argc, char *argv[])
{
    int i=0, j,k, ntaxa, nchars, avg_cnt, num_genomes, num_homologs, d_level, max_breaks;
    //long seed, seed2;
    double boltz, relax;
    BOOL use_basic_anneal, guess_order, use_pos_score, *use_taxa;
    string *genome_files, ortho_file, outfile;
    //ifstream seedin;
    //ofstream seedout;
    Random_Gen *mygen;
    
  //Classes
  Exchange current_exchange;
  Genome **list_of_genomes;
  Read_Genome *read_genome;
  
  Read_WGX_Data read_homologs;
  Track_anneal *curr_annealer;
 
  if (argc>3)
  {
      use_taxa=0;
      parse_args (argc, argv, &current_exchange, num_genomes, genome_files, ortho_file, d_level, max_breaks, outfile, boltz, relax, use_basic_anneal, guess_order, use_pos_score, use_taxa);
 
  
    //Gets the data
    list_of_genomes=new Genome * [num_genomes];
	for(i=0; i<num_genomes; i++) {
		read_genome=new Read_Genome();
		list_of_genomes[i]=read_genome->get_genome(genome_files[i]);
		delete read_genome;
	}
	
	global_genomes=new Clade(num_genomes, list_of_genomes);
	current_exchange.set_num_taxa(num_genomes);
      current_exchange.set_max_breaks(max_breaks);
      cout<<"Set max breaks to "<<max_breaks<<endl;
      
      if (guess_order==TRUE) {cout<<"Pre-optimizing initial order\n"; current_exchange.set_guess_order();}
      if (use_pos_score==TRUE) {
          current_exchange.set_use_break_pos();
          cout<<"Counting number of pillar positions with breaks for scoring\n";
      }
	//Clear the list we made, as the Clade object has its own copy
	for(i=0; i<num_genomes; i++)
		delete list_of_genomes[i];
	delete[] list_of_genomes;

	global_homologs=read_homologs.get_data(ortho_file, d_level, num_homologs, global_genomes);
      cout<<"Found "<<global_homologs->get_num_homologs()<<" homologous pillars in homolog file "<<ortho_file<<endl;


      //Make and log random number generator
      //seedin.open("./seed.txt");
      //seedin>>seed>>seed2;
      //seedin.close();
      //setall(seed, seed2);
      mygen=new Random_Gen;
        
      if (use_taxa !=0) {
          for(i=0; i<current_exchange.get_num_taxa(); i++) {
              if (use_taxa[i] == FALSE) {
                  current_exchange.mark_taxa_unused(i);
                  cout<<"Including genome number "<<i<<" in scoring function\n";
              }
          }
      }

      //cout<<"Begining global optimization\n"<<flush;
      if (use_basic_anneal==TRUE)
          curr_annealer = new Track_anneal(&current_exchange, 1);
      else
          curr_annealer = new Track_anneal_density(&current_exchange, 1);
      curr_annealer->set_relax(relax);
      curr_annealer->set_boltz(boltz);
      curr_annealer->set_temp_fact(2.0);
      curr_annealer->set_save_last(0);

      curr_annealer->begin_sim();
    
  
      //curr_annealer->best->score=curr_annealer->best->point->curr_model->find_appropriate_ln_like();
      cout<<"Final score: "<<curr_annealer->best->score<<endl;
	  curr_annealer->best->point->the_tracks->change_order(curr_annealer->best->point->homolog_order);
	  curr_annealer->best->point->the_tracks->update_tracking();
	

	  
      for(i=0; i<curr_annealer->best->point->the_homologs->get_num_homologs(); i++)
          cout<<curr_annealer->best->point->homolog_order[i]<<endl;
      
      curr_annealer->best->point->the_tracks->print_homolog_file(outfile);

	  //curr_annealer->best->point->the_tracks->print_all_tracks();

   
   
      //Output seed for next run
      //seedout.open("./seed.txt");
      //getsd(&seed, &seed2);
      //seedout<<seed<<"\t"<<seed2<<endl;
      //seedout.close();
      delete mygen;
      
	  delete global_homologs;
	  delete global_genomes;
      return(0);
  } 
  else
    {
        cerr<<"Usage: minimize_track_breaks -d:d_level -g:<genome file> -g:<genome file> -o:<ortholog file> -n:<new ortholog file> (-m:<Max Breaks>) (-b:<Boltz const>) (-r:<relax time>) (-infer_init_order) (-simpleanneal) (-usepositionscore) (-t:<TAXANAME>) (-t:<TAXANAME>)\n";
      return(-1);
    }
}  //End main




void parse_args(int argc, char** argv, Exchange *curr_exchange, int &num_genomes, string *&genome_files,
				string &ortho_file, int &d_level, int &max_breaks, string &outfile, double &boltz, double &relax, BOOL &use_basic_anneal, BOOL &guess_order, BOOL &use_pos_score, BOOL *&use_taxa)
{
	int i, j, cnt_genome=0, taxa_index;
    string dupl_level, relax_t, boltz_p, m_breaks, genome_name;
	LKMODEL the_model=DUPL_NOSTATE;
    ifstream genomein;
    map<string, int> taxa_hash;
    boltz=0.005;
    relax=200.0;

	num_genomes=0;
    max_breaks=0;
    use_basic_anneal=FALSE;
    guess_order=FALSE;
    use_pos_score=FALSE;
    
    cout<<"In parse: guess is "<<guess_order<<endl;

	for (i=0; i<argc; i++)
		if ((argv[i][1] == 'g') || (argv[i][i] == 'G')) num_genomes++;

	genome_files = new string [num_genomes];
	 
	for(i=1; i<argc; i++) {
		switch(argv[i][1]) {
		case 'd':
        case 'D':
                dupl_level=argv[i];
                dupl_level=dupl_level.substr(3, dupl_level.length()-3);
                d_level=string_to_int(dupl_level.c_str());
                break;
        case 'g':
		case 'G':
            genome_files[cnt_genome]=argv[i];
            genome_files[cnt_genome]=genome_files[cnt_genome].substr(3, genome_files[cnt_genome].length()-3);
            
			cnt_genome++;
		break;
	
		case 'n':
		case 'N':
            outfile=argv[i];
            outfile=outfile.substr(3, outfile.length()-3);
			
			break;
		case 'o':
		case 'O':
            ortho_file=argv[i];
            ortho_file=ortho_file.substr(3, ortho_file.length()-3);
        break;
        case 'r':
        case 'R':
                relax_t=argv[i];
                relax_t=relax_t.substr(3, relax_t.length()-3);
                relax=string_to_float(relax_t.c_str());
        break;
        case 'm':
        case 'M':
                m_breaks=argv[i];
                m_breaks=m_breaks.substr(3, m_breaks.length()-3);
                max_breaks=string_to_int(m_breaks.c_str());
                break;
        case 's':
        case 'S':
                use_basic_anneal=TRUE;
                break;
        case 'i':
        case 'I':
                cout<<"Found i: resetting guess to "<<guess_order<<endl;
                guess_order=TRUE;
                break;
        case 'u':
        case 'U':
                use_pos_score=TRUE;
                break;
        case 'b':
        case 'B':
                boltz_p=argv[i];
                boltz_p=boltz_p.substr(3, boltz_p.length()-3);
                boltz=string_to_float(boltz_p.c_str());
        break;
        case 't':
        case 'T':
                if (use_taxa ==0) {
                    use_taxa=new BOOL [num_genomes];
                    for(taxa_index=0; taxa_index<num_genomes; taxa_index++) {
                        use_taxa[taxa_index] =FALSE;
                        genomein.open(genome_files[taxa_index].c_str());
                        genomein>>genome_name;
                        genomein.close();
                        taxa_hash[genome_name]=taxa_index;
                        cout<<"Genome "<<genome_name<<" found in "<<genome_files[taxa_index]<<" and will be number "<<taxa_index<<endl;
                    }
                    
                }
                
                genome_name=argv[i];
                genome_name=genome_name.substr(3,genome_name.length()-3);
                use_taxa[taxa_hash[genome_name]]=TRUE;
                cout<<"Genome "<<genome_name<<" is number "<<taxa_hash[genome_name]<<" and will be scored\n";
                break;
        }

	}
	curr_exchange->set_model(the_model);
  
}
