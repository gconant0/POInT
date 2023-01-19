#include "sim_data.h"
#include <iostream>
#include <iomanip>

using namespace::std;

Simulate_data::Simulate_data()
{
  cerr<<"Default Constructor called for base class Simulate_data\n";
}


Simulate_data::Simulate_data(Tree *ctree, Exchange *cexchange, Like_model *cmodel)
{
  int i;

  curr_tree=ctree;
  curr_exchange=cexchange;
  curr_model=cmodel;

  dataset =new Sequence_dataset(curr_exchange->get_num_taxa(),curr_exchange->get_num_sites()) ;

  for(i=0; i<curr_exchange->get_num_taxa(); i++)
    (*dataset)[i].Assign_name((*curr_tree)[i]->get_name());
 
}

void Simulate_data::make_dataset()
{
  int i;

  for(i=0; i<curr_exchange->get_num_localities(); i++)
      get_next_site(i, (int)ignuin(0, curr_exchange->get_num_rates()-1));
    
}


void Simulate_data::make_dataset(double **rate_probs, int num_probs, int *used_rates)
{
	int i, rate, choose_prob;
	double so_far, val;
	
	for(i=0; i<curr_exchange->get_num_localities(); i++) {
		choose_prob=ignuin(0, num_probs-1);
		val=ranf();
		used_rates[i]=choose_prob;
		
		rate=0;
		so_far = rate_probs[choose_prob][rate];
		while (so_far <= val) {
			rate++;
			so_far += rate_probs[choose_prob][rate];
			
		}
		//cout<<i<<"\t"<<rate<<"\t"<<val<<"\tSF:"<<so_far<<endl;
		get_next_site(i, rate);
	}
}


//Simulate_data protected functions

void Simulate_data::get_next_site(int site_num, int rate)
{
  ascend_tree(get_root_site(), site_num, curr_tree->find_root(), rate);
}  //End Simulate_data::get_next_site




void Simulate_data::ascend_tree(int curr_site, int site_num, Branch *curr_branch, int rate)
{
  int new_site=0;
  double p_so_far, ran_p;


  if (curr_branch->get_brnlen()>FLOAT_TOL)
    {
      ran_p=ranf();
      
      p_so_far=curr_branch->get_trpb(rate, curr_site, new_site);
      
      while(p_so_far<ran_p && new_site<curr_exchange->get_condlike_size()-1)
	{
	  new_site++;
	  p_so_far+=curr_branch->get_trpb(rate, curr_site, new_site);
	}
     
    }
  else
    new_site=curr_site;

  if (curr_branch->is_tip()==TRUE)
      write_data(new_site, site_num, curr_branch);
 
  else
    {
      ascend_tree(new_site, site_num, curr_branch->get_child(0), rate);
      ascend_tree(new_site, site_num, curr_branch->get_child(1), rate);
    }
   
}  //End Simulate_data::ascend_tree



Simulate_Nucleotide::Simulate_Nucleotide() : Simulate_data()
{
}


Simulate_Nucleotide::Simulate_Nucleotide(Tree *ctree, Exchange *cexchange, Like_model *cmodel) : Simulate_data
(ctree, cexchange, cmodel)
{
}




void Simulate_Nucleotide::write_data(int site, int site_num, Branch *taxa)
{
  (*dataset)[taxa->get_taxa_id()].Assign_site(site_num, site);  
}  //End Simulate_Codon::write_data




int Simulate_Nucleotide::get_root_site()
{
  int i,  nuc_num;
  double p_so_far, ran_p;
 
  ran_p=ranf();
  
  nuc_num=0;
  p_so_far=curr_model->root_freq(0);
  
  while(p_so_far<ran_p && nuc_num<curr_exchange->get_condlike_size()-1)
    {
      nuc_num++;
      p_so_far+=curr_model->root_freq(nuc_num);
    }

  return(nuc_num);
}







Simulate_Codon::Simulate_Codon() : Simulate_data()
{
}


Simulate_Codon::Simulate_Codon(Tree *ctree, Exchange *cexchange, Like_model *cmodel, Genetic_code *ccode) : Simulate_data
(ctree, cexchange, cmodel)
{
  curr_code=ccode;
}




void Simulate_Codon::write_data(int site, int site_num, Branch *taxa)
{
  (*dataset)[taxa->get_taxa_id()].Assign_site(3*site_num, curr_code->get_pos_n(site, 0));
  (*dataset)[taxa->get_taxa_id()].Assign_site(3*site_num+1, curr_code->get_pos_n(site, 1));
  (*dataset)[taxa->get_taxa_id()].Assign_site(3*site_num+2, curr_code->get_pos_n(site, 2));
}  //End Simulate_Codon::write_data




int Simulate_Codon::get_root_site()
{
  int i,  codon_num;
  double p_so_far, ran_p;
 
  ran_p=ranf();
  
  codon_num=0;
  p_so_far=curr_model->root_freq(0);
  
  while(p_so_far<ran_p && codon_num<curr_exchange->get_condlike_size()-1)
    {
      codon_num++;
      p_so_far+=curr_model->root_freq(codon_num);
    }

  while (curr_code->is_stop(codon_num)==TRUE )
    {
      ran_p=ranf();
      
      codon_num=0;
      p_so_far=curr_model->root_freq(0);
      
      while(p_so_far<ran_p && codon_num<curr_exchange->get_condlike_size()-1)
	{
	  codon_num++;
	  p_so_far+=curr_model->root_freq(codon_num);
	}
    }

  return(codon_num);
}




Simulate_Dupl::Simulate_Dupl() : Simulate_data()
{
}


Simulate_Dupl::Simulate_Dupl(Tree *ctree, Exchange *cexchange, Like_model *cmodel) : Simulate_data(ctree, cexchange, cmodel)
{

}
 

void Simulate_Dupl::write_data(int site, int site_num, Branch *taxa)
{
	 (*dataset)[taxa->get_taxa_id()].Assign_site(site_num, site); 
}


int Simulate_Dupl::get_root_site()
{
	return(loss_state_to_dupl(BOTH_PRESENT));

}




Simulate_from_PhyloMatrix ::Simulate_from_PhyloMatrix () : Simulate_data()
{
}


Simulate_from_PhyloMatrix ::Simulate_from_PhyloMatrix (Tree *ctree, Exchange *cexchange, Like_model *cmodel, Phylo_Matrix *cmatrix) : Simulate_data(ctree, cexchange, cmodel)
{
    curr_model_matrix=cmatrix;
}


void Simulate_from_PhyloMatrix::write_data(int site, int site_num, Branch *taxa)
{
    (*dataset)[taxa->get_taxa_id()].Assign_site(site_num, site);
}


void Simulate_from_PhyloMatrix::ascend_tree(int curr_site, int site_num, Branch *curr_branch, int rate)
{
    int new_site=0;
    double p_so_far, ran_p;
    
    
    if (curr_branch->get_brnlen()>FLOAT_TOL)
    {
        ran_p=ranf();
        
        p_so_far=curr_branch->get_trpb(rate, curr_site, new_site);
        
        while(p_so_far<ran_p && new_site<curr_model_matrix->get_num_states()-1)
        {
            new_site++;
            p_so_far+=curr_branch->get_trpb(rate, curr_site, new_site);
        }
        
    }
    else
    new_site=curr_site;
    
    if (curr_branch->is_tip()==TRUE)
    write_data(new_site, site_num, curr_branch);
    
    else
    {
        ascend_tree(new_site, site_num, curr_branch->get_child(0), rate);
        ascend_tree(new_site, site_num, curr_branch->get_child(1), rate);
    }
    
}  //End Simulate_from_PhyloMatrix::ascend_tree



int Simulate_from_PhyloMatrix::get_root_site()
{
    int i;
    double p_so_far, ran_p, step;
    
    ran_p=ranf();
    
    step = 1.0/(double)curr_model_matrix->get_num_root_states();
    i=0;
    p_so_far=0;
    
    while((i<curr_model_matrix->get_num_root_states()-1) &&(ran_p >p_so_far)) {
        p_so_far+=step;
        i++;
    }
    
    return(curr_model_matrix->get_nth_root_state(i)->get_state_id());
    
}
