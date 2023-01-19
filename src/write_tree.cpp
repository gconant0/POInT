#include <iostream>
#include <sstream>
#include <fstream>
#include "string.h"
#include "write_tree.h"

//#define DEBUG

using namespace::std;

void Write_Tree::write_tree(string output_file,  string prog_name, Tree *ctree, Exchange *cexchange)

{
    write_tree(output_file.c_str(), prog_name.c_str(), ctree, cexchange);
}

void Write_Tree::write_tree(const char *output_file, const char *prog_name, AA_matrix **the_mats,
			  Tree *ctree, Exchange *cexchange)

{

    treefile = new ofstream(output_file);
    
	curr_exchange=cexchange;
	tree_object=ctree;
	strcpy(filename, output_file);

	strcpy(program_name, prog_name);
	the_matrices=the_mats;
	write_descript_string();
    
    //treefile->open(filename);
    
	write_tree_string();
    //treefile->close();
    delete treefile;
}


void Write_Tree::write_tree(const char *output_file, const char *prog_name, Tree *ctree, Exchange *cexchange)
{

    treefile = new ofstream(output_file);
    
	curr_exchange=cexchange;
	tree_object=ctree;
	strcpy(filename,output_file);

	strcpy(program_name, prog_name);
	strcpy(matrix_file, "\0");
	write_descript_string();
    
    //treefile->open(filename);
    
	write_tree_string();
    //treefile->close();
    delete treefile;
}

void Write_Tree::write_tree_to_string(string &tree_string, string prog_name, Tree *ctree, Exchange *cexchange)
{
    ostringstream stringout;
    
    treefile = &stringout;
    
    curr_exchange=cexchange;
    tree_object=ctree;
    
    strcpy(program_name, prog_name.c_str());
    strcpy(matrix_file, "\0");
    dlines=0;
    //write_descript_string();
    
    
    write_tree_string();
   
    //cout<<"Tryed to write "<<stringout.str()<<" to string\n";
    
    tree_string=stringout.str();
    
    //delete treefile;
}


void Write_Tree::write_descript_string()
{
  int i, j, k, pos;
  char model_name[200], dummy[10], temp_num[40];

  
  switch(curr_exchange->get_model()) {
	case JC_NUCLEOTIDE:
		dlines=3;
		strcpy(model_name, "Jukes-Cantor 1969 model");
		break;
	case K2P_NUCLEOTIDE:
		dlines=4;
		strcpy(model_name, "Kimura 2-Parameter (1980) model");
		break;
	case HKY_NUCLEOTIDE:
		dlines=9;
		strcpy(model_name, "Hasegawa-Kishino-Yano (1985) model");
		break;
	case MG_94_JC:
	case MG_94_K2P:
	case MG_94_HKY:
		dlines=9+(curr_exchange->get_num_rates()*curr_exchange->get_num_p_non_syn());
		strcpy(model_name, "Muse-Gaut/Goldman-Yang 1994 model");
		break;
	case C_00_JC:
	case C_00_K2P:
	case C_00_HKY:
		dlines=9+(curr_exchange->get_num_rates()*
			(2*curr_exchange->get_num_nonsyn_patterns()+curr_exchange->curr_groups->get_num_groups()));
		strcpy(model_name, "Amino acid groups model");
		break;
	case LCAP_JC:
	case LCAP_K2P:
	case LCAP_HKY:
		dlines=9+curr_exchange->get_num_rates()*(curr_exchange->get_num_aa_props()+1)*curr_exchange->get_num_nonsyn_patterns();
		strcpy(model_name, "Linear Comb. of Amino Acid Props. model\0");
		break;
	case AAMULMAT_JC:
	case AAMULMAT_K2P:
	case AAMULMAT_HKY:
		dlines=9+curr_exchange->get_num_rates()*curr_exchange->get_num_matrices()+curr_exchange->get_num_matrices()*curr_exchange->get_num_p_non_syn();
		strcpy(model_name, "Log-Odds AA Sub. property model w/ multiple mats\0");
		break;
	case DUPL:
		dlines=3;
		strcpy(model_name, "Three-state duplicate-loss model");
		break;
	case DUPL_FIX:
		dlines=4;
		strcpy(model_name, "Four-state duplicate loss model with fixed duplicates");
		break;
	case DUPL_PARALLEL:
		dlines=4;
		strcpy(model_name, "Five-state duplicate loss model with parallel losses");
		break;
	case DUPL_PARALLEL_2_RATE:
		dlines=5;
		strcpy(model_name, "Five-state duplicate loss model with convergence and 2 rates of losses");
		break;
	case DUPL_PARALLEL_FIX:
		dlines=5;
		strcpy(model_name, "Six-state duplicate loss model with parallel losses and fixed duplicates");
		break;
	case DUPL_PARALLEL_FIX_SUBF:
		dlines=5;
		strcpy(model_name, "Six-state duplicate loss model with subfunctionalizing losses and fixed duplicates");
		break;
	case DUPL_SUBF_3_RATE:
		dlines=7;
		strcpy(model_name, "Six-state duplicate loss model with subfunction. losses, fixed dupls, and all rates");
		break;
	case DUPL_SUBF_ONLY:
		dlines=6;
		strcpy(model_name, "Six-state duplicate loss model with only subfunctionaliztion");
		break;
	case DUPL_2_RATE_NOSUBF:
		dlines = 6;
		strcpy(model_name, "Six-state duplicate loss model without subfunction. losses and with fixed dupls and all rates");
		break;
	case DUPL_SUBF_ALL_STATES:
		dlines=5;
		strcpy(model_name, "Six-state duplicate loss model with subfunctionalization and all states modeled");
		break;
	case DUPL_NOSTATE:
		dlines=4;
		strcpy(model_name, "Three-state duplicate-loss model with inferred tracking");
		break;
	case DUPL_FIX_NOSTATE:
		dlines=5;
		strcpy(model_name, "Four-state duplicate loss model with fixed duplicates and inferred tracking");
		break;
	case DUPL_PARALLEL_NOSTATE:
		dlines=5;
		strcpy(model_name, "Five-state duplicate loss model with parallel losses and inferred tracking");
		break;
	case DUPL_PARALLEL_2_RATE_NOSTATE:
		dlines=6;
		strcpy(model_name, "Five-state duplicate loss model with convergence, 2 rates of losses, and inferred tracking");
		break;
	case DUPL_PARALLEL_FIX_NOSTATE:
		dlines=6;
		strcpy(model_name, "Six-state duplicate loss model with parallel losses, fixed duplicates and inferred tracking");
		break;
	case DUPL_PARALLEL_FIX_SUBF_NOSTATE:
		dlines=6;
		strcpy(model_name, "Six-state duplicate loss model with subfunctionalizing losses, fixed duplicates and inferred tracking");
		break;
	case DUPL_SUBF_3_RATE_NOSTATE:
		dlines=8;
		strcpy(model_name, "Six-state duplicate loss model with subfunction. losses, fixed dupls, all rates, and inferred tracking");
		break;
	case DUPL_SUBF_3_RATE_ALLSTATES_NOSTATE:
		dlines = 8;
		strcpy(model_name, "Six-state duplicate loss model with subfunction. losses, fixed dupls, all rates, all states and inferred tracking");
		break;
	case DUPL_2_RATE_NOSUBF_NOSTATE:
		dlines = 7;
		strcpy(model_name, "Six-state duplicate loss model without subfunction. losses and with fixed dupls, all rates, and inferred tracking");
		break;
	case DUPL_2_RATE_NOSUBF_ALLSTATES_NOSTATE:
		dlines = 7;
		strcpy(model_name, "Six-state duplicate loss model without subfunction. losses and with fixed dupls, all rates, all states and inferred tracking");
		break;
	case DUPL_SUBF_ONLY_NOSTATE:
		dlines=7;
		strcpy(model_name, "Six-state duplicate loss model with only subfunctionaliztion and inferred tracking");
		break;
	case DUPL_SLOW_LOSS_CONV_FIX:
		dlines=7;
		strcpy(model_name, "Six-state model with slow losses from the convergent and pseudo-fixed duplicates");
		break;
	case DUPL_SLOW_LOSS_CONV_FIX_NOSTATE:
		dlines=8;
		strcpy(model_name, "Six-state model with slow losses from the convergent and pseudo-fixed duplicates and inferred tracking");
		break;
	case SNP_EVOL_MODEL:
		dlines=5;
		strcpy(model_name, "SNP model with differing rates of SNP loss/gain");
		break;
      case ORTHO_PARSIMONY:
          dlines=2;
          strcpy(model_name, "Single origin parsimony model for orthologs");
          break;
  }
  
  if (curr_exchange->using_generic_site_rate_probs() == TRUE) 
	dlines += curr_exchange->get_num_rates();
  
  cout<<"Writing tree for model: "<<model_name<<endl;

  switch (curr_exchange->get_model()) {
	case MG_94_JC:
	case C_00_JC:
	case LCAP_JC:
	case AAMULMAT_JC:
		dlines-=6;
		strcat(model_name, " (with Juke-Cantor nuc. subs)\0");
		break;
	case MG_94_K2P:
	case C_00_K2P:
	case LCAP_K2P:
	case AAMULMAT_K2P:
		dlines-=1;
		strcat(model_name, " (with Kimura 2-param. nuc. subs)\0");
		break;
	case MG_94_HKY:
	case C_00_HKY:
	case LCAP_HKY:
	case AAMULMAT_HKY:
		strcat(model_name, " (with Hasegawa-Kishino-Yano nuc. subs)\0");
		break;
  }
   
  cout<<"Model name: "<<model_name<<" lines: "<<dlines<<endl<<flush<<flush;
  descript=new char*[dlines];
  for(i=0; i<dlines; i++)
	descript[i]=new char[120];
  
  strcpy(descript[0], program_name);
  strcpy(descript[1], model_name);
  pos=2;

  switch (curr_exchange->get_model()) {
	case HKY_NUCLEOTIDE:
	case MG_94_HKY:
	case C_00_HKY:
	case LCAP_HKY:
	case AAMULMAT_HKY:	
		if (curr_exchange->fixed_basefreq()==FALSE)
			strcpy(descript[2], "Estimated Base frequencies: ");      
		else
			strcpy(descript[2], "Observed Base frequencies: ");
  
		for(i=0; i<4; i++)
			{
				strcpy(descript[3+i], "\t\t   ");      
    
				double_to_string (temp_num, 39, 5, curr_exchange->return_basefreq(i));
				dummy[0]=num_to_base(i);
				dummy[1] = '\0';
				strcat(descript[3+i], dummy);
				strcat(descript[3+i], " = ");
				strcat(descript[3+i], temp_num);
			}
		pos=7;
	case K2P_NUCLEOTIDE:
	case MG_94_K2P:
	case C_00_K2P:
	case LCAP_K2P:
	case AAMULMAT_K2P:
		strcpy(descript[pos], "Estimated ti/tv ratio = ");      
  
		double_to_string (temp_num, 39, 5, curr_exchange->get_obs_trs_trv());
		strcat(descript[pos], temp_num);
		strcat(descript[pos], " (kappa = ");
		double_to_string (temp_num, 39, 5, curr_exchange->get_trs_trv());
		strcat(descript[pos], temp_num);
		strcat(descript[pos], ")\0");
		pos++;
		break;
	}
 	


  switch(curr_exchange->get_model()) {
	case MG_94_JC:
	case MG_94_K2P:
	case MG_94_HKY:
		for(i=0; i<curr_exchange->get_num_p_non_syn(); i++) {
			for(j=0; j<curr_exchange->get_num_rates(); j++) {
				strcpy(descript[pos], "Prob. non-synonymous substitution");
				if(curr_exchange->get_num_p_non_syn() > 1) {
					strcat(descript[pos], " ");
					int_to_string(temp_num, 39, i);
					strcat(descript[pos], temp_num);
					strcat(descript[pos], ": ");
				}	
				else
					strcat(descript[pos], ": ");
				double_to_string (temp_num, 39, 5, curr_exchange->get_p_non_syn(i, j));
				strcat(descript[pos], temp_num);
				pos++;
			}
		}
		break;
	case C_00_JC:
	case C_00_K2P:
	case C_00_HKY:
		for(i=0; i<curr_exchange->get_num_nonsyn_patterns(); i++) {
			for(k=0; k<curr_exchange->get_num_rates(); k++) {
				j=0;
				while((*tree_object)[j]->get_nonsyn_pattern() != i) j++;

				strcpy(descript[pos], "Prob. Intra-group substitution");
				if(curr_exchange->get_num_nonsyn_patterns() > 1) {
					strcat(descript[pos], " ");
					int_to_string(temp_num, 39, i);
					strcat(descript[pos], temp_num);
					strcat(descript[pos], ": ");
				}	
				else
					strcat(descript[pos], ": ");
				double_to_string (temp_num, 39, 5, 
				curr_exchange->get_p_inter_group((*tree_object)[j]->get_p_intergroup_num(), k));
				strcat(descript[pos++], temp_num);
				strcpy(descript[pos], "Prob. Extra-group substitution");
				if(curr_exchange->get_num_nonsyn_patterns() > 1) {
					strcat(descript[pos], " ");
					int_to_string(temp_num, 39, i);
					strcat(descript[pos], temp_num);
					strcat(descript[pos], ": ");
				}	
				else
					strcat(descript[pos], ": ");
					double_to_string (temp_num, 39, 5, 
							curr_exchange->get_p_non_syn((*tree_object)[j]->get_p_nonsyn_num(), k));
			
				strcat(descript[pos++], temp_num);	
			}
		}
		for (i=0; i<curr_exchange->curr_groups->get_num_groups(); i++)
		{
			strcpy(descript[pos], "Groups: ");

			int_to_string(temp_num, 39, i);
			strcat(descript[pos], temp_num);
			strcat(descript[pos], ":\t");

			for(j=0; j<20; j++)
				if (curr_exchange->curr_groups->get_group(j)==i)
				{
					dummy[0]=num_to_aa(j);
					strcat(descript[pos], dummy);
					strcat(descript[pos], " ");
				}
			pos++;
		}

		break;
	case LCAP_JC:
	case LCAP_K2P:
	case LCAP_HKY:
		for(i=0; i<+curr_exchange->get_num_nonsyn_patterns(); i++) {
			for(k=0; k<curr_exchange->get_num_rates(); k++) {
				j=0;
				while((*tree_object)[j]->get_nonsyn_pattern() != i) j++;

				strcpy(descript[pos], "Chemical Composition factor");
				if(curr_exchange->get_num_nonsyn_patterns() > 1) {
					strcat(descript[pos], " ");
					int_to_string(temp_num, 39, i);
					strcat(descript[pos], temp_num);
					strcat(descript[pos], ": ");
				}	
				else
					strcat(descript[pos], ": ");
				double_to_string (temp_num, 39, 5, 
					curr_exchange->get_aa_prop_fac((*tree_object)[j]->get_aa_prop_num(CHEM_COMP), k, CHEM_COMP));
				strcat(descript[pos++], temp_num);


				strcpy(descript[pos], "Polarity Factor");
				if(curr_exchange->get_num_nonsyn_patterns() > 1) {
					strcat(descript[pos], " ");
					int_to_string(temp_num, 39, i);
					strcat(descript[pos], temp_num);
					strcat(descript[pos], ": ");
				}	
				else
					strcat(descript[pos], ": ");
				double_to_string (temp_num, 39, 5, 
				curr_exchange->get_aa_prop_fac((*tree_object)[j]->get_aa_prop_num(POLARITY), k, POLARITY));	
				strcat(descript[pos++], temp_num);
			
				strcpy(descript[pos], "Volume factor");    
				if(curr_exchange->get_num_nonsyn_patterns() > 1) {
					strcat(descript[pos], " ");
					int_to_string(temp_num, 39, i);
					strcat(descript[pos], temp_num);
					strcat(descript[pos], ": ");
				}	
				else
					strcat(descript[pos], ": ");
				double_to_string (temp_num, 39, 5, 
				curr_exchange->get_aa_prop_fac((*tree_object)[j]->get_aa_prop_num(VOLUME), k, VOLUME));
				strcat(descript[pos++], temp_num);
				
				strcpy(descript[pos], "Iso-electric point factor");
				if(curr_exchange->get_num_nonsyn_patterns() > 1) {
					strcat(descript[pos], " ");
					int_to_string(temp_num, 39, i);
					strcat(descript[pos], temp_num);
					strcat(descript[pos], ": ");
				}	
				else
					strcat(descript[pos], ": ");
				double_to_string (temp_num, 39, 5, 
				curr_exchange->get_aa_prop_fac((*tree_object)[j]->get_aa_prop_num(ISO_ELEC), k, ISO_ELEC));
				strcat(descript[pos++], temp_num);
			
				strcpy(descript[pos], "Hydropathy point factor");
				if(curr_exchange->get_num_nonsyn_patterns() > 1) {
					strcat(descript[pos], " ");
					int_to_string(temp_num, 39, i);
					strcat(descript[pos], temp_num);
					strcat(descript[pos], ": ");
				}	
				else
					strcat(descript[pos], ": ");
				double_to_string (temp_num, 39, 5, 
				curr_exchange->get_aa_prop_fac((*tree_object)[j]->get_aa_prop_num(HYDROPATHY), k, HYDROPATHY));
				strcat(descript[pos++], temp_num);
			
				strcpy(descript[pos], "Scaling factor");
				if(curr_exchange->get_num_nonsyn_patterns() > 1) {
					strcat(descript[pos], " ");
					int_to_string(temp_num, 39, i);
					strcat(descript[pos], temp_num);
					strcat(descript[pos], ": ");
				}	
				else
					strcat(descript[pos], ": ");
				double_to_string (temp_num, 39, 5, 
				curr_exchange->get_aa_prop_fac((*tree_object)[j]->get_aa_prop_num(SCALING), k, SCALING));
				strcat(descript[pos++], temp_num);
			}
		}
		break;
	case AAMULMAT_JC:
	case AAMULMAT_K2P:
	case AAMULMAT_HKY:
		for(i=0; i<curr_exchange->get_num_matrices(); i++) {
			for(j=0; j<curr_exchange->get_num_p_non_syn(); j++) {
				strcpy(descript[pos], "Prob. non-syn. substitution ignoring log-odds matrix: ");
				double_to_string (temp_num, 39, 5, curr_exchange->get_matrix_coeff(i,j));
				strcat(descript[pos++], temp_num);
			}
		}
		for(i=0; i<curr_exchange->get_num_matrices(); i++) {			
			strcpy(descript[pos], "Log-Odds Matrix ");
			int_to_string(temp_num, 39, i);
			strcat(descript[pos], temp_num);
			strcat(descript[pos], " used: ");
			strcat(descript[pos++], the_matrices[i]->get_matrix_name());
		}
		break;	
	case DUPL_SUBF_3_RATE:
	case DUPL_SUBF_3_RATE_NOSTATE:	
	case DUPL_SUBF_3_RATE_ALLSTATES_NOSTATE:
		strcpy(descript[pos], "Relative rate of duplicate fixation after converg.: ");
		double_to_string (temp_num, 39, 5, curr_exchange->get_fix_rate_scale());
		strcat(descript[pos++], temp_num);
	case DUPL_SUBF_ONLY:
	case DUPL_2_RATE_NOSUBF:
	case DUPL_SUBF_ONLY_NOSTATE:

	case DUPL_2_RATE_NOSUBF_NOSTATE:
	case DUPL_2_RATE_NOSUBF_ALLSTATES_NOSTATE:
	case DUPL_SLOW_LOSS_CONV_FIX:
	case DUPL_SLOW_LOSS_CONV_FIX_NOSTATE:
		strcpy(descript[pos], "Relative rate of duplicate loss after converg.: ");
		double_to_string (temp_num, 39, 5, curr_exchange->get_loss_rate_scale());
		strcat(descript[pos++], temp_num);
	case DUPL_FIX:
	case DUPL_PARALLEL_FIX:
	case DUPL_PARALLEL_FIX_SUBF:
	case DUPL_SUBF_ALL_STATES:
	case DUPL_SUBF_ALL_STATES_NOSTATE:
	case DUPL_FIX_NOSTATE:
	case DUPL_PARALLEL_FIX_NOSTATE:
	case DUPL_PARALLEL_FIX_SUBF_NOSTATE:
		cout<<"Position to write fixation: "<<pos<<endl;
		strcpy(descript[pos], "Instantaneous rate of duplicate fixation: ");
		double_to_string (temp_num, 39, 5, curr_exchange->get_dupl_fix_rate());
		strcat(descript[pos++], temp_num);

	case DUPL_PARALLEL:
	case DUPL_PARALLEL_NOSTATE:
	case DUPL_PARALLEL_2_RATE:
	case DUPL_PARALLEL_2_RATE_NOSTATE:
		if ((curr_exchange->get_model() != DUPL_FIX) && (curr_exchange->get_model() != DUPL_FIX_NOSTATE)) { 
			strcpy(descript[pos], "Instantaneous rate of convergence to single duplicate: ");
			double_to_string (temp_num, 39, 5, curr_exchange->get_dupl_parallel_rate());
			strcat(descript[pos++], temp_num);
		}
		break;
	case SNP_EVOL_MODEL:
		strcpy(descript[pos], "Ratio of SNP gains to losses: ");
		double_to_string(temp_num, 39, 5, curr_exchange->get_snp_gains_to_losses());
		strcat(descript[pos++], temp_num);
		strcpy(descript[pos], "Ratio of SNP losses to SNP state changes: ");
		double_to_string(temp_num, 39, 5, curr_exchange->get_snp_loss_to_switch());
		strcat(descript[pos++], temp_num);
		break;
	}

	//Second pass for model params
	switch(curr_exchange->get_model()) {
		case DUPL_PARALLEL_2_RATE:
		case DUPL_PARALLEL_2_RATE_NOSTATE:
			strcpy(descript[pos], "Relative rate of duplicate loss after converg.: ");
			double_to_string (temp_num, 39, 5, curr_exchange->get_loss_rate_scale());
			strcat(descript[pos++], temp_num);
		break;
		case DUPL_SLOW_LOSS_CONV_FIX:
		case DUPL_SLOW_LOSS_CONV_FIX_NOSTATE:
			strcpy(descript[pos], "Relative rate of duplicate loss after fixation: ");
			double_to_string (temp_num, 39, 5, curr_exchange->get_fix_loss_rate());
			strcat(descript[pos++], temp_num);
			break;
	}


	switch(curr_exchange->get_model()) {
	case DUPL_NOSTATE:
	case DUPL_FIX_NOSTATE:
	case DUPL_PARALLEL_NOSTATE:
	case DUPL_PARALLEL_FIX_NOSTATE:
	case DUPL_PARALLEL_FIX_SUBF_NOSTATE:
	case DUPL_SUBF_ALL_STATES_NOSTATE:
	case DUPL_SUBF_ONLY_NOSTATE:
	case DUPL_SUBF_3_RATE_ALLSTATES_NOSTATE:
	case DUPL_SUBF_3_RATE_NOSTATE:
	case DUPL_2_RATE_NOSUBF_NOSTATE:
	case DUPL_2_RATE_NOSUBF_ALLSTATES_NOSTATE:
	case DUPL_PARALLEL_2_RATE_NOSTATE:
	case DUPL_SLOW_LOSS_CONV_FIX_NOSTATE:
		strcpy(descript[pos], "Track switch between genes probability: ");
		double_to_string (temp_num, 39, 5, curr_exchange->get_strand_switch_prob());
		strcat(descript[pos++], temp_num);
		break;

		
  }
  
  if(curr_exchange->using_generic_site_rate_probs() == TRUE) {
	  for(i=0; i<curr_exchange->get_num_rates(); i++) {
		strcpy(descript[pos], "Rate ");
		int_to_string(temp_num, 39, i);
		strcat(descript[pos], temp_num);
		strcat(descript[pos], " probability: ");
		double_to_string (temp_num, 39, 5, curr_exchange->get_rate_prob(i));
		strcat(descript[pos++], temp_num);
	  }
		
  }

    if (curr_exchange->get_model() != ORTHO_PARSIMONY) {
        strcpy(descript[pos], "Final ln likelihood: ");
        double_to_string (temp_num, 39, 4, curr_exchange->get_saved_lnL());
        strcat(descript[pos], temp_num);
    }
}


void Write_Tree::decend_clade(Branch *current)
{
  cerr<<"Wrong decend_clade\n";
}






void Write_Tree::get_sibling(Branch *sib)
{
  
}


void Write_Nexus_Tree::write_tree_string()
{
  int i=0,j ;

  
  *treefile<<"#NEXUS"<<endl<<endl<<"Begin trees; \n [!\n";
    for (i=0; i<dlines; i++) {
#ifdef DEBUG
        cout<<">"<<descript[i]<<endl;
#endif
        *treefile<<">"<<descript[i]<<endl;
    }
  
  *treefile<<"]\n	Translate\n";
	
  for(i=0; i<curr_exchange->get_num_taxa(); i++)
    { 
	  j=0;
	  while((*tree_object)[j]->get_taxa_id() != i) j++;
	  *treefile<<"\t"<<i+1<<" \'"<<(*tree_object)[j]->get_name()<<"\'";
	  
	  if (i+1!=curr_exchange->get_num_taxa())
			*treefile<<",\n";
	  else
			*treefile<<"\n\t;\n";
	  
    }
  
  *treefile<<"\ntree PAUP_1 = [&R] "<<flush;
  decend_clade(tree_object->find_root());
 
  *treefile<<";\nEnd;\n";
 
}



void Write_Nexus_Tree::decend_clade(Branch *current)
{

  if (current->is_tip()==TRUE)
    *treefile<<current->get_taxa_id()+1<<":"<<current->expect_subs_site();
  else
    {
      *treefile<<"(";
      decend_clade(current->get_child(0));
    }
  if (current!=tree_object->find_root())
    {
      *treefile<<",";
      get_sibling(current->get_sibling());
      *treefile<<"):"<<current->get_parent()->expect_subs_site();
    }
 
}



void Write_Nexus_Tree::get_sibling(Branch *sib)
{
  if (sib->is_tip()==TRUE)
    *treefile<<sib->get_taxa_id()+1<<":"<<sib->expect_subs_site();
  else
    {
      *treefile<<"(";
      decend_clade(sib->get_child(0));
    }

}



void Write_GCC_Tree::decend_clade(Branch *current)
{

	if (current->is_tip()==TRUE) {
		*treefile<<current->get_taxa_id()+1<<":"<<current->expect_subs_site();
		write_branch_params(current);		
	}
  else
    {
      *treefile<<"(";
      decend_clade(current->get_child(0));
    }
  if (current!=tree_object->find_root())
    {
      *treefile<<",";
      get_sibling(current->get_sibling());
      *treefile<<"):"<<current->get_parent()->expect_subs_site();
	  write_branch_params(current->get_parent());
    }
 
}



void Write_GCC_Tree::get_sibling(Branch *sib)
{
	if (sib->is_tip()==TRUE) {
		*treefile<<sib->get_taxa_id()+1<<":"<<sib->expect_subs_site();
		write_branch_params(sib);
	}
	else
    {
      *treefile<<"(";
      decend_clade(sib->get_child(0));
    }

}



void Write_GCC_Tree::write_branch_params(Branch *current)
{
	int i;

	switch(curr_exchange->get_model())
		{
		case MG_94_JC:
		case MG_94_K2P:
		case MG_94_HKY:	
			*treefile<<":"<<current->get_p_nonsyn_num()+1;
			break;
		case C_00_JC:
		case C_00_K2P:
		case C_00_HKY:
			*treefile<<":"<<current->get_p_nonsyn_num()+1;
			*treefile<<":"<<current->get_p_intergroup_num()+1;
			break;
		case LCAP_JC:
		case LCAP_K2P:
		case LCAP_HKY:
			for(i=0; i<curr_exchange->get_num_aa_props()+1; i++)
				*treefile<<":"<<current->get_aa_prop_num(i)+1;
			break;
		case AAMULMAT_JC:
		case AAMULMAT_K2P:
		case AAMULMAT_HKY:
			for(i=0; i<curr_exchange->get_num_matrices(); i++)
				*treefile<<":"<<current->get_matrix_coeff_num(i)+1;
			break;
	}
}



void Write_Phylip_Tree::write_tree_string()
{
    int i=0,j ;
   
    decend_clade(tree_object->find_root());
   
   
}



void Write_Phylip_Tree::decend_clade(Branch *current)
{
    
    if (current->is_tip()==TRUE)
        *treefile<<current->get_name()<<":"<<current->expect_subs_site();
    else
    {
        *treefile<<"(";
        decend_clade(current->get_child(0));
    }
    if (current!=tree_object->find_root())
    {
        *treefile<<",";
        get_sibling(current->get_sibling());
        *treefile<<"):"<<current->get_parent()->expect_subs_site();
    }
    
}



void Write_Phylip_Tree::get_sibling(Branch *sib)
{
    if (sib->is_tip()==TRUE)
        *treefile<<sib->get_name()<<":"<<sib->expect_subs_site();
    else
    {
        *treefile<<"(";
        decend_clade(sib->get_child(0));
    }
    
}
