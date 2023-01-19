#include "search_trees.h"

using namespace::std;

Tree_searcher::Tree_searcher(Exchange *cexchange)
{
	start=end=-1;
	local_model=0;
	save_num=1;
	num_saved=0;
	curr_code=0;
	curr_exchange=cexchange;
    the_matrix=0;
	allocate();
}



Tree_searcher::Tree_searcher(Exchange *cexchange, int n)
{
	start=end=-1;
	local_model=0;
	save_num=n;
	num_saved=0;
	curr_code=0;
	curr_exchange=cexchange;
    the_matrix=0;
	allocate();
}


Tree_searcher::Tree_searcher(Exchange *cexchange, Genetic_code *ccode)
{
	start=end=-1;
	local_model=0;
	save_num=1;
	num_saved=0;
	curr_code=ccode;
	curr_exchange=cexchange;
    the_matrix=0;
	allocate();
}


Tree_searcher::Tree_searcher(Exchange *cexchange, Genetic_code *ccode, int n)
{
	start=end=-1;
	local_model=0;
	save_num=n;
	num_saved=0;
	curr_code=ccode;
	curr_exchange=cexchange;
    the_matrix=0;
	allocate();
}


Tree_searcher::Tree_searcher(Exchange *cexchange, Clade *cgenomes, WGD_Data *chomologs)
{
	start=end=-1;
	local_model=0;
	save_num=1;
	num_saved=0;
	the_genomes=cgenomes;
	the_homologs=chomologs;
	curr_exchange=cexchange;
    the_matrix=0;
	allocate();
}


Tree_searcher::Tree_searcher(Exchange *cexchange, Clade *cgenomes, WGD_Data *chomologs, char * post_file)
{
	start=end=-1;
	local_model=0;
	save_num=1;
	num_saved=0;
	the_genomes=cgenomes;
	the_homologs=chomologs;
	curr_exchange=cexchange;
	strcpy(post_prob_file, post_file);
    the_matrix=0;
	allocate();
}


Tree_searcher::Tree_searcher(Exchange *cexchange, Clade *cgenomes, WGX_Data *chomologs, Phylo_Matrix *cmatrix)
{
    start=end=-1;
    local_model=0;
    save_num=1;
    num_saved=0;
    the_genomes=cgenomes;
    the_WGX_homologs=chomologs;
    the_matrix=cmatrix;
    curr_exchange=cexchange;
    allocate();
}


Tree_searcher::Tree_searcher(Exchange *cexchange, Clade *cgenomes, WGX_Data *chomologs, Phylo_Matrix *cmatrix, int n)
{
    start=end=-1;
    local_model=0;
    save_num=n;
    num_saved=0;
    the_genomes=cgenomes;
    the_WGX_homologs=chomologs;
    the_matrix=cmatrix;
    curr_exchange=cexchange;
    allocate();
}

void Tree_searcher::write_out_trees(int num, char *partial_file_name, const char *prog_name)
{
	int i, j;
	char filename[200], temp_num[40];
	Write_Nexus_Tree *writeout_tree;
  

	if (num>num_saved) {
		cerr<<"ERROR: "<<num<<" treefiles were requested when there are only "<<num_saved<<" in memory\n";
		num=num_saved;
	}
	for(i=0; i<num; i++) {
	//	for(j=0; j<curr_exchange->get_num_branches(); j++)
	//		(*best_trees[i])[j]->set_expect_subs_site((*best_trees[i])[j]->get_brnlen());

		writeout_tree=new Write_Nexus_Tree();
		strcpy(filename, partial_file_name);
		int_to_string (temp_num, 39, i);
		strcat(filename, temp_num);
		strcat(filename, ".tre");

		writeout_tree->write_tree(filename, prog_name, best_trees[i], best_exchanges[i]);
		delete writeout_tree;
	}

	cout<<"Wrote "<<num<<" treefiles\n";
}


void Tree_searcher::write_post_probs()
{
	Exchange *new_exchange;
	Tree *new_tree;

	new_tree=new Tree(best_exchanges[0], curr_exchange->is_rooted_tree());
	new_exchange=new Exchange();

	(*new_tree)=(*best_trees[0]);
	(*new_exchange)=(*best_exchanges[0]);

	switch(best_exchanges[0]->get_model())
      {

		case DUPL_NOSTATE:
			 local_model = new Dupl_NoState_model (new_exchange, new_tree, the_genomes, the_homologs);
			 break;   
		case DUPL_FIX_NOSTATE:
			local_model = new Dupl_Fix_NoState_model(new_exchange, new_tree, the_genomes, the_homologs);
			break;
		case DUPL_PARALLEL_NOSTATE:
			local_model = new Dupl_Parallel_NoState_model(new_exchange, new_tree, the_genomes, the_homologs);
			break;
		case DUPL_PARALLEL_FIX_NOSTATE:
			local_model = new Dupl_Fix_Parallel_NoState_model(new_exchange, new_tree, the_genomes, the_homologs);
			break;
		case DUPL_PARALLEL_FIX_SUBF_NOSTATE:
			local_model = new Dupl_Fix_Parallel_SubF_NoState_model(new_exchange, new_tree, the_genomes, the_homologs);
			break;  
		case DUPL_2_RATE_NOSUBF_NOSTATE:
			  local_model = new Dupl_2_Rate_NoSubF_NoState_model (new_exchange, new_tree, the_genomes, the_homologs);
			  break;
	}

	//The branch lengths get overwritten as the like_model is build--fix
	(*new_tree)=(*best_trees[0]);
	(*new_exchange)=(*best_exchanges[0]);

	local_model->recalculate_transprobs();

	cout<<"Likelihood for post probs is "<<local_model->find_appropriate_ln_like()<<endl;

	local_model->print_tracking_probs(post_prob_file);

	delete new_tree;


}

Tree_searcher::~Tree_searcher()
{
	int i;

	if (local_model !=0)
		delete local_model;

	if (save_num > 0) {
	
		//First delete the objects themselves
		for(i=0; i<save_num; i++) {
			if (best_trees[i] != 0)
				delete best_trees[i];
			if(best_exchanges[i] != 0)
				delete best_exchanges[i];
		}

		delete[] best_trees;
		delete[] best_exchanges;
	}
}


void Tree_searcher::allocate() 
{
	int i;

	treenum=0;
	best_scores=new double[save_num];
	best_trees=new Tree*[save_num];
	best_exchanges=new Exchange* [save_num];
	best_num=new int [save_num];
    
	for(i=0; i<save_num; i++) {
		best_trees[i]=0;
		best_exchanges[i]=0;
		best_scores[i]=0;
		best_num[i]=-1;
	}
}



void Tree_searcher::evaluate_tree(Tree *tree_to_run, Exchange *exchange_to_run)

//Takes an input tree, optimizes it w/ the Powell code
//and saves it if it is one of the total save_num found
//so far.
{
	int i, j, k, tempnum;
	double temp_score, *temp_params=0;

    
    if (the_matrix !=0) temp_params = new double [the_matrix->get_num_params()];
	//Pointer variables to allow us to insert trees in the lists
	Exchange *tempexchange;
	Tree *temptree;

	if ((start == -1) || ((treenum >= start) && (treenum <= end))) {
		eval_tree=tree_to_run;
		eval_exchange=exchange_to_run;
		//cout<<"Eval excahgne has "<<(*eval_exchange->get_dataset())[0].Sequence_size()<<" sites"<<endl;

		for(i=0; i<curr_exchange->get_num_branches(); i++) 
			(*eval_tree)[i]->set_expect_subs_site(0.3);

		switch(eval_exchange->get_model())
		  {

			case LCAP_JC:
				eval_model = new LCAP_JC_model (eval_exchange, eval_exchange->get_dataset(), eval_tree, curr_code);
				break;
			case LCAP_K2P:
				eval_model = new LCAP_K2P_model (eval_exchange, eval_exchange->get_dataset(), eval_tree, curr_code);
				break;
			case LCAP_HKY:
				eval_model = new LCAP_HKY_model (eval_exchange, eval_exchange->get_dataset(), eval_tree, curr_code);
				break;
			case C_00_JC:
				eval_model = new C_00_JC_model (eval_exchange, eval_exchange->get_dataset(), eval_tree, curr_code);
				break;
			case C_00_K2P:
				eval_model = new C_00_K2P_model (eval_exchange, eval_exchange->get_dataset(), eval_tree, curr_code);
				break;
			case C_00_HKY:
				eval_model = new C_00_HKY_model (eval_exchange, eval_exchange->get_dataset(), eval_tree, curr_code);
				break;
			case MG_94_JC:
				eval_model = new MG_94_JC_model (eval_exchange, eval_exchange->get_dataset(), eval_tree, curr_code);
				break;
			case MG_94_K2P:
				eval_model = new MG_94_K2P_model (eval_exchange, eval_exchange->get_dataset(), eval_tree, curr_code);
				break;
			case MG_94_HKY:
				eval_model = new MG_94_HKY_model (eval_exchange, eval_exchange->get_dataset(), eval_tree, curr_code);
				break;
			//Ignore for the momement
			/*case AAMULMAT_JC:
				eval_model = new MultiMatrix_JC_model (eval_exchange, eval_exchange->get_dataset(), eval_tree, curr_code, the_matrices);
				break;
			case AAMULMAT_K2P:
				eval_model = new MultiMatrix_K2P_model (eval_exchange, eval_exchange->get_dataset(), eval_tree, curr_code, the_matrices);
				break;
			case AAMULMAT_HKY:
				eval_model = new MultiMatrix_HKY_model (eval_exchange, eval_exchange->get_dataset(), eval_tree, curr_code, the_matrices);
				break;
			*/
			case DUPL:
				 eval_model = new Dupl_model (eval_exchange, eval_exchange->get_dataset(), eval_tree);
				 break;   
			case DUPL_FIX:
				eval_model = new Dupl_Fix_model(eval_exchange, eval_exchange->get_dataset(), eval_tree);
				break;
			case DUPL_PARALLEL:
				eval_model = new Dupl_Parallel_model(eval_exchange, eval_exchange->get_dataset(), eval_tree);
				break;
			case DUPL_PARALLEL_FIX:
				eval_model = new Dupl_Fix_Parallel_model(eval_exchange, eval_exchange->get_dataset(), eval_tree);
				break;
			case DUPL_PARALLEL_FIX_SUBF:
				eval_model = new Dupl_Fix_Parallel_SubF_model(eval_exchange, eval_exchange->get_dataset(), eval_tree);
				break;
			case DUPL_SLOW_LOSS_CONV_FIX:
				eval_model = new Dupl_Slow_Loss_Con_Fix_model(eval_exchange, eval_exchange->get_dataset(), eval_tree);
				break;
			case DUPL_NOSTATE:
				 eval_model = new Dupl_NoState_model (eval_exchange, eval_tree, the_genomes, the_homologs);
				 break;   
			case DUPL_FIX_NOSTATE:
				eval_model = new Dupl_Fix_NoState_model(eval_exchange, eval_tree, the_genomes, the_homologs);
				break;
			case DUPL_PARALLEL_NOSTATE:
				eval_model = new Dupl_Parallel_NoState_model(eval_exchange, eval_tree, the_genomes, the_homologs);
				break;
			case DUPL_PARALLEL_FIX_NOSTATE:
				eval_model = new Dupl_Fix_Parallel_NoState_model(eval_exchange, eval_tree, the_genomes, the_homologs);
				break;
			case DUPL_PARALLEL_FIX_SUBF_NOSTATE:
				eval_model = new Dupl_Fix_Parallel_SubF_NoState_model(eval_exchange,  eval_tree, the_genomes, the_homologs);
				break; 
			case DUPL_SUBF_ONLY_NOSTATE:
				eval_model = new Dupl_SubF_Only_NoState_model(eval_exchange,  eval_tree, the_genomes, the_homologs);
				break; 
			case DUPL_SUBF_3_RATE_NOSTATE:
				eval_model = new Dupl_SubF_3_Rate_NoState_model(eval_exchange,  eval_tree, the_genomes, the_homologs);
				break;
			case DUPL_2_RATE_NOSUBF_NOSTATE:
				eval_model = new Dupl_2_Rate_NoSubF_NoState_model(eval_exchange,  eval_tree, the_genomes, the_homologs);
				break;
			case DUPL_SLOW_LOSS_CONV_FIX_NOSTATE:
				eval_model= new Dupl_Slow_Loss_Con_Fix_NoState_model(eval_exchange, eval_tree, the_genomes, the_homologs);
				break;
              //This model is created in setup to give access to its functions
                  //case DUPL_ARBITRARY:
                  //eval_model = new Ploidy_Like_model(eval_exchange, eval_tree, the_genomes, the_WGX_homologs, the_matrix);
                //break;
		}
        setup();
		eval_powell=new Powell();
		
        //eval_exchange->set_saved_lnL(eval_model->find_appropriate_ln_like());
        
		//cout<<"Begining global optimization: there are "<<current_exchange.get_num_params()<<" parameters\n"<<flush;
	   //cout<<"Exchng parms: "
		//<<" DP: "<<eval_exchange->get_dupl_parallel_rate()<<" DF: "<<eval_exchange->get_dupl_fix_rate()<<" SP: "<<eval_exchange->get_strand_switch_prob()<<" CL: "<<eval_exchange->get_loss_rate_scale()
		//<<" FRS: "<<eval_exchange->get_fix_rate_scale()<<" FLR: "<<eval_exchange->get_fix_loss_rate()<<endl;
		

        
	   // if (treenum > 10)
			eval_exchange->set_saved_lnL(eval_powell->Init_min(eval_model, eval_exchange, FALSE));
		cout<<"Evalulated tree.  lnL: "<<eval_exchange->get_saved_lnL()<<": "<<treenum<<endl;
		

		for(i=0; i<eval_exchange->get_num_branches(); i++)
				(*eval_tree)[i]->set_expect_subs_site(eval_model->get_expect_sub_site((*eval_tree)[i]->get_brnlen()));


		if (num_saved > 0) {
		//See if this tree is one of save_num best
				i=0;
				while((i<num_saved) && (best_scores[i]>eval_exchange->get_saved_lnL())) i++;
				if(i<num_saved) {
				//This tree belongs at position i in the list
					temp_score=best_scores[i];
					temptree=best_trees[i];
					tempexchange=best_exchanges[i];
					tempnum=best_num[i];
                    
                    if (the_matrix !=0) {
                        for(j=0; j<the_matrix->get_num_params(); j++) temp_params[j]=the_matrix->get_scaled_param(j);
                    }
				
					//Add an element if needed
					if (num_saved < save_num) {
						best_scores[num_saved]=best_scores[num_saved-1];
						best_trees[num_saved]=best_trees[num_saved-1];
						best_exchanges[num_saved]=best_exchanges[num_saved-1];
						best_num[num_saved]=best_num[num_saved-1];
                        
                        if (the_matrix !=0) {
                            for(j=0; j<the_matrix->get_num_params(); j++) best_params[num_saved][j]=best_params[num_saved-1][j];
                        }
                        
					}
					else
						//Clear the last element
					{
						delete best_exchanges[num_saved-1];
						delete best_trees[num_saved-1];
					}


					best_scores[i]=eval_exchange->get_saved_lnL();
					best_exchanges[i]=new Exchange();
					(*best_exchanges[i])=(*eval_exchange);
					//best_trees[i]=new Tree(best_exchanges[i], best_exchanges[i]->is_rooted_tree());
					//(*best_trees[i])=(*eval_tree);
                    save_tree(i, eval_tree);
                    best_num[i]=treenum;

                    if (the_matrix !=0) {
                        for(j=0; j<the_matrix->get_num_params(); j++) best_params[i][j]=the_matrix->get_scaled_param(j);
                    }

					for(j=num_saved-1; j>i+1; j--) {
						best_scores[j]=best_scores[j-1];
						best_trees[j]=best_trees[j-1];
						best_exchanges[j]=best_exchanges[j-1];
						best_num[j]=best_num[j-1];
                        
                        if (the_matrix !=0) {
                            for(k=0; k<the_matrix->get_num_params(); k++) best_params[j][k]=best_params[j-1][k];
                        }

					}
					if (i+1 < save_num) {
						best_scores[i+1]=temp_score;
						best_trees[i+1]=temptree;
						best_exchanges[i+1]=tempexchange;
						best_num[i+1]=tempnum;
                        
                        if (the_matrix !=0) {
                            for(k=0; k<the_matrix->get_num_params(); k++) best_params[i+1][k]=temp_params[k];
                        }
					}

					if (num_saved<save_num) num_saved++;
				}
			
				else if (num_saved<save_num) {
					//This list is incomplete: add this tree to it
					best_exchanges[num_saved]=new Exchange();
					(*best_exchanges[num_saved])=(*eval_exchange);
					//best_trees[num_saved]=new Tree(best_exchanges[num_saved], best_exchanges[num_saved]->is_rooted_tree());
					//(*best_trees[num_saved])=(*eval_tree);
                    save_tree(num_saved, eval_tree);
					best_scores[num_saved]=eval_exchange->get_saved_lnL();
					best_num[num_saved]=treenum;
                    if (the_matrix !=0) {
                        for(j=0; j<the_matrix->get_num_params(); j++) best_params[num_saved][j]=the_matrix->get_scaled_param(j);
                    }
                    
					num_saved++;
				}


		}
		else {
		//This is the first tree we've found: save it
			best_exchanges[0]=new Exchange();
			(*best_exchanges[0])=(*eval_exchange);
			//best_trees[0]=new Tree(best_exchanges[0], best_exchanges[num_saved]->is_rooted_tree());
			//(*best_trees[0])=(*eval_tree);
            save_tree(0, eval_tree);
			best_scores[0]=eval_exchange->get_saved_lnL();
			best_num[0]=treenum;
            if (the_matrix !=0) {
                for(j=0; j<the_matrix->get_num_params(); j++) best_params[0][j]=the_matrix->get_scaled_param(j);
            }
            
			num_saved++;
		}
		cout<<"Best scores: ";
		for(i=0; i<save_num; i++)
			cout<<"\t"<<best_scores[i]<<", "<<best_num[i];
		cout<<endl;

		//Delete used objects
		delete eval_powell;
        if (temp_params!=0)
            delete[] temp_params;
		//delete eval_model;
	}
	treenum++;
}

void Tree_searcher::save_tree(int tree_num, Tree * source_tree)
{
    best_trees[tree_num]=new Tree(best_exchanges[tree_num], curr_exchange->is_rooted_tree());
    (*best_trees[tree_num])=(*source_tree);
}

Exhaustive_Tree_searcher::Exhaustive_Tree_searcher(Exchange *cexchange) : Tree_searcher(cexchange)
{
}
	
Exhaustive_Tree_searcher::Exhaustive_Tree_searcher(Exchange *cexchange, int n) :
Tree_searcher(cexchange, n)
{
}


Exhaustive_Tree_searcher::Exhaustive_Tree_searcher(Exchange *cexchange, Genetic_code *ccode) :
Tree_searcher(cexchange, ccode)
{
}


Exhaustive_Tree_searcher::Exhaustive_Tree_searcher(Exchange *cexchange, Genetic_code *ccode, int n) :
Tree_searcher(cexchange, ccode, n)
{
}

Exhaustive_Tree_searcher::Exhaustive_Tree_searcher(Exchange *cexchange, Clade *cgenomes, WGD_Data *chomologs) :
Tree_searcher(cexchange, cgenomes, chomologs)
{
}

Exhaustive_Tree_searcher::Exhaustive_Tree_searcher(Exchange *cexchange, Clade *cgenomes, WGD_Data *chomologs, char * post_file) :
Tree_searcher(cexchange, cgenomes, chomologs, post_file)
{
}




void Exhaustive_Tree_searcher::start_search()
{
	assemble_trees(0, 0, 0);
}



void Exhaustive_Tree_searcher::assemble_trees(Tree *partial_tree, Exchange *partial_exchange, int num_added)
{
	int i, j, num_brns, childnum, brncnt;
	Branch *new_brn, *new_parent, *old_parent, *old_sibling, **new_added_brns;
	Tree *new_tree;
	Exchange *new_exchange;

	if (num_added == 0) 
	{
		//Setup the tip branches now--indexed 0 through num_taxa-1
		new_exchange=new Exchange();
		(*new_exchange)=(*curr_exchange);
		new_tree=new Tree(new_exchange, new_exchange->is_rooted_tree());
	
		for(i=0; i<new_exchange->get_num_taxa(); i++) {
			(*new_tree)[i]->set_name((*new_exchange->get_dataset())[i].Sequence_name());
			(*new_tree)[i]->set_tip(TRUE);
			(*new_tree)[i]->set_taxa_id(i);
		}
	
		(*new_tree)[0]->initialized();
		(*new_tree)[1]->initialized();
		new_tree->set_as_siblings((*new_tree)[0], (*new_tree)[1]);
		new_parent=(*new_tree)[new_exchange->get_num_taxa()];
		new_parent->initialized();
		new_tree->set_as_parent_child(new_parent, (*new_tree)[0], 0);
		new_tree->set_as_parent_child(new_parent, (*new_tree)[1], 1);
		new_tree->set_root(new_parent);
		new_added_brns=new Branch *[3];

		new_added_brns[0]=(*new_tree)[0];
		new_added_brns[1]=(*new_tree)[1];
		new_added_brns[2]=new_parent;

		//Recursively assemble all possible trees
		assemble_trees(new_tree, new_exchange, 2);

		//Clean up
		delete new_tree;
		delete new_exchange;
		delete[] new_added_brns;
	}
	else if (num_added < partial_exchange->get_num_taxa())
	{
		num_brns=2*num_added-1;

	

		//Add the new tip to every possible branch
		for(i=0; i<num_brns; i++) {
			
			new_exchange=new Exchange();
			(*new_exchange)=(*partial_exchange);
			new_tree=new Tree(new_exchange, new_exchange->is_rooted_tree());
				
			//Copy the tree as assembled so far into a new object
			(*new_tree)=(*partial_tree);
					
			new_added_brns=new Branch*[2*(num_added+1)-1];
			j=brncnt=0;
			while(j<new_exchange->get_num_branches()) {
				if ((*new_tree)[j]->is_uninitialized() == FALSE) 
					new_added_brns[brncnt++]=(*new_tree)[j];
				j++;
			}
			if (brncnt != num_brns)
				cerr<<"ERROR: Did not reassemble array correctly\n";
				
			//For unrooted trees, we don't add to the null branch or the root
			if ((new_exchange->is_rooted_tree() == TRUE) || 
			((new_added_brns[i] != new_tree->find_null_branch_id()) && (new_added_brns[i] != new_tree->find_root())) ) {
			
				new_brn=(*new_tree)[num_added];
				new_brn->initialized();

                
				//The next "free" branch is after all n taxa + s-1 of the already added branches
				new_parent=(*new_tree)[num_added-1+new_exchange->get_num_taxa()];
					
				if (new_parent->is_uninitialized() == FALSE)						
					cerr<<"ERROR: Bad ordering of interior branches!\n";


				new_parent->initialized();

				old_sibling=new_added_brns[i]->get_sibling();
				old_parent=new_added_brns[i]->get_parent();

				//Could be adding at the root
				if (old_parent != 0) {
					if (old_sibling==old_parent->get_child(0)) childnum=1;
					else childnum=0;
				
					new_tree->set_as_parent_child(old_parent, new_parent, childnum);
					new_tree->set_as_siblings(old_sibling, new_parent);
				}

				new_tree->set_as_parent_child(new_parent, new_brn, 0);
				new_tree->set_as_parent_child(new_parent, new_added_brns[i], 1);
				new_tree->set_as_siblings(new_brn, new_added_brns[i]);					
				new_tree->set_root(new_parent);
				new_added_brns[num_brns]=new_brn;
				new_added_brns[num_brns+1]=new_parent;				
					

				//Recursively assemble all possible trees
				assemble_trees(new_tree, new_exchange, num_added+1);

				//Clean up
				delete new_tree;
				delete new_exchange;
				delete[] new_added_brns;	
			}
			else {
				delete new_tree;
				delete new_exchange;
				delete[] new_added_brns;
			}
		}

	}
    else {
		//We have a complete tree--optimize it
       
		evaluate_tree(partial_tree, partial_exchange);
    }

}


Exhaustive_Tree_PhyloMat_searcher::Exhaustive_Tree_PhyloMat_searcher(Exchange *cexchange, Clade *cgenomes, WGX_Data *chomologs, Phylo_Matrix *cmatrix) :
Tree_searcher(cexchange, cgenomes, chomologs, cmatrix)
{
    PhyloMat_allocate();
}

Exhaustive_Tree_PhyloMat_searcher::Exhaustive_Tree_PhyloMat_searcher(Exchange *cexchange, Clade *cgenomes, WGX_Data *chomologs, Phylo_Matrix *cmatrix, int n) :
Tree_searcher(cexchange, cgenomes, chomologs, cmatrix, n)
{
    PhyloMat_allocate();
}

void Exhaustive_Tree_PhyloMat_searcher::start_search()
{
    assemble_trees(0, 0, 0);
}



void Exhaustive_Tree_PhyloMat_searcher::assemble_trees(Tree_Ex *partial_tree, Exchange *partial_exchange, int num_added)
{
    int i, j, num_brns, childnum, brncnt;
    Branch_Ex *new_brn, *new_parent, *old_parent, *old_sibling, **new_added_brns;
    Tree_Ex *new_tree;
    Exchange *new_exchange;
    
    if (num_added == 0)
    {
        //Setup the tip branches now--indexed 0 through num_taxa-1
        new_exchange=new Exchange();
        (*new_exchange)=(*curr_exchange);
        new_tree=new Tree_Ex(new_exchange, the_matrix, new_exchange->is_rooted_tree());
        
        for(i=0; i<new_exchange->get_num_taxa(); i++) {
            (*new_tree)[i]->set_name((*new_exchange->get_dataset())[i].Sequence_name());
            (*new_tree)[i]->set_tip(TRUE);
            (*new_tree)[i]->set_taxa_id(i);
        }
        
        (*new_tree)[0]->initialized();
        (*new_tree)[1]->initialized();
        new_tree->set_as_siblings((*new_tree)[0], (*new_tree)[1]);
        new_parent=new_tree->get_nth_branch(new_exchange->get_num_taxa());
        new_parent->initialized();
        new_tree->set_as_parent_child(new_parent, (*new_tree)[0], 0);
        new_tree->set_as_parent_child(new_parent, (*new_tree)[1], 1);
        new_tree->set_root(new_parent);
        new_added_brns=new Branch_Ex *[3];
        
        new_added_brns[0]=new_tree->get_nth_branch(0);
        new_added_brns[1]=new_tree->get_nth_branch(1);
        new_added_brns[2]=new_parent;
        
        //Recursively assemble all possible trees
        assemble_trees(new_tree, new_exchange, 2);
        
        //Clean up
        delete new_tree;
        delete new_exchange;
        delete[] new_added_brns;
    }
    else if (num_added < partial_exchange->get_num_taxa())    {
        num_brns=2*num_added-1;
        
        //Add the new tip to every possible branch
        for(i=0; i<num_brns; i++) {
            
            new_exchange=new Exchange();
            (*new_exchange)=(*partial_exchange);
            new_tree=new Tree_Ex(new_exchange, the_matrix, new_exchange->is_rooted_tree());
            
            //Copy the tree as assembled so far into a new object
            (*new_tree)=(*partial_tree);
            
            new_added_brns=new Branch_Ex*[2*(num_added+1)-1];
            j=brncnt=0;
            while(j<new_exchange->get_num_branches()) {
                if ((*new_tree)[j]->is_uninitialized() == FALSE)
                    new_added_brns[brncnt++]=new_tree->get_nth_branch(j);
                j++;
            }
            if (brncnt != num_brns)
                cerr<<"ERROR: Did not reassemble array correctly\n";
            
            //For unrooted trees, we don't add to the null branch or the root
            if ((new_exchange->is_rooted_tree() == TRUE) ||
                ((new_added_brns[i] != new_tree->find_null_branch_id()) && (new_added_brns[i] != new_tree->find_root())) ) {
                
                new_brn=new_tree->get_nth_branch(num_added);
                new_brn->initialized();
                
                //The next "free" branch is after all n taxa + s-1 of the already added branches
                new_parent=new_tree->get_nth_branch(num_added-1+new_exchange->get_num_taxa());
                
                if (new_parent->is_uninitialized() == FALSE)
                    cerr<<"ERROR: Bad ordering of interior branches!\n";
                
                //cout<<"Adding branch "<<new_brn->get_brn_num()<<": "<<new_brn->get_name()
                //<<" to "<<new_added_brns[i]->get_brn_num()<<" : "<<new_added_brns[i]->get_name()<<endl;

                new_parent->initialized();
                
                if (new_added_brns[i]->get_sibling() != 0)
                    old_sibling=new_tree->get_nth_branch(new_added_brns[i]->get_sibling()->get_brn_num());
                else
                    old_sibling=0;
                if (new_added_brns[i]->get_parent() !=0)
                    old_parent=new_tree->get_nth_branch(new_added_brns[i]->get_parent()->get_brn_num());
                else
                    old_parent=0;
                
                //Could be adding at the root
                if (old_parent != 0) {
                    if (old_sibling==old_parent->get_child(0)) childnum=1;
                    else childnum=0;
                    
                    new_tree->set_as_parent_child(old_parent, new_parent, childnum);
                    new_tree->set_as_siblings(old_sibling, new_parent);
                }
                
                new_tree->set_as_parent_child(new_parent, new_brn, 0);
                new_tree->set_as_parent_child(new_parent, new_added_brns[i], 1);
                new_tree->set_as_siblings(new_brn, new_added_brns[i]);
                new_tree->set_root(new_parent);
                new_added_brns[num_brns]=new_brn;
                new_added_brns[num_brns+1]=new_parent;
                
                
                //Recursively assemble all possible trees
                assemble_trees(new_tree, new_exchange, num_added+1);
                
                //Clean up
                delete new_tree;
                delete new_exchange;
                delete[] new_added_brns;
            }
            else {
                delete new_tree;
                delete new_exchange;
                delete[] new_added_brns;
            }
        }
        
    }
    else {
        //We have a complete tree--optimize it, but first save the Tree_Ex object for copying later
       // for(i=0; i<partial_exchange->get_num_branches(); i++) {
         //   cout<<i<<" S: ";
         //   if ((*partial_tree)[i]->get_sibling() !=0) cout<<(*partial_tree)[i]->get_sibling()->get_brn_num();
         //   else cout<<"NONE";
         //   cout<<" P: ";
         //   if ( (*partial_tree)[i]->get_parent() !=0) cout<<(*partial_tree)[i]->get_parent()->get_brn_num();
         //       else cout<<"NONE";
          //  cout<<endl;
        //}
        
        full_tree=partial_tree;
        full_tree->diganose_tree((*full_tree)[0]);
        evaluate_tree(partial_tree, partial_exchange);
    }
    
}

void Exhaustive_Tree_PhyloMat_searcher::PhyloMat_allocate()
{
    int i;
    
    cout<<"Allocating searcher for PhyloMat trees"<<endl;
    
    best_params =new double* [save_num];
    best_trees_Ex=new Tree_Ex*[save_num];
   
    for(i=0; i<save_num; i++) {
        best_trees_Ex[i]=0;
        best_params[i]=new double [the_matrix->get_num_params()];
    }
    
}


void Exhaustive_Tree_PhyloMat_searcher::write_out_trees(int num, char *partial_file_name, const char *prog_name)
{
    int i, j;
    char filename[200], temp_num[40];
    Write_Tree_Arb_Model *writeout_tree;

    
    
    if (num>num_saved) {
        cerr<<"ERROR: "<<num<<" treefiles were requested when there are only "<<num_saved<<" in memory\n";
        num=num_saved;
    }
    for(i=0; i<num; i++) {
        //	for(j=0; j<curr_exchange->get_num_branches(); j++)
        //		(*best_trees[i])[j]->set_expect_subs_site((*best_trees[i])[j]->get_brnlen());
        
        writeout_tree=new Write_Tree_Arb_Model();
        strcpy(filename, partial_file_name);
        int_to_string (temp_num, 39, i);
        strcat(filename, temp_num);
        strcat(filename, ".tre");
        
        for(j=0; j<the_matrix->get_num_params(); j++)
            the_matrix->set_param(j, best_params[i][j]);
        
        writeout_tree->write_tree(filename, prog_name, the_matrix, 0, best_trees[i], best_exchanges[i]);
        delete writeout_tree;
    }
    
    cout<<"Wrote "<<num<<" treefiles\n";
}


void Exhaustive_Tree_PhyloMat_searcher::save_tree(int tree_num, Tree * source_tree)
{
    
    cout<<"Saving "<<tree_num<<" EXchange has lnL: "<<best_exchanges[tree_num]->get_saved_lnL()
        <<" Rooted val is "<<curr_exchange->is_rooted_tree()<<" total to save: "<<save_num<<endl;
    best_trees_Ex[tree_num]=new Tree_Ex(best_exchanges[tree_num], the_matrix, curr_exchange->is_rooted_tree());
    (*best_trees_Ex[tree_num])=(*full_tree);
    best_trees[tree_num]=best_trees_Ex[tree_num];
}

void Exhaustive_Tree_PhyloMat_searcher::setup()
{
    my_model = new Ploidy_Like_model(eval_exchange, eval_tree, the_genomes, the_WGX_homologs, the_matrix);
    my_model->the_tracks->update_tracking();
    eval_model=my_model;
}
