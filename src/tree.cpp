//Copyright 1999-2002 Gavin Conant

#include "tree.h"
#include <iostream>
#include <math.h>
#include <string.h>
#include <float.h>
#include <string>

using namespace::std;

//#define DEBUG


Branch::Branch()
{
  cerr<<"Call to default constructor of Branch class\n";
}

Branch::Branch(Exchange *cexchange)
{
	int i,j;
  
	curr_exchange=cexchange;
    brn_num=-1;
    
    param_set_id=0;
    
	strcpy(name, "NONE");
	tip=FALSE;
    zero_fixed=FALSE;
	taxa_id=-1;
	uninitialized=TRUE;
    extern_trpb=FALSE;
	basefreqs[0]=basefreqs[1]=basefreqs[2]=basefreqs[3]=0.25;
	brlen=0.0;
	expect_numsubs_site=0.0;
	parent=children[0]=children[1]=sibling=0;
	pns_num=0; 
	pitg_num=0;
	nonsyn_pattern=0;
	which_nonsyn_patt=&Branch::pns_num;
	has_name=FALSE;
    pruned=FALSE;
	dist_to_root=0;

	for(i=0; i<curr_exchange->get_num_aa_props()+1; i++)
		aa_prop_num[i]=0;

	
	if(curr_exchange->get_num_matrices() != 0)
		matrix_coeff_nums=new int[curr_exchange->get_num_matrices()];
	else
		matrix_coeff_nums=0;


	transprobs=new double**[curr_exchange->get_num_rates()];
	transprobs_prime=new double**[curr_exchange->get_num_rates()];
	transprobs_double_prime=new double**[curr_exchange->get_num_rates()];

	condprobs=new long double [curr_exchange->get_condlike_size()];
	
#ifdef _OPEN_MP_VERSION_
	condprobs_locale=new long double * [curr_exchange->get_num_open_mp_threads()];
	for(i=0; i<curr_exchange->get_num_open_mp_threads(); i++)
		condprobs_locale[i]=new long double [curr_exchange->get_condlike_size()];
#endif
	
	for(i=0; i<curr_exchange->get_num_rates(); i++)
	{
		transprobs[i]=new double*[curr_exchange->get_condlike_size()+1];
		transprobs_prime[i]=new double*[curr_exchange->get_condlike_size()];
		transprobs_double_prime[i]=new double*[curr_exchange->get_condlike_size()];

		  for(j=0; j<curr_exchange->get_condlike_size()+1; j++)
			transprobs[i][j]=new double [curr_exchange->get_condlike_size()+1];
		  for (j=0; j<curr_exchange->get_condlike_size(); j++) {
			  transprobs_prime[i][j]=new double[curr_exchange->get_condlike_size()];
			  transprobs_double_prime[i][j]=new double[curr_exchange->get_condlike_size()];
			}
    }


	//Hold derivatives of the transition probablities
	
	
}  //End Branch::Branch(Exchange *cexchange)

Branch::Branch(Exchange *cexchange, int b_num)
{
    int i,j;
    
    curr_exchange=cexchange;
    brn_num=b_num;
    
    strcpy(name, "NONE");
    tip=FALSE;
    zero_fixed=FALSE;
    taxa_id=-1;
    param_set_id=0;
    uninitialized=TRUE;
    extern_trpb=FALSE;
    basefreqs[0]=basefreqs[1]=basefreqs[2]=basefreqs[3]=0.25;
    brlen=0.0;
    expect_numsubs_site=0.0;
    parent=children[0]=children[1]=sibling=0;
    pns_num=0;
    pitg_num=0;
    nonsyn_pattern=0;
    which_nonsyn_patt=&Branch::pns_num;
    has_name=FALSE;
    pruned=FALSE;
    dist_to_root=0;
    
    for(i=0; i<curr_exchange->get_num_aa_props()+1; i++)
        aa_prop_num[i]=0;
    
    
    if(curr_exchange->get_num_matrices() != 0)
        matrix_coeff_nums=new int[curr_exchange->get_num_matrices()];
    else
        matrix_coeff_nums=0;
    
    
    transprobs=new double**[curr_exchange->get_num_rates()];
    transprobs_prime=new double**[curr_exchange->get_num_rates()];
    transprobs_double_prime=new double**[curr_exchange->get_num_rates()];
    
    condprobs=new long double [curr_exchange->get_condlike_size()];
    
#ifdef _OPEN_MP_VERSION_
    condprobs_locale=new long double * [curr_exchange->get_num_open_mp_threads()];
    for(i=0; i<curr_exchange->get_num_open_mp_threads(); i++)
        condprobs_locale[i]=new long double [curr_exchange->get_condlike_size()];
#endif
    
    for(i=0; i<curr_exchange->get_num_rates(); i++)
    {
        transprobs[i]=new double*[curr_exchange->get_condlike_size()+1];
        transprobs_prime[i]=new double*[curr_exchange->get_condlike_size()];
        transprobs_double_prime[i]=new double*[curr_exchange->get_condlike_size()];
        
        for(j=0; j<curr_exchange->get_condlike_size()+1; j++)
            transprobs[i][j]=new double [curr_exchange->get_condlike_size()+1];
        for (j=0; j<curr_exchange->get_condlike_size(); j++) {
            transprobs_prime[i][j]=new double[curr_exchange->get_condlike_size()];
            transprobs_double_prime[i][j]=new double[curr_exchange->get_condlike_size()];
        }
    }
    
    
    //Hold derivatives of the transition probablities
    
    
}  //End Branch::Branch(Exchange *cexchange)


Branch  & Branch::operator=(Branch & assign_from)
{
  int i, j, k;
 
  set_name(assign_from.get_name());
 
  brlen=assign_from.get_brnlen();
  expect_numsubs_site=assign_from.expect_subs_site();
	has_name=assign_from.branch_has_name();

  uninitialized=assign_from.is_uninitialized();
    param_set_id=assign_from.get_param_set();

  tip=assign_from.is_tip();
    zero_fixed=assign_from.is_fixed_zero();
  if (tip == TRUE) 
	  children[0]=children[1]=0;
  
    pruned=assign_from.is_pruned();
  taxa_id=assign_from.get_taxa_id();
  pns_num=assign_from.get_p_nonsyn_num();
  pitg_num=assign_from.get_p_intergroup_num();
  for(i=0; i<curr_exchange->get_num_aa_props()+1; i++)
	  aa_prop_num[i]=assign_from.get_aa_prop_num(i);

    if (assign_from.has_extern_trpb() == FALSE) {
        for(i=0; i<curr_exchange->get_num_rates(); i++)
            for(j=0; j<curr_exchange->get_condlike_size()+1; j++)
                for(k=0; k<curr_exchange->get_condlike_size()+1; k++)
                    transprobs[i][j][k]=assign_from.get_trpb(i,j,k);
    }
  //Don' think this is right--MS VC complains w/o it
  return(*this);
}

double Branch::get_brlen_ci(BOOL upper)
{
	if (upper ==TRUE)
		return(brlen_upper);
	else
		return(brlen_lower);
}


void Branch::set_brlen_ci(BOOL upper, double val)
{
	if (upper == TRUE)
		brlen_upper=val;
	else
		brlen_lower=val;
}

void Branch::set_aa_prop_num(int aanum, AA_PROPERTIES prop)
{
	aa_prop_num[curr_exchange->get_prop_index_num(prop)]=aanum;
}

 
void Branch::set_matrix_coeff_num(int matrix, int matrix_num)
{
	if (matrix_coeff_nums != 0)
		matrix_coeff_nums[matrix]=matrix_num;
	else
		cerr<<"Error: request to set matrix coeff when matrices are not in use\n";

}

void Branch::set_nonsyn_pattern(int pat)
{
	which_nonsyn_patt=&Branch::nonsyn_pattern;
	nonsyn_pattern=pat;

}


double Branch::get_trpb(int rate, int start, int end)
{
	if (end<curr_exchange->get_condlike_size()) {
		if (start<curr_exchange->get_condlike_size())
			return(transprobs[rate][start][end]);
		else	
			return(transprobs[rate][curr_exchange->get_condlike_size()][end]);
	}
	else
		return(transprobs[rate][start][curr_exchange->get_condlike_size()]);
} //End Branch::get_trpb


int Branch::get_aa_prop_num(AA_PROPERTIES prop)
{
	return(aa_prop_num[curr_exchange->get_prop_index_num(prop)]);
}


int Branch::get_aa_prop_num(int prop)
{
	return(aa_prop_num[prop]);
}



int Branch::get_matrix_coeff_num(int matrix)
{
	if (matrix_coeff_nums != 0)
		return(matrix_coeff_nums[matrix]);
    else {
		cerr<<"Error: request to get matrix coeff when matrices are not in use\n";
        return(-1);
    }
}



Branch * Branch::get_child(int childnum)
{
  if (childnum==0)
    return(children[0]);
  else
    return(children[1]);
}




void Branch::set_child(Branch *newchild, int childnum)
{
  if (childnum==0)
    children[0]=newchild;
  else
    children[1]=newchild;
}


void Branch::null_child(int childnum)
{
  if (childnum==0)
    children[0]=0;
  else
    children[1]=0;
}


void Branch::name_this_branch()
{
	int new_len;
	if (tip == TRUE) has_name=TRUE;
	else if (strcmp(name, "ROOT") == 0) {
		has_name=TRUE;
		children[0]->name_this_branch();
		children[1]->name_this_branch();
		//cout<<"Root is named\n";
	}
	else {
		if (children[0]->branch_has_name() == FALSE) children[0]->name_this_branch();
		if (children[1]->branch_has_name() == FALSE) children[1]->name_this_branch();
		strcpy(name, "\0");
		strcat(name, children[0]->get_name());
		strcat(name, "+");
		strcat(name, children[1]->get_name());
		//cout<<"Set name to "<<name<<endl;
		has_name=TRUE;
	}
}


void Branch::set_dist_to_root()
//SHOULD ONLY BE CALLED BY TREE OBJECT FROM ROOT
{
	if (parent !=0)
		dist_to_root=parent->get_dist_to_root()+expect_numsubs_site;
	else
		dist_to_root=expect_numsubs_site;
	
	//cout<<"Name: "<<name<<" Dist: "<<dist_to_root<<endl;
	
	if (tip == FALSE) {
		children[0]->set_dist_to_root();
		children[1]->set_dist_to_root();
	}
		
}

void Branch::set_param_set(int n)
{
    cout<<"Assiging set it "<<n<<" to  brnach num "<<brn_num<<endl;
    param_set_id=n;
}

Branch::~Branch()
{
  int i, j;

  
  for(i=0; i<curr_exchange->get_num_rates(); i++)
    {
      for(j=0; j<curr_exchange->get_condlike_size()+1; j++)
		delete[] transprobs[i][j];
	  for(j=0; j<curr_exchange->get_condlike_size(); j++) {
			delete[] transprobs_prime[i][j];
			delete[] transprobs_double_prime[i][j];
		}

      delete[] transprobs[i];
	  delete[] transprobs_prime[i];
	  delete[] transprobs_double_prime[i];
    }
  
 
  


  delete[] transprobs;
  delete[] condprobs;

#ifdef _OPEN_MP_VERSION_
	for (i=0; i<curr_exchange->get_num_open_mp_threads(); i++) delete[] condprobs_locale[i];
	delete[] condprobs_locale;
#endif
	
  delete[] transprobs_prime;
  delete[] transprobs_double_prime;

  if(matrix_coeff_nums != 0)
	  delete[] matrix_coeff_nums;

}





Tree::Tree()
{
  //cerr<<"Call to default constructor of Tree class\n";
}



Tree::Tree(Exchange *cexchange)
{
  int i;
  curr_exchange=cexchange;

  //We flag that this tree contains its own local data
  tree_is_local_mem=TRUE;
   null_root_brlns=TRUE;
   is_rooted_tree=FALSE;

	start_extra_brn=new Branch_list;
	extra_brn=start_extra_brn;
	start_tips=new Branch_list;
	tips=start_tips;

	start_tips->next=0;
	start_extra_brn->next=0;

  if (curr_exchange->get_num_taxa() <= 3)
    three_taxa_tree=TRUE;
  else
    three_taxa_tree=FALSE;

  prune_root=0;
  
  //Declares a new array of branches large enough for the tree
  tree = new Branch * [curr_exchange->get_num_branches()];
  
  for (i=0; i<curr_exchange->get_num_branches(); i++)
    tree[i]=new Branch(curr_exchange);
}

Tree::Tree(Exchange *cexchange, BOOL rooted)
{
  int i;
  curr_exchange=cexchange;

  //We flag that this tree contains its own local data
  tree_is_local_mem=TRUE;
   is_rooted_tree=rooted;
   if(is_rooted_tree==TRUE)
	   null_root_brlns=FALSE;
   else
	   null_root_brlns=TRUE;

	start_extra_brn=new Branch_list;
	extra_brn=start_extra_brn;
	start_tips=new Branch_list;
	tips=start_tips;

	start_tips->next=0;
	start_extra_brn->next=0;

  if (curr_exchange->get_num_taxa() <= 3)
    three_taxa_tree=TRUE;
  else
    three_taxa_tree=FALSE;

  prune_root=0;
	constrain_brn=0;
  
  //Declares a new array of branches large enough for the tree
  tree = new Branch * [curr_exchange->get_num_branches()];
  
  for (i=0; i<curr_exchange->get_num_branches(); i++)
    tree[i]=new Branch(curr_exchange, i);
}




Tree::Tree(Exchange *cexchange, Branch *new_root)
{
   int i;
   curr_exchange=cexchange;

   //This tree uses memory allocated in another tree object--
   //We won't dealloc it in the destructor
   tree_is_local_mem=FALSE;
	null_root_brlns= FALSE;
	start_tips=0;
	start_extra_brn=0;

   prune_root=0;
   

   root=new_root;
   if (new_root->is_tip() == FALSE) {
     if (!((root->get_child(0)->is_tip() ==TRUE) && (root->get_child(1)->is_tip() ==TRUE))) {
       set_root(root);
       if ((find_null_branch()->get_child(0)->is_tip() ==TRUE) &&
	   (find_null_branch()->get_child(0)->is_tip() ==TRUE) && find_null_branch()->get_sibling()->is_tip()==TRUE)
	 three_taxa_tree=TRUE;
       else
	 three_taxa_tree=FALSE;
     }
   }
  
}



Branch * Tree::operator[] (int element)
{
  if (element<curr_exchange->get_num_branches())
    return(tree[element]);
  else
    {      cerr<<"Invalid Branch number\n";
      return(dummy);
    }  
}



Tree & Tree::operator= (Tree & assign_from)
{
  int i, j;

  for (i=0; i<curr_exchange->get_num_branches(); i++)
  {
	
    *(tree[i])=*(assign_from[i]);
    
    if (assign_from[i]!=assign_from.find_root())
    {   
		//Added a bunch of checks so that we can copy trees that are partially assembled:
		//i.e. where some branches do not have parents/siblings/children yet
		j=0;
      
		if(assign_from[i]->get_sibling() != 0 ){
			while (assign_from[j]!=assign_from[i]->get_sibling()) j++;
			set_as_siblings(tree[i], tree[j]);
		}

		if (tree[i]->is_tip()==FALSE)
		{	
			j=0;
			
			if(assign_from[i]->get_child(0) != 0) {
				while (assign_from[j]!=assign_from[i]->get_child(0)) j++;
				set_as_parent_child(tree[i], tree[j], 0);
			}

			j=0;
			if(assign_from[i]->get_child(1) != 0) {
				while (assign_from[j]!=assign_from[i]->get_child(1)) j++;
				set_as_parent_child(tree[i], tree[j], 1);
			}
		}
		j=0;
      
		if (assign_from[i]->get_parent() != 0) {
			while (assign_from[j]!=assign_from[i]->get_parent()) j++;
			if (assign_from[i]->get_parent()->get_child(0)==assign_from[i])
				set_as_parent_child(tree[j], tree[i], 0);
			else
				set_as_parent_child(tree[j], tree[i], 1);
		}
    }
   
    else
    {
      j=0;
	
	  if (assign_from[i]->get_child(0) != 0) {
			while (assign_from[j]!=assign_from[i]->get_child(0)) j++;
			set_as_parent_child(tree[i], tree[j], 0);
	  }
      j=0;
      
	  if (assign_from[i]->get_child(1) != 0) {
			while (assign_from[j]!=assign_from[i]->get_child(1)) j++;
			set_as_parent_child(tree[i], tree[j], 1);
	  }
      tree[i]->null_parent();
      tree[i]->null_sibling();
    }
  }
  set_root(tree[0]);
  //Don't think this is right--MS VC complains w/o it
  return(assign_from);
}



void Tree::describe_tree()
{
  //Prints out the siblings of all of the tip nodes in the tree array structure

  for (int i=0; i<curr_exchange->get_num_branches(); i++)
    if (tree[i]!=root)
      cout<<tree[i]->get_name()<<" "<<" (Sibling: "<<tree[i]->get_sibling()->get_name()<<"): " 
	 <<" Len: "<<tree[i]->get_brnlen()<<endl;
    else
      cout<<"ROOT: "<<tree[i]->get_name()<<" Len: "<<tree[i]->get_brnlen()<<endl;
  
}  //End Tree::describe_tree



Branch* Tree::initialize_branch()
{
	int start_pos;
	Branch * start;
	
	start_pos=curr_exchange->get_num_taxa();
	
	start=tree[start_pos];
	
	while(start->is_uninitialized()==FALSE && start_pos<curr_exchange->get_num_branches())
		start=tree[++start_pos];
	
	if (start_pos==curr_exchange->get_num_branches())
		cerr<<"No empty branch avaliable for request\n";
		start->initialized();
		return(start);
}  //End Tree::initialize_branch




Branch * Tree::get_pott_branch_num(Branch *start, int &curr_num, int target_num)
//This function finds the branch with post-order tree transversal number ==
//target_num + (original curr_num +1).  Thus, if called from the tree root with
//curr_num == -1, it will find the target_num +1 th branch in the post-order transversal
{
  Branch *ret=0;

 
  //First progress to the left tip
  if (start->is_tip() == FALSE)
    ret=get_pott_branch_num(start->get_child(0), curr_num, target_num);
 
  //Have we already found the branch--then return
  if (ret != 0)
    return(ret);
  
  //Number this branch
  curr_num++;
 
  //This branch is the target
  if (curr_num == target_num)
    return(start);

  //Check if this branch is a left child (root check performed to avoid
  //seg fault if check parent of root)
  if ((start != root) && (start->get_parent()->get_child(0)==start))
    {
      //If sib is tip, number it and check if it is the target
      if (start->get_sibling()->is_tip() == TRUE)
	{
	  curr_num++;
	   if (curr_num == target_num)
	     return(start->get_sibling());
	   else
	     return(0);
	}
      //If our sib is not a tip, number that subtree
      else
	return(get_pott_branch_num(start->get_sibling(), curr_num, target_num));
    }
  //If we are at a right child just return
  else
    return(ret);
}



int Tree::count_branches_above(Branch *start)
{
  int retval=0;
  if (start->is_tip() == FALSE)
    return(1+count_branches_above(start->get_child(0))+
      count_branches_above(start->get_child(1)));
  else
    return(1);
}


Branch * Tree::tree_bisect_reconnect(Branch * split_loc, int tree_1_join, int tree_2_join, int base_tree)
  //Does a tree bisection and reconnection operation.  Split occurs at split_loc.  Post-oder tree numbering is used to identify the two branches numbered tree_1_join and tree_2_join.  Tree_1_join is the "smaller" of the trees 
  //(i.e. the tree "above" split loc)
{
  int join1num=-1, join2num=-1;
  BOOL old_null_val;
  Branch *join1, *join2, *retval=0;
  Tree *newtree;

  old_null_val=null_root_brlns;
  null_root_brlns=FALSE;
  newtree=split_tree(split_loc);

  if (newtree != 0) {
   
    join1 = newtree->get_pott_branch_num(newtree->find_root(), join1num, tree_1_join);
    join2 = get_pott_branch_num(root, join2num, tree_2_join);
        
    if ((join1 == 0) || (join2 == 0)) 
      cerr<<"Tree bisection and reconnection failed: invalid join point on one or both subtrees\n";
    else
      retval=join_trees(newtree, this, join1, join2, base_tree);
    
    delete newtree;
  }
  null_root_brlns=old_null_val;
  if(retval != 0)
    set_root(retval);
  return(retval);
}




void Tree::re_root_tree(Branch * root_from)
{
  int root_from_child_num;


  if (root_from != root) {
      if (root->is_tip() != TRUE) {
          if (!((root->get_child(0)->is_tip() ==TRUE) && (root->get_child(1)->is_tip() ==TRUE))) {
              if (three_taxa_tree == TRUE) {
                  if (root_from != find_null_branch()->get_sibling()) {
                      old_sibling=root_from->get_sibling();
                      old_parent=root_from->get_parent();
                      single_branch=root->get_child(1);
            
                      if (old_parent->get_child(0)==root_from)
                          root_from_child_num=0;
                      else
                          root_from_child_num=1;
            
                      set_as_parent_child(root, root_from, 1);
                      set_as_siblings(root_from, find_null_branch());
            
                      set_as_parent_child(old_parent, single_branch, root_from_child_num);
                      set_as_siblings(single_branch, old_sibling);
                  }
              }
              else {
                  if (find_null_branch_id()->get_sibling()!=root_from) {
                      single_branch=root->get_child(1);
                      pair_branch1=find_null_branch()->get_child(0);
                      pair_branch2=find_null_branch()->get_child(1);
                      old_sibling=root_from->get_sibling();
                      old_parent=root_from->get_parent();
              
              
                      find_null_branch()->get_child(0)->null_parent();
                      find_null_branch()->get_child(1)->null_parent();
                      find_null_branch()->get_sibling()->null_parent();
              
                      set_as_parent_child(root, root_from, 1);
                      set_as_siblings(root_from, find_null_branch());
              
            
                      if (old_parent == find_null_branch()) {
                          set_as_parent_child(find_null_branch(), old_sibling, 0);
                          set_as_parent_child (find_null_branch(), single_branch, 1);
                          set_as_siblings(old_sibling, single_branch);
                      }
                      else {
                          current=invert_clade(find_null_branch(), old_parent, old_sibling);
              
                          if (current==single_branch) {
                              set_as_parent_child(current, pair_branch1, 0);
                              set_as_parent_child(current, pair_branch2, 1);
                          }
                          else {
                              if (current==pair_branch1)
                                  set_as_parent_child(current, pair_branch2, 0);
                              else
                                  set_as_parent_child(current, pair_branch1, 0);
                              set_as_parent_child(current, single_branch, 1);
                          }
              
                          set_as_siblings(current->get_child(0), current->get_child(1));
                      }
                  }
                  set_root(root_from);
              }
              diganose_tree(root);
          }
      }
  }
}


void Tree::set_root(Branch *start)
  //This function determines the root of the tree, given some branch on the tree (start)
{
  Branch *loc, *new_null;
 
    if (start == 0) loc = tree[0];
    else
        loc=start;

    
 

  while(loc->get_parent()!=0)
    loc=loc->get_parent();


  root=loc; 
	//cout<<"Root 0: "<<loc->get_child(0)->get_taxa_id()<<": "<<loc->get_child(0)->get_name()<<endl;
	//cout<<"Root 1: "<<loc->get_child(1)->get_taxa_id()<<": "<<loc->get_child(1)->get_name()<<endl;
	
 
  if (root->is_tip() == FALSE) {

    if (!((root->get_child(0)->is_tip() == TRUE) && (root->get_child(1)->is_tip() == TRUE)))
      {
	if (find_null_branch_id()->is_tip()==TRUE)
	  {
	    new_null=find_null_branch_id()->get_sibling();
	    set_as_parent_child(root, find_null_branch_id(), 1);
	    set_as_parent_child(root, new_null, 0);
	  }
	
	if (null_root_brlns == TRUE) {
	 
	  find_null_branch_id()->get_sibling()->set_brnlen(find_null_branch_id()->get_sibling()->get_brnlen()+
							   find_null_branch_id()->get_brnlen());
	  
	  find_null_branch_id()->get_sibling()->set_expect_subs_site(find_null_branch_id()->get_sibling()->expect_subs_site()+
								     find_null_branch_id()->expect_subs_site());
	  
	  find_null_branch_id()->set_brnlen(0.0);
	  find_null_branch_id()->set_expect_subs_site(0.0);      

	 
	}	
      }
    //Moves root branch length to one of it's children
    if (((root->get_brnlen() != 0.0) || (root->expect_subs_site() != 0.0)) && (null_root_brlns == TRUE)) {
      find_null_branch_id()->get_sibling()->set_brnlen(find_null_branch_id()->get_sibling()->get_brnlen()+
						       root->get_brnlen());
      find_null_branch_id()->get_sibling()->set_expect_subs_site(find_null_branch_id()->get_sibling()->expect_subs_site()+
								 root->expect_subs_site());
      
      root->set_brnlen(0.0);
      root->set_expect_subs_site(0.0);
      
    }
    while (loc->is_tip()==FALSE)
      loc=loc->get_child(0);
    
    left_tip=loc;
      
      //cout<<"Left tip name: "<<left_tip->get_name()<<endl;
  
  }

}




void Tree::set_rooted_tree(BOOL rooted)
{
  is_rooted_tree=rooted;
}



Branch * Tree::get_leftmost_tip()
{
  return(left_tip);
}




Branch * Tree::find_left_tip(Branch * start)
{
  while (start->is_tip()==FALSE)
    start=start->get_child(0);
  return(start);
}


void Tree::prune_and_reconnect(Branch *subtree_root, Branch *reconnect_loc, double split_fac)
{
	int i, child_num;	
	Branch *old_par, *old_sib, *left, *right, *reconnect_sib, *reconnect_par, *orig_null=0;
	extra_brn=start_extra_brn;
	tips=start_tips;
	diganose_tree(root);
	get_tips(subtree_root, subtree_root);
	tips=start_tips;
	while (tips->element != 0) {
		i=0;
		while (tree[i]->get_taxa_id() != tips->element->get_taxa_id()) i++;

		if (tree[i]->get_parent()->get_parent() != 0) {
			if ((reconnect_loc==tree[i]) || (reconnect_loc==tree[i]->get_parent()))
				reconnect_loc=tree[i]->get_parent()->get_parent();
			add_to_branch_list(extra_brn, tree[i]);
			add_to_branch_list(extra_brn, tree[i]->get_parent());
			old_par=tree[i]->get_parent()->get_parent();
			old_sib=tree[i]->get_sibling();

			if (old_par->get_child(0)==tree[i]->get_parent())
				child_num=0;
			else
				child_num=1;
			set_as_parent_child(old_par, old_sib,child_num);
			set_as_siblings(old_par->get_child(0), old_par->get_child(1));
		}
		else {
			if ((reconnect_loc==tree[i]) || (reconnect_loc == tree[i]->get_parent()))
				reconnect_loc=tree[i]->get_sibling();
			add_to_branch_list(extra_brn, tree[i]);
			add_to_branch_list(extra_brn, tree[i]->get_parent());
			
			old_sib=tree[i]->get_sibling();
			old_sib->null_parent();
			old_sib->null_sibling();
			root=old_sib;
		}
		tips=tips->next;
	}


	
	extra_brn=start_extra_brn;
	build_tree(subtree_root, subtree_root, left, right);
	old_sib=left;
	if (old_sib->get_brnlen() == 0)
		old_sib->set_brnlen(subtree_root->get_sibling()->get_brnlen()*split_fac);
	old_par=extra_brn->element;
	old_par->set_tip(FALSE);
	old_par->set_taxa_id(-1);
	
	if (reconnect_loc->get_brnlen() == 0) {
		reconnect_loc->set_brnlen(old_par->get_brnlen());
		old_par->set_brnlen(0);
		if (reconnect_loc->get_parent() == 0)
			if (reconnect_loc->get_child(0)->get_brnlen() == 0)
				orig_null=reconnect_loc->get_child(0);
	}
	
	

	reconnect_sib=reconnect_loc->get_sibling();
	reconnect_par=reconnect_loc->get_parent();
	if (reconnect_par != 0) {
		if (reconnect_par->get_child(0) == reconnect_loc)
			child_num=0;
		else
			child_num=1;
	}


	set_as_siblings(reconnect_loc, old_sib);
	reconnect_loc->set_brnlen(reconnect_loc->get_brnlen()*split_fac);
	old_par->set_brnlen(reconnect_loc->get_brnlen()*(1.0-split_fac));
	set_as_parent_child(old_par, reconnect_loc, 0);
	set_as_parent_child(old_par, old_sib, 1);
	if (orig_null != 0)
	{
		orig_null->set_brnlen(old_par->get_child(0)->get_brnlen());
		old_par->get_child(0)->set_brnlen(0);
	}
	
	if (reconnect_sib != 0)
		set_as_siblings(reconnect_sib, old_par);
	if (reconnect_par != 0)
		set_as_parent_child(reconnect_par, old_par, child_num);
	else {
		old_par->null_sibling();
		old_par->null_parent();
	}

	set_root(reconnect_loc);
	diganose_tree(root);
	diganose_tree(root);

}


void Tree::prune_from_tree(Branch *prune_branch)
{
  int sib_num;
  Branch *new_parent;
 
  if (prune_root==0) {
      if (prune_branch==find_root() || prune_branch==find_null_branch_id())
          cerr<<"Cannot prune the root or null branch from the tree\n";
      else if (prune_branch!=find_null_branch_id()->get_sibling()) {
          prune_root=prune_branch->get_parent();
          
          new_parent=prune_root->get_parent();

          if (new_parent->get_child(0)==prune_root)
            sib_num=0;
          else
            sib_num=1;
          
          set_as_parent_child(new_parent, prune_branch->get_sibling(), sib_num);
          set_as_siblings(new_parent->get_child(0), new_parent->get_child(1));
          prune_branch->null_sibling();
          if(prune_root->get_child(0)==prune_branch)
            prune_root->null_child(1);
          else
            prune_root->null_child(0);
      }

      else {
	  
          prune_root=find_null_branch_id();
          set_as_parent_child(find_root(), find_null_branch_id()->get_child(0), 0);
          set_as_parent_child(find_root(), find_root()->get_child(0)->get_sibling(), 1);
          

          set_as_parent_child(prune_root, prune_branch, 0);
          prune_root->null_child(1);
          prune_root->null_parent();
          prune_branch->null_sibling();
	}
     
  }
  else
    cerr<<"An existing part of the tree has already been pruned--reattach before\ncontinuing\n";
}  //End Tree::prune_from_tree





void Tree::add_back_prune(Branch *add_to)
{  int child_num, prune_child;

  if (prune_root==0)
    cerr<<"Must prune before adding pruned piece back to tree\n";
  else
    {
      if (add_to==find_root())
	cerr<<"Can't add pruned subtree to the root of current tree\n";
      else
	{
	  if (add_to->get_parent()->get_child(0)==add_to)
	    child_num=0;
	  else
	    child_num=1;

	  set_as_parent_child(add_to->get_parent(), prune_root, child_num);
	  set_as_siblings(add_to->get_parent()->get_child(0), add_to->get_parent()->get_child(1));

	  if (prune_root->get_child(0)==0)
	    prune_child=0;
	  else
	    prune_child=1;

	  set_as_parent_child(prune_root, add_to, prune_child);
	  set_as_siblings(prune_root->get_child(0), prune_root->get_child(1));
	  prune_root=0;
	  set_root(tree[0]);
	}
    }

}  //End Tree::add_back_prune



void Tree::remove_tip(Branch * taxa)
{
	int i, sib_num, old_id;
	Branch *old_parent, *old_sib, *new_parent, *child1, *child2;
	
	if(taxa->is_tip() == TRUE) {	
		old_id=taxa->get_taxa_id();

		if (taxa==find_null_branch_id()->get_sibling()) {
			root=taxa->get_sibling();
			
			root->null_parent();
			root->null_sibling();
			
			if(root->get_child(0)->is_tip() == TRUE) {
				child1=root->get_child(0);
				child2=root->get_child(1);
				root->set_child(child2, 0);
				root->set_child(child1, 1);
			}
			
			root->get_child(1)->set_brnlen(root->get_child(1)->get_brnlen() + root->get_child(0)->get_brnlen());
			root->get_child(0)->set_brnlen(0.0);
		}
		else {
			old_sib=taxa->get_sibling();
			old_parent=taxa->get_parent();
			new_parent=old_parent->get_parent();
			if (new_parent->get_child(0) == old_parent) sib_num=0;
			else sib_num=1;
			
			set_as_parent_child(new_parent, old_sib, sib_num);
			set_as_siblings(new_parent->get_child(0), new_parent->get_child(1));

            old_parent->set_pruned();
            
			if (taxa->get_parent() == find_null_branch_id()) {
				if (old_sib->is_tip() == TRUE) {
					old_sib->set_brnlen(old_sib->get_brnlen() + old_sib->get_sibling()->get_brnlen());
					old_sib->get_sibling()->set_brnlen(0.0);
					
					if (old_sib == find_root()->get_child(0)) {
						find_root()->set_child(old_sib, 1);
						find_root()->set_child(old_sib->get_sibling(), 0);
					}
					
				}
				else {
					old_sib->get_sibling()->set_brnlen(old_sib->get_sibling()->get_brnlen() + old_sib->get_brnlen());
					old_sib->set_brnlen(0.0);
				}
			}
		}
		
		for(i=0; i<curr_exchange->get_num_taxa(); i++)
			if (tree[i]->get_taxa_id() > old_id) tree[i]->set_taxa_id(tree[i]->get_taxa_id() -1);
		taxa->set_taxa_id(curr_exchange->get_num_taxa()-1);
        taxa->set_pruned();
		
	}
	else
		cerr<<"ERROR: Cannot remove non-tip branches\n";
	
}


void Tree::copy_tree(Tree *intree) 
{
	int i, j;
	Branch *new_branch, *new_root;
	
	for (i=0; i<curr_exchange->get_num_taxa(); i++)
	{
		tree[i]->set_tip(TRUE);
		tree[i]->set_taxa_id(i);
		tree[i]->initialized();
	}
	
	new_root=initialize_branch();
	
	ascend_assign_clade (intree->find_root(),new_root);
	
	set_root(tree[0]);
}


void Tree::diganose_tree(Branch *start)
{
	BOOL error=FALSE; 
  if (start->get_sibling()==0 && start!=root)
    cout<<"Taxa "<<start->get_taxa_id()<<" sibling is null\n"<<flush;
  if (start->get_parent()==0 && start!=root)
    cout<<"Taxa "<<start->get_taxa_id()<<" parent is null\n"<<flush;
  if (start->get_child(0)==0 && start->is_tip()==FALSE) 
    cout<<"Taxa "<<start->get_taxa_id()<<" child 0 is null\n"<<flush;
  if (start->get_child(1)==0 && start->is_tip()==FALSE)
    cout<<"Taxa "<<start->get_taxa_id()<<" child 1 is null\n"<<flush;
  if (((start->get_child(0) != 0) || (start->get_child(1) != 0)) && (start->is_tip() ==TRUE)) {
    cout<<"Taxa "<<start->get_taxa_id()<<" is tip with child\n"<<flush;
	error=TRUE;
  }
  if((start->is_tip() == TRUE)  && (strcmp(start->get_name(), "NONE") == 0))
  {
	cout<<"Tip taxa "<<start->get_taxa_id()<<" has invalid name\n"<<flush;
	error=TRUE;
  }
	//if (((start != root) && (start != find_null_branch()) && start->get_brnlen() == 0))
	//	cerr<<"Zero-length branch "<<start<<" is not root or null\n";
  
  if ((start->is_tip()==FALSE) && (error == FALSE))
    {
      diganose_tree(start->get_child(0));
      diganose_tree(start->get_child(1));
    }
  else if (error == TRUE)
	cout<<"Exiting\n";
}


int Tree::find_max_tree_depth()
{
	int i, retval=0, cnt;
	Branch *start_brn;

	for(i=0; i<curr_exchange->get_num_branches(); i++) {
		if (tree[i]->is_tip() == TRUE) {
			cnt=1;
			start_brn=tree[i];
			while (start_brn != root) {
				start_brn=start_brn->get_parent();
				cnt++;
			}
			if (cnt > retval) retval=cnt;
		}
	}
	return(retval);
}


BOOL Tree::is_root(Branch *the_brn)
{
  if (the_brn->get_parent()==0)
    return(TRUE);
  else
    return(FALSE);
}




Tree::~Tree()
{
	int i;
	Branch_list *temp;
	
	while (start_extra_brn != 0)
	{
		temp=start_extra_brn->next;
		delete start_extra_brn;
		start_extra_brn=temp;
	}

	while (start_tips != 0)
	{
		temp=start_tips->next;
		delete start_tips;
		start_tips=temp;
	}

  if (tree_is_local_mem==TRUE) {
    for (i=0; i<curr_exchange->get_num_branches(); i++)
      delete tree[i];
    
    delete[] tree;
  }
}
//End of Public Functions





void Tree::set_as_siblings(Branch *sib1, Branch *sib2)
{
  sib1->set_sibling(sib2);
  sib2->set_sibling(sib1);
}  //end Tree::set_as_siblings




void Tree::set_as_parent_child(Branch *parent, Branch *child, int child_num)
{
  parent->set_child(child, child_num);
  child->set_parent(parent);
} //End Tree::set_as_parent_child


int Tree::other_child(int child_num)
{
  if (child_num==0)
    return(1);
  else
    return(0);
}


//Private Functions:
void Tree::get_tips(Branch *start, Branch *myroot)
{

	//First progress to the left tip
	if (start->is_tip() == FALSE)
		get_tips(start->get_child(0), myroot);
	else 
		add_to_branch_list(tips, start);	
	
 
	//Check if this branch is a left child (root check performed to avoid
	//seg fault if check parent of root)
	if ((start != myroot) && (start->get_parent()->get_child(0)==start))
    {
	  //If sib is tip, number it and check if it is the target
      if (start->get_sibling()->is_tip() == TRUE)
		add_to_branch_list(tips, start->get_sibling());	
      //If our sib is not a tip, number that subtree
      else {
		get_tips(start->get_sibling(), myroot);
		return;
	  }
    }
  //If we are at a right child--just return
  else
    return;
}



void Tree::build_tree(Branch *tree_loc, Branch *myroot, Branch *&child1, Branch *&child2)
{
	Branch *new_brn, *new_sib;

	//First progress to the left tip
	if (tree_loc->is_tip() == FALSE) {
		build_tree(tree_loc->get_child(0), myroot, child1, child2);
		new_brn=extra_brn->element;
		extra_brn=extra_brn->next;
		(*new_brn)=(*tree_loc);
		set_as_parent_child(new_brn, child1, 0);
		set_as_parent_child(new_brn, child2, 1);
	}	
	else {
		new_brn=extra_brn->element;
		
		extra_brn=extra_brn->next;
		(*new_brn)=(*tree_loc);	
	}
 
	//Check if this branch is a left child (root check performed to avoid
	//seg fault if check parent of root)
	if ((tree_loc != myroot) && (tree_loc->get_parent()->get_child(0)==tree_loc))
    {
		if (tree_loc->get_sibling()->is_tip() == TRUE)
		{
			new_sib=extra_brn->element;
			extra_brn=extra_brn->next;
			(*new_sib)=(*tree_loc->get_sibling());		  
		}
		//If our sib is not a tip, traverse that subtree
		else {
			build_tree(tree_loc->get_sibling()->get_child(0), myroot, child1, child2);
			new_sib=extra_brn->element;
			extra_brn=extra_brn->next;
			(*new_sib)=(*tree_loc->get_sibling());
			set_as_parent_child(new_sib, child1, 0);
			set_as_parent_child(new_sib, child2, 1);
		}
		set_as_siblings(new_brn, new_sib);
		child1=new_brn;
		child2=new_sib;
		return;
    }
	else if (tree_loc == myroot) {
		child1=new_brn;
		return;
	}
  //If we are at a right child--just return
	else
		return;	
}



Branch* Tree::invert_clade(Branch *new_parent, Branch *old_parent, Branch *old_sibling)
{
     Branch *old_par, *old_sib;

     old_par=old_parent->get_parent();
     old_sib=old_parent->get_sibling();

     set_as_parent_child(new_parent, old_parent, 0);
     set_as_parent_child(new_parent, old_sibling, 1);
     set_as_siblings(old_parent, old_sibling);


     if (old_par!=0)
         return(invert_clade(old_parent, old_par, old_sib));
     else
         return(old_parent);

} //End Tree::invert_clade





Tree * Tree::split_tree(Branch *split_point)
  //This function takes a rooted tree, splits it at branch 
  //split_point and returns a pointer to one of the two new trees
  //(The other remains in the original tree pointer
{
  int old_par_child_num;
  Tree *new_tree;
  Branch *old_par, *old_sib, *old_par_sib;

  if (split_point != root) {
    old_par=split_point->get_parent();
    old_sib=split_point->get_sibling();
  
    if (old_par != root) { 
      old_par_sib=old_par->get_sibling();
      if (old_par->get_parent()->get_child(0) == old_par)
	old_par_child_num=0;
      else
	old_par_child_num=1;
      
      set_as_parent_child(old_par->get_parent(), old_sib, old_par_child_num);
      set_as_siblings(old_sib, old_par_sib);
      extra_branch = old_par;
      set_root(old_sib);
   
    }
    else 
      {
	old_sib->null_parent();
	old_sib->null_sibling();
	extra_branch = root;
	set_root(old_sib);
      }
    
    split_point->null_parent();
    split_point->null_sibling();
    new_tree=new Tree(curr_exchange, split_point);

    return(new_tree);
  }
  else
    return(0);

}



Branch * Tree::join_trees(Tree *tree1, Tree *tree2, Branch *join1, Branch *join2, int choosetree)
  //Takes pointers to two trees, with a branch on each to join the trees at.
  //Reroots the non-base tree at its join point, and then joins them
{
  int old_par_child_num;
  double brn_swap;
  Branch *tree_root, *add_to, *old_par, *old_sib, *old_join_null=0, *old_join_root=0;
  Tree *basetree, *jointree;


  //Descide which tree is the base
  if (choosetree == 0) {
    basetree=tree1;
    jointree=tree2;
    tree_root=join2;
    add_to=join1;
  }  
  else {
    basetree=tree2;
    jointree=tree1;
    tree_root=join1;
    add_to=join2;
  }

  
  if (jointree->find_null_branch() != 0) 
	old_join_null=jointree->find_null_branch();
  old_join_root=jointree->find_root();
  
  //Re-roots the non-base tree at it's join point
  if ((tree_root != jointree->find_root()) && (tree_root != jointree->find_null_branch_id()))
      jointree->re_root_tree(tree_root);

 
  //Finds the "real" root, whose child is tree_root
  if (jointree->find_root() != tree_root)
    tree_root=tree_root->get_parent();

  //Joins the two trees
  if (add_to != basetree->find_root())
    {
      old_sib=add_to->get_sibling();
      old_par=add_to->get_parent();
      
      if (old_par->get_child(0) == add_to)
	old_par_child_num=0;
      else
	old_par_child_num=1;
     
      set_as_siblings(add_to, tree_root);
      set_as_parent_child(extra_branch, add_to, 0);
      set_as_parent_child(extra_branch, tree_root, 1);
      set_as_siblings(extra_branch, old_sib);
      set_as_parent_child(old_par, extra_branch, old_par_child_num);
     
	  if (add_to->get_brnlen() == 0) {
		add_to->set_brnlen(extra_branch->get_brnlen());
		extra_branch->set_brnlen(0);
	  }

      basetree->set_root(add_to);
    
	}
  else {
    extra_branch->null_parent();
    extra_branch->null_sibling();
   

    set_as_parent_child(extra_branch, add_to, 0);
    set_as_parent_child(extra_branch, tree_root, 1);
    set_as_siblings(add_to, tree_root);
	if (add_to->get_child(0) != 0)
		if (add_to->get_child(0)->get_brnlen() == 0)
		{
			add_to->get_child(0)->set_brnlen(extra_branch->get_brnlen());
			extra_branch->set_brnlen(0);
		}
	
    basetree->set_root(add_to);
  }
	if (old_join_root->get_brnlen()==0)
	{
		old_join_root->set_brnlen(basetree->find_root()->get_brnlen());
		basetree->find_root()->set_brnlen(0);
	}
	if (old_join_null != 0)
		if (old_join_null->get_brnlen() == 0) {
			old_join_null->set_brnlen(basetree->find_null_branch()->get_brnlen());
			basetree->find_null_branch()->set_brnlen(0);
		}
	if((extra_branch->get_brnlen() == 0) &&(extra_branch != basetree->find_root()) && (extra_branch != basetree->find_null_branch()))
	{
		extra_branch->set_brnlen(basetree->find_null_branch()->get_brnlen());
		basetree->find_null_branch()->set_brnlen(0);
	}
	
  root=basetree->find_root();
  left_tip=basetree->get_leftmost_tip();

  diganose_tree(basetree->find_root());
  if ((basetree->find_root()->get_sibling() !=0 ) || (basetree->find_root()->get_parent() != 0))
    cerr<<"Error joining tree--root is not bottom: sib: "<<basetree->find_root()->get_sibling()<<" Parent: "
	<<basetree->find_root()->get_parent()<<"\n";
  return(add_to);
}



void Tree::set_num_nonsyn_params()
{
	int i, j;
	BOOL used;
	

	max_p_nonsyn=0;
	max_evol_pattern=new int [curr_exchange->get_num_nonsyn_params()];
	for(i=0; i<curr_exchange->get_num_nonsyn_params(); i++)
		max_evol_pattern[i]=0;


	switch(curr_exchange->get_model()) {
	case MG_94_JC:
	case MG_94_K2P:
	case MG_94_HKY:
		for(i=0; i<curr_exchange->get_num_branches(); i++) 
			if(tree[i]->get_p_nonsyn_num()+1>max_p_nonsyn) {
				max_p_nonsyn=tree[i]->get_p_nonsyn_num()+1;
				max_evol_pattern[0]=tree[i]->get_p_nonsyn_num();
			}
		
		break;
	case C_00_JC:
	case C_00_K2P:
	case C_00_HKY:
		for(i=0; i<curr_exchange->get_num_branches(); i++) {
			if(tree[i]->get_p_nonsyn_num()>max_evol_pattern[0])
				max_evol_pattern[0]=tree[i]->get_p_nonsyn_num();
			
			if(tree[i]->get_p_intergroup_num()>max_evol_pattern[1])
				max_evol_pattern[1]=tree[i]->get_p_intergroup_num();


			if(tree[i]->get_p_nonsyn_num()+1>max_p_nonsyn)
				max_p_nonsyn=tree[i]->get_p_nonsyn_num()+1;
			if(tree[i]->get_p_intergroup_num()+1>max_p_nonsyn)
				max_p_nonsyn=tree[i]->get_p_intergroup_num()+1;
		}
		break;
	case LCAP_JC:
	case LCAP_K2P:
	case LCAP_HKY:
			for(i=0; i<curr_exchange->get_num_branches(); i++) {
				for(j=0; j<curr_exchange->get_num_aa_props()+1; j++) {
					if(tree[i]->get_aa_prop_num(j)>max_evol_pattern[j])
						max_evol_pattern[j]=tree[i]->get_aa_prop_num(j);

					if(tree[i]->get_aa_prop_num(j) +1 > max_p_nonsyn)
						max_p_nonsyn=tree[i]->get_aa_prop_num(j) +1;
				}
			}
		break;	
	case AAMULMAT_JC:
	case AAMULMAT_K2P:
	case AAMULMAT_HKY:
		for(i=0; i<curr_exchange->get_num_branches(); i++) {
				for(j=0; j<curr_exchange->get_num_matrices(); j++) {
					if(tree[i]->get_matrix_coeff_num(j)>max_evol_pattern[j])
						max_evol_pattern[j]=tree[i]->get_aa_prop_num(j);

					if(tree[i]->get_matrix_coeff_num(j) +1 > max_p_nonsyn)
						max_p_nonsyn=tree[i]->get_matrix_coeff_num(j) +1;
				}
			}
			
		break;
	}


	get_num_distinct_evol_patterns();
	cout<<"Identified "<<num_patterns<<" distinct patterns of non-synonymous substitutions\n";
    curr_exchange->set_num_p_nonsyn(max_p_nonsyn, num_patterns);



	for(i=0; i<curr_exchange->get_num_branches(); i++)
		cout<<"Branch "<<i<<" Name: "<<tree[i]->get_name()<<" PATERN: "<<
		tree[i]->get_nonsyn_pattern()<<" PNS#: "<<tree[i]->get_p_nonsyn_num()
		<<" PIG:#: "<<tree[i]->get_p_intergroup_num()<<" AA0#:"<<tree[i]->get_aa_prop_num(0)<<endl;

	switch (curr_exchange->get_model()) {
	case C_00_JC:
	case C_00_K2P:
	case C_00_HKY:
		for(i=0; i<max_p_nonsyn; i++)
		{
			if (max_evol_pattern[0] >= i)
				used=TRUE;
			else
				used=FALSE;

			curr_exchange->set_p_non_syn_used(i, used);
			
			if (max_evol_pattern[1] >= i)
				used=TRUE;
			else
				used=FALSE;

			curr_exchange->set_p_inter_group_used(i, used);
		
		}
		break;
	case LCAP_JC:
	case LCAP_K2P:
	case LCAP_HKY:
		for(i=0; i<max_p_nonsyn; i++)
		{
			for(j=0; j<curr_exchange->get_num_nonsyn_params(); j++) {
				if (max_evol_pattern[j] >= i)
					used=TRUE;
				else
					used=FALSE;

				curr_exchange->set_aa_property_used(i, j, used);
			}
		
		}
		break;
	case AAMULMAT_JC:
	case AAMULMAT_K2P:
	case AAMULMAT_HKY:
		for(i=0; i<max_p_nonsyn; i++)
		{
			for(j=0; j<curr_exchange->get_num_nonsyn_params(); j++) {
				if (max_evol_pattern[j] >= i)
					used=TRUE;
				else
					used=FALSE;

				curr_exchange->set_matrix_coeff_used(i, j, used);
			}
		
		}	
		break;



	}



	delete[] max_evol_pattern;

}


void Tree::get_num_distinct_evol_patterns()
{
	int i, j, *patt, nonsyn_match;
	BOOL matches, same;
	EPattern_list *listpoint;

	num_patterns=1;
	the_patterns=0;

	if(curr_exchange->get_num_nonsyn_params() == 1) {
		num_patterns=max_p_nonsyn;
	}
	else {
		patt=new int[curr_exchange->get_num_nonsyn_params()];

		for(i=0; i<curr_exchange->get_num_branches(); i++) {
			switch(	curr_exchange->get_model()) {
			case C_00_JC:
			case C_00_K2P:
			case C_00_HKY:
				patt[0]=tree[i]->get_p_nonsyn_num();
				patt[1]=tree[i]->get_p_intergroup_num();
				break;
			case LCAP_JC:
			case LCAP_K2P:
			case LCAP_HKY:
				for(j=0; j<curr_exchange->get_num_nonsyn_params(); j++)
					patt[j]=tree[i]->get_aa_prop_num(j);
				break;
			case AAMULMAT_JC:
			case AAMULMAT_K2P:
			case AAMULMAT_HKY:
				for(j=0; j<curr_exchange->get_num_nonsyn_params(); j++)
					patt[j]=tree[i]->get_matrix_coeff_num(j);
				break;
				
				break;

			}
			if(the_patterns != 0) {
				listpoint=the_patterns;

				the_patterns=start_patt;

				matches=FALSE;
				while(the_patterns != 0) {
					same=TRUE;
					for(j=0; j<curr_exchange->get_num_nonsyn_params(); j++) {
						if (the_patterns->pattern[j] != patt[j])
							same=FALSE;
					}
					if (same == TRUE) {
						matches = TRUE;
						nonsyn_match=the_patterns->pattern_branch;
					}
				
					the_patterns=the_patterns->next;
				}

				the_patterns=listpoint;
				if (matches == FALSE) {
					num_patterns++;
					tree[i]->set_nonsyn_pattern(num_patterns-1);
					the_patterns->next=new EPattern_list;
					the_patterns->next->last=the_patterns;
					the_patterns=the_patterns->next;
					the_patterns->next=0;
					the_patterns->pattern=new int [curr_exchange->get_num_nonsyn_params()];
					the_patterns->pattern_branch=i;
					for(j=0; j<curr_exchange->get_num_nonsyn_params(); j++)
						the_patterns->pattern[j]=patt[j];

				}
				else
					tree[i]->set_nonsyn_pattern(tree[nonsyn_match]->get_nonsyn_pattern());
			}
			else {
				the_patterns=new EPattern_list;
				the_patterns->pattern=new int[curr_exchange->get_num_nonsyn_params()];
				the_patterns->pattern_branch=i;
				start_patt=the_patterns;
				the_patterns->last=the_patterns->next=0;
				for(j=0; j<curr_exchange->get_num_nonsyn_params(); j++)
					the_patterns->pattern[j]=patt[j];
				tree[i]->set_nonsyn_pattern(num_patterns-1);
			}

		}

		the_patterns=start_patt;

		while(the_patterns != 0)
		{
			delete[] the_patterns->pattern;
			listpoint=the_patterns->next;
			delete the_patterns;
			the_patterns=listpoint;
		}

	}


}

void Tree::ascend_assign_clade (Branch *other_branch, Branch *my_pos)
{
	int i;
	Branch *new_left, *new_right;
	
	cout<<"Initialiing "<<other_branch->get_name()<<endl;
	cout<<"My tip: "<<my_pos->is_tip()<<endl;
	
	if (other_branch->get_child(0)->is_tip() == FALSE) {
		new_left=initialize_branch();
	}
	else {
		for(i=0; i<curr_exchange->get_num_taxa(); i++)
			if (tree[i]->get_taxa_id() == other_branch->get_child(0)->get_taxa_id())
				new_left=tree[i];
	}
	
	if (other_branch->get_child(1)->is_tip() == FALSE) {
		new_right=initialize_branch();
	}
	else {
		for(i=0; i<curr_exchange->get_num_taxa(); i++)
			if (tree[i]->get_taxa_id() == other_branch->get_child(1)->get_taxa_id())
				new_right=tree[i];
	}
	
	cout<<"Old left: "<<other_branch->get_child(0)->get_taxa_id()<<": "<<other_branch->get_child(0)->get_name()<<endl;
	cout<<"Old right: "<<other_branch->get_child(1)->get_taxa_id()<<": "<<other_branch->get_child(1)->get_name()<<endl;
	cout<<"New left: "<<new_left->get_taxa_id()<<": "<<new_left->get_name()<<endl;
	cout<<"New right: "<<new_right->get_taxa_id()<<": "<<new_right->get_name()<<endl<<flush;
	
	
	set_as_parent_child(my_pos, new_left, 0);
	set_as_parent_child(my_pos, new_right, 1);
	set_as_siblings(new_left, new_right);
	
	if (new_left->is_tip() == FALSE)
		ascend_assign_clade(other_branch->get_child(0), new_left);
	
	if (new_right->is_tip() == FALSE)
		ascend_assign_clade(other_branch->get_child(1), new_right);
	
	
}

void add_to_branch_list(Branch_list *&list, Branch *new_ele)
{
	list->element=new_ele;
	if (list->next==0) {
		list->next=new Branch_list;
		list=list->next;
		list->next=0;
	}
	else
		list=list->next;
	list->element=0;		
}


void  merge_trees (Tree *tree1, Tree *tree2, Exchange *curr_exchange1,  Exchange *curr_exchange2, Tree *&new_tree, Exchange *&new_exchange)
{
    int i, j, brn_cnt=0,  my_child, my_sib, root1, root2;
    Branch **lookup_table, *my_match, *new_root;
    
    new_exchange=new Exchange();
    
    (*new_exchange)=(*curr_exchange1);
    new_exchange->set_num_taxa(curr_exchange1->get_num_taxa() + curr_exchange2->get_num_taxa());
    
    
    lookup_table = new Branch *[curr_exchange1->get_num_branches() +curr_exchange2->get_num_branches()];
    
    new_tree = new Tree(new_exchange, TRUE);
    
    
    for(i=0; i<curr_exchange1->get_num_branches(); i++) {
        (*(*new_tree)[brn_cnt])=(*(*tree1)[i]);
        lookup_table[brn_cnt]=(*tree1)[i];
        (*new_tree)[brn_cnt]->null_child(0);
        (*new_tree)[brn_cnt]->null_child(1);
        (*new_tree)[brn_cnt]->null_parent();
        (*new_tree)[brn_cnt]->null_sibling();
        
        brn_cnt++;
    }
    
    for(i=0; i<curr_exchange2->get_num_branches(); i++) {
        (*(*new_tree)[brn_cnt])=(*(*tree2)[i]);
        lookup_table[brn_cnt]=(*tree2)[i];
        (*new_tree)[brn_cnt]->null_child(0);
        (*new_tree)[brn_cnt]->null_child(1);
        (*new_tree)[brn_cnt]->null_parent();
        (*new_tree)[brn_cnt]->null_sibling();
        
        brn_cnt++;
    }
    
    
    for(i=0; i<new_exchange->get_num_branches(); i++) {
        
        if ((*new_tree)[i]->is_uninitialized() == FALSE) {
            my_match = lookup_table[i];
            
            if ((*new_tree)[i]->is_tip() == FALSE) {
                for(j=0; j<2; j++) {
                    if ((*new_tree)[i]->get_child(j)==0) {
                        my_child=0;
                        while(lookup_table[my_child] != my_match->get_child(j)) my_child++;
                        
                        new_tree->set_as_parent_child((*new_tree)[i], (*new_tree)[my_child], j);
                        
                    }
                }
            }
            
            if (((*new_tree)[i]->get_sibling() ==0) && (my_match->get_parent() != 0)) {
                my_sib=0;
                while(lookup_table[my_sib] != my_match->get_sibling()) my_sib++;
                
                new_tree->set_as_siblings((*new_tree)[i], (*new_tree)[my_sib]);
            }
        }
        else {
            new_root = (*new_tree)[i];
            new_root->initialized();
            new_root->set_brnlen(0.0);
        }
    }
    
    
    root1=0;
    while(lookup_table[root1] != tree1->find_root()) root1++;
    
    root2=0;
    while(lookup_table[root2] != tree2->find_root()) root2++;
    
    new_tree->set_as_parent_child(new_root, (*new_tree)[root1], 0);
    new_tree->set_as_parent_child(new_root, (*new_tree)[root2], 1);
    delete lookup_table;
    
}

void double_tree (Tree *tree, Exchange *curr_exchange, Tree *&new_tree, Exchange *&new_exchange, int level)
{
    int i, j, brn_cnt=0,  my_child, my_sib, root_cnt=0, *brn_level, *brn_backref, **brn_forward_ref, my_root, last_root;
    string name;
    Branch ***lookup_table, *my_match, **new_roots, *last_root_brn;
    
    new_exchange=new Exchange();
    
    (*new_exchange)=(*curr_exchange);
    new_exchange->set_num_taxa(level*curr_exchange->get_num_taxa() );
    
    new_roots= new Branch *[level-1];
    
    lookup_table = new Branch **[level];
    
    for(j=0; j<level; j++)
        lookup_table[j] = new Branch *[curr_exchange->get_num_branches()];
    
    new_tree = new Tree(new_exchange, TRUE);
    brn_level = new int [new_exchange->get_num_branches()];
    brn_backref = new int [new_exchange->get_num_branches()];
    
    brn_forward_ref=new int * [level];
    for(j=0; j<level; j++) brn_forward_ref[j]=new int [curr_exchange->get_num_branches()];
    
    for(j=0; j<level; j++) {
        for(i=0; i<curr_exchange->get_num_branches(); i++) {
            (*(*new_tree)[brn_cnt])=(*(*tree)[i]);
            if ((*new_tree)[brn_cnt]->is_tip() == TRUE)
                (*new_tree)[brn_cnt]->set_taxa_id(j*curr_exchange->get_num_taxa() + (*tree)[i]->get_taxa_id() );
            (*new_tree)[brn_cnt]->set_brn_id(brn_cnt);
            brn_level[brn_cnt]=j;
            brn_backref[brn_cnt]=i;
            brn_forward_ref[j][i]=brn_cnt;
            lookup_table[j][i]=(*tree)[i];
            (*new_tree)[brn_cnt]->null_child(0);
            (*new_tree)[brn_cnt]->null_child(1);
            (*new_tree)[brn_cnt]->null_parent();
            (*new_tree)[brn_cnt]->null_sibling();
            
            brn_cnt++;
        }
    }
    
    
    for(i=0; i<new_exchange->get_num_branches(); i++) {
        
        if ((*new_tree)[i]->is_uninitialized() == FALSE) {
            my_match = lookup_table[brn_level[i]][brn_backref[i]];
            
            if ((*new_tree)[i]->is_tip() == FALSE) {
                
                for(j=0; j<2; j++) {
                    //cout<<"Branch "<<i<<" looking for child "<<j<<endl;
                    if ((*new_tree)[i]->get_child(j)==0) {
                        my_child=0;
                        while(lookup_table[brn_level[i]][my_child] != my_match->get_child(j)) my_child++;
                        
                        //cout<<"Child j is "<<my_match->get_child(j)->get_name()<<endl;
                        
                        new_tree->set_as_parent_child((*new_tree)[i], (*new_tree)[brn_forward_ref[brn_level[i]][my_child]], j);
                        
                    }
                }
            }
            
            if (((*new_tree)[i]->get_sibling() ==0) && (my_match->get_parent() != 0)) {
                my_sib=0;
                while(lookup_table[brn_level[i]][my_sib] != my_match->get_sibling()) my_sib++;
                
                new_tree->set_as_siblings((*new_tree)[i], (*new_tree)[brn_forward_ref[brn_level[i]][my_sib]]);
            }
        }
        else {
            new_roots[root_cnt] = (*new_tree)[i];
            new_roots[root_cnt]->initialized();
            new_roots[root_cnt]->set_brnlen(0.0);
            new_roots[root_cnt]->null_parent();
            new_roots[root_cnt]->null_sibling();
            new_roots[root_cnt]->set_brn_id(i);
            root_cnt++;
        }
    }
    
    last_root=0;
    while(lookup_table[0][last_root] != tree->find_root()) last_root++;
    
    last_root_brn=(*new_tree)[brn_forward_ref[0][last_root]];
    
    (*new_tree)[brn_forward_ref[0][last_root]]->set_name("GT1ROOT");
    for(j=1; j<level; j++) {
        my_root=0;
        while(lookup_table[j][my_root] != tree->find_root()) my_root++;
        

        (*new_tree)[brn_forward_ref[j][my_root]]->set_name("GT2Root");
        new_tree->set_as_siblings(last_root_brn, (*new_tree)[brn_forward_ref[j][my_root]]);
        new_tree->set_as_parent_child(new_roots[j-1], last_root_brn, 0);
        new_tree->set_as_parent_child(new_roots[j-1], (*new_tree)[brn_forward_ref[j][my_root]], 1);
        //last_root=my_root;
        new_roots[j-1]->set_name("Root");
        last_root_brn=new_roots[j-1];
        
    }
    
    //cout<<"New root p "<<new_roots[level-1]<<" levle "<<level<<endl;
    new_tree->set_root(new_roots[level-2]);
#if 0
    for(i=0; i<new_exchange->get_num_branches(); i++) {
        cout<<"Brn# "<<i<<" = "<<(*new_tree)[i]->get_brn_num()<<" Name: "<<(*new_tree)[i]->get_name();
        if ((*new_tree)[i]->get_sibling() != 0) cout<<" sib: "<<(*new_tree)[i]->get_sibling()->get_name();
        if ((*new_tree)[i]->get_parent() != 0) cout<<" par: "<<(*new_tree)[i]->get_parent()->get_name();
        cout<<endl;
     }
#endif
    
    for(j=0; j<level; j++) {
          delete[] lookup_table[j];
          delete[] brn_forward_ref[j];
    }
    
    delete[] brn_forward_ref;
    delete[] lookup_table;
    delete[] brn_level;
    delete[] brn_backref;
    delete[] new_roots;
    
    
}


