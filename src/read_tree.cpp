#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <cctype>

#include "read_tree.h"

using namespace::std;
//#define DEBUG


Branch* Read_Tree::initialize_branch ()
  //This function is used to create new internal branches.  It finds the first 
  //"unused" branch in the current tree array, gives it the next unused id number 
  //and returns a pointer to it.

{
 
  return(curr_tree->initialize_branch());
}  //End Read_PAUP_Tree::initialize_branch




Branch * Read_Tree::add_to_branch(Branch *add_to, Branch *add)
  //This function inserts a new branching along add_to
{
  int add_to_child_num;
  Branch *new_par, *add_to_par;

  new_par=initialize_branch();
 
  if (add_to->get_parent()!=0)
    {
      if (add_to->get_parent()->get_child(0)==add_to)
	add_to_child_num=0;
      else
	add_to_child_num=1;
      
      add_to_par=add_to->get_parent();

      curr_tree->set_as_parent_child(add_to_par, new_par, add_to_child_num);
      curr_tree->set_as_siblings(new_par, add_to_par->get_child(curr_tree->other_child(add_to_child_num)));
      curr_tree->set_as_parent_child(new_par, add_to, add_to_child_num);
      curr_tree->set_as_parent_child(new_par, add, curr_tree->other_child(add_to_child_num));
      curr_tree->set_as_siblings(add_to, add);

      
    }
  else
    {
      curr_tree->set_as_parent_child(new_par, add_to->get_child(0), 0);
      curr_tree->set_as_parent_child(new_par, add_to->get_child(1), 1);
      curr_tree->set_as_parent_child(add_to, new_par, 0);
      curr_tree->set_as_parent_child(add_to, add, 1);
      curr_tree->set_as_siblings(add, new_par);
    }
  return(new_par);
} //End  Read_Tree::add_to_branch



Tree * Read_PAUP_Tree::create_tree_from_file (Exchange *cexchange, Sequence_dataset *curr_data)
{
	return(create_tree_from_file(cexchange, curr_data, cexchange->is_rooted_tree()));
}


Tree * Read_PAUP_Tree::create_tree_from_file (Exchange *cexchange, Sequence_dataset *curr_data, BOOL rooted)

  //This function receives the name of file containing a tree description in Paup 4.0 format.
  //An array of the structure Branch is created with 2n-1 branches (for n taxa).  This includes 
  //a 0 length branch at the root specified by Paup.
  //A partially recursive method is used to read in the tree description:
  //   1. The first tip is read in either alone or by Tree::new_subtree if it is a nested branch
  //      a) Tree::new_subtree will return the ancestor of the tip that is read when it is called
  //   2. Tree::get_sibling looks to see if the next part of the tree is a sibling or another
  //      new subtree.  If a sibling, Tree::new_sibling is called, otherwise Tree::new_subtree   
  //   3. Tree::new_parent is used by get_sibling, new_subtree and this function to 
  //      create parents
  //This process repeats until a semicolon is found in the tree description for rooted trees
  //Unrooted trees are read in three parts corrisponding the the three branches at the base.
  //One of these base branches is converted into a peudo-root by division in half

{
    Branch  *curbranch, *unroot_sib;
    ifstream treein;
    char dump;
    string myline, myname, kappa, param;
    int i,j=0,k, place=0, level=0, curid, group, num_groups=0, rate_num=0, cnt_pns=0;
    double param_val;
    //BOOL rooted=TRUE;
    Groups_list *start, *the_list=0, *temp;
    std::size_t found;

    curr_exchange=cexchange;

    curr_tree=new Tree(curr_exchange, rooted);

    treein.open(curr_exchange->get_treefile());
    if (treein.fail()) {
        cerr<<"Cannot find file "<<curr_exchange->get_treefile()<<endl;
        return(0);
    }
    else {

      //treein.getline(line, 399);
      getline(treein, myline);
      std::transform(myline.begin(), myline.end(), myline.begin(), ::tolower);
        
      while (myline.find("translate") == std::string::npos) {
            myline="";
            getline(treein, myline);
            std::transform(myline.begin(), myline.end(), myline.begin(), ::tolower);
      
            if (myline.find("assumed nucleotide frequencies") != std::string::npos) {
                for (i=0; i<4; i++) {
                    treein>>myline;
                    treein>>myline;
                    treein>>myline;
                    treein>>myline;
                    if (curr_exchange->fixed_basefreq() == FALSE)
                            curr_exchange->set_basefreqs (string_to_float(myline.c_str()), i);
                }
                treein.get(dump);
            }
      
            if (myline.find("estimated base frequencies") != std::string::npos) {
                for (i=0; i<4; i++) {
                    getline(treein, myline);
                    if (curr_exchange->fixed_basefreq() == FALSE)
                        curr_exchange->set_basefreqs (string_to_float(myline.c_str()), i);
                }
            }
                

                
            if ((myline.find("transition/transversion ratio = ")!= std::string::npos) ||
               (myline.find("Estimated ti/tv ratio =") != std::string::npos))  {
                if (myline.find("(kappa =")!= std::string::npos) {
                    kappa="";
                    found=myline.find("(");
                    kappa=myline.substr(found+1);
                    found =kappa.find(")");
                    kappa=kappa.substr(0, found-1);
                
                    curr_exchange->set_trs_trv(string_to_float(kappa.c_str()));
                }
                else
                    curr_exchange->set_obs_trs_trv(string_to_float(myline.c_str()));
            }


            if (myline.find("prob. intra-group substitution")!= std::string::npos) {
                found=myline.find(":");
                param=myline.substr(found+1);
                curr_exchange->set_p_inter_group(0, string_to_float(param.c_str()));
            }

            if ((myline.find("prob. extra-group substitution") != std::string::npos) ||
                (myline.find("prob. non-synonymous substitution")!= std::string::npos)) {
                found=myline.find(":");
                param=myline.substr(found+1);
                if (cnt_pns >= curr_exchange->get_num_rates())
                    curr_exchange->set_p_non_syn(0, string_to_float(param.c_str()));
                else {
                    curr_exchange->set_p_non_syn(0, cnt_pns, string_to_float(param.c_str()));
                    cnt_pns++;
                }
            }

            if (myline.find("chemical composition factor")!= std::string::npos) {
                found=myline.find(":");
                param=myline.substr(found+1);
                curr_exchange->set_aa_prop_fact(0, rate_num, CHEM_COMP, string_to_float(param.c_str()));
            }

            if (myline.find("Polarity Factor")!= std::string::npos) {
                found=myline.find(":");
                param=myline.substr(found+1);
                curr_exchange->set_aa_prop_fact(0, rate_num, POLARITY, string_to_float(param.c_str()));
            }
            if (myline.find("volume factor")!= std::string::npos) {
                found=myline.find(":");
                param=myline.substr(found+1);
                curr_exchange->set_aa_prop_fact(0, rate_num, VOLUME, string_to_float(param.c_str()));
            }
            if (myline.find("iso-electric point factor") != std::string::npos) {
                found=myline.find(":");
                param=myline.substr(found+1);
                cout<<"Found iso: "<<rate_num<<": "<<string_to_float(param.c_str())<<endl;
                curr_exchange->set_aa_prop_fact(0, rate_num, ISO_ELEC, string_to_float(param.c_str()));
            }
            if (myline.find("hydropathy point factor")!= std::string::npos) {
                found=myline.find(":");
                param=myline.substr(found+1);
                param_val=string_to_float(param.c_str());
                //Nasty hack to avoid NaN problems in like_rate_partition
                if (param_val < -2)
                      param_val=-2;
                curr_exchange->set_aa_prop_fact(0, rate_num, HYDROPATHY, param_val);
            }
            if (myline.find("scaling factor")!= std::string::npos) {
                found=myline.find(":");
                param=myline.substr(found+1);
                curr_exchange->set_aa_prop_fact(0, rate_num, SCALING, string_to_float(param.c_str()));
                if (rate_num < curr_exchange->get_num_rates())
                    rate_num++;
            }
            if (myline.find("instantaneous rate of duplicate fixation") != std::string::npos) {
                found=myline.find(":");
                param=myline.substr(found+1);
                curr_exchange->set_dupl_fix_rate( string_to_float(param.c_str()));
            }
            if (myline.find("instantaneous rate of convergence to single duplicate") != std::string::npos) {
                found=myline.find(":");
                param=myline.substr(found+1);
                curr_exchange->set_dupl_parallel_rate( string_to_float(param.c_str()));
            }
            if (myline.find("relative rate of duplicate fixation after converg.")!= std::string::npos) {
                found=myline.find(":");
                param=myline.substr(found+1);
                curr_exchange->set_fix_rate_scale( string_to_float(param.c_str()));

            }
            if (myline.find("relative rate of duplicate loss after converg.")!= std::string::npos) {
                found=myline.find(":");
                param=myline.substr(found+1);
                curr_exchange->set_loss_rate_scale( string_to_float(param.c_str()));
            }
            if (myline.find("track switch between genes probability:") != std::string::npos) {
                found=myline.find(":");
                param=myline.substr(found+1);
                curr_exchange->set_strand_switch_prob(string_to_float(param.c_str()));
            }
            if (myline.find("ratio of snp losses to snp state changes:")!= std::string::npos) {
                found=myline.find(":");
                param=myline.substr(found+1);
                curr_exchange->set_snp_loss_to_switch(string_to_float(param.c_str()));
            }
            if (myline.find("ratio of snp gains to losses:") != std::string::npos) {
                found=myline.find(":");
                param=myline.substr(found+1);
                curr_exchange->set_snp_gains_to_losses(string_to_float(param.c_str()));
            }

            if (myline.find("groups:")!= std::string::npos) {
                found=myline.find(": ");
                param=myline.substr(found+2);
                //cout<<"Partial group line |"<<param<<"|"<<endl;
                found=param.find(":");
                
                param=param.substr(0, found);
                //cout<<"After : clean: "<<found<<": "<<param<<endl;
                group=string_to_int(param.c_str());

                if (group+1 >num_groups) num_groups=group+1;
                std::cout<<"Found group "<<group<<endl;

                temp=the_list;
                the_list=new Groups_list;
                the_list->next=0;
                for(j=0; j<20; j++)
                    the_list->aa[j]=0;

                if (temp == 0)
                    start=the_list;
                else {
                    the_list->last=temp;
                    temp->next=the_list;
                }
                i=0;
                found=myline.find(": ");
                param=myline.substr(found+2);
                while(i < param.length()) {
                    //cout<<"Checking "<<param[i]<<" for group "<<group<<endl;
                    if (((param[i]>=65 && param[i]<=90) || (param[i]>=97 && param[i]<=122)) && (is_aa(param[i])==TRUE))
                        the_list->aa[readchar_to_aa(param[i])]=1;
                    i++;
                }

            }
      }
      
      if(num_groups != 0) {
          curr_exchange->curr_groups=new Amino_acid_group(num_groups);

          the_list=start;
          for(i=0; i<num_groups; i++)
            {
              //cout<<" Current list: "<<the_list<<endl;
              for(j=0; j<20; j++)
            if(the_list->aa[j]==1)
              curr_exchange->curr_groups->assign_to_group(j, i);
              the_list=the_list->next;
            }
          the_list=start;
          while (the_list !=0) {
              temp=the_list->next;
              delete the_list;
              the_list=temp;
            }
          for(i=0; i<num_groups; i++) {
              cout<<"Group "<<i<<": ";
              for(j=0; j<20; j++)
                  if(curr_exchange->curr_groups->get_group(j) == i)
                      cout<<num_to_aa(j)<<" ";
              cout<<endl;
          }
      }
      
      
      //Uses Paup's translate table to name the tip branches

        for (i=0; i<curr_exchange->get_num_taxa(); i++) {
            treein>>curid>>myname;
            cout<<i<<": Read id as "<<curid<<" and the name as "<<myname<<endl;
            found=myname.find(",");
            if (found != std::string::npos) myname=myname.substr(0, found);
            found=myname.find(";");
            if (found != std::string::npos) myname=myname.substr(0, found);
            
            found=myname.find("\'");
            if (found != std::string::npos) {
                //cout<<"We have a quote:. position; "<<found<<endl;
                if (found <1) myname=myname.substr(found+1);
                else myname=myname.substr(0, found-1);
            }
            found=myname.find("\'");
            if (found != std::string::npos) myname=myname.substr(0, found);
            
            #if defined (DEBUG)
            cout<<"Read "<<myname<<" for taxa "<<i<<endl;
            #endif

            (*curr_tree)[i]->set_name(myname.c_str());
            (*curr_tree)[i]->initialized();

            if (curr_exchange->have_data()==TRUE) {
              j=0;
              
              while (strcmp((*curr_data)[j].Sequence_name(), myname.c_str())!=0 && j<curr_exchange->get_num_taxa())
              j++;
             
              if (j==curr_exchange->get_num_taxa())
                  cerr<<"Number of taxa exceeded and sequence "<<myname<<" still not found in sequence file\n";
              
              (*curr_tree)[i]->set_taxa_id(j);
            }
            else
                (*curr_tree)[i]->set_taxa_id(i);
        }
      
        //treein.get(dump);
        //treein.get(dump);
        treein>>myline;
        cout<<"Next line character was "<<myline<<endl;
        treein>>myline;
        treein>>myline;
        treein>>myline;
        treein>>myline;
      
        if (myline.find("[&R]") == std::string::npos)
            curr_tree->set_rooted_tree(FALSE);
        else
            curr_tree->set_rooted_tree(TRUE);
		
#if defined (DEBUG)
		if (curr_tree->rooted_tree()==TRUE) cout<<"ROOTED TREE\n";
		else cout<<"UNROOTED TREE\n";
        cout<<"Last line read was "<<myline<<endl;
#endif
		
        nest_level=0;
        //****************ROOTED*****************
        //Handles rooted trees
        if (curr_tree->rooted_tree()==TRUE) {
			treein>>rel;
	  
			curbranch=make_subtree(treein);
	  
#if defined (DEBUG)
			cout<<curbranch->get_name()<<" Length: "<<curbranch->expect_subs_site()<<endl;
#endif
		 
      
        } //End Rooted Section
      
      
      //****************UNROOTED*****************
      //Handles unrooted trees by finding the three siblings at the base and joining them
        else {
            treein>>rel;
	  
            curbranch=make_base_unrooted(treein);
	  
#if defined (DEBUG)
            cout<<curbranch->get_name()<<" Length: "<<curbranch->get_brnlen()<<endl;
#endif
	  
        } //End Unrooted Section
      
      
        curr_tree->set_root((*curr_tree)[0]);
      
        return(curr_tree);
    }
} //End Read_PAUP_Tree::create_tree_from_file



Branch* Read_PAUP_Tree::get_tip(ifstream &treein)
{
  int curid;
  double blen;
  Branch *new_tip;

  curid=rel-48;

  //Handles Taxa IDs greater than 9
  treein>>rel1;

  while(rel1>=48 && rel1<=57)
    {
      curid=curid*10+(rel1-48);
      treein>>rel1;
    }
  
  new_tip=(*curr_tree)[curid-1];
  new_tip->set_tip(TRUE);
  treein>>blen;
  new_tip->set_expect_subs_site(blen);
  
  return(new_tip);

}  //End Read_PAUP_Tree::get_tip



Branch* Read_PAUP_Tree::get_interior(ifstream &treein)
{
  double blen;
  Branch *curbranch;
  
  treein>>rel1;
  treein>>rel1;

  curbranch=initialize_branch();
  curbranch->set_tip(FALSE);
  treein>>blen;
  curbranch->set_expect_subs_site(blen);
  
  return(curbranch);
}  //End Read_PAUP_Tree::get_interior




Branch* Read_PAUP_Tree::make_subtree(ifstream &treein)
{
  int i;
  double brn;
  Branch *sibling1, *sibling2, *parent;
  
  if (rel=='(')
    nest_level++;
  
  treein>>rel;
  
  if (rel=='(')
    sibling1=make_subtree(treein);                                     //Recursion
  else
    sibling1=get_tip(treein);

  treein>>rel;
  treein>>rel;
 
  if (rel=='(')
    sibling2=make_subtree(treein);                                 //Recursion
  else
    sibling2=get_tip(treein);

#if defined (DEBUG)
  cout<<" Sibling 1: "<<sibling1->get_name()<<" Sibling 2: "<<sibling2->get_name()<<" Nest level: "<<nest_level<<" \n"<<flush;
#endif 

  curr_tree->set_as_siblings(sibling1, sibling2);
  
  if (nest_level>1)
    parent=get_interior(treein);
  else
    {
      parent=initialize_branch();
      parent->set_tip(FALSE);
	  if (curr_tree->rooted_tree() == FALSE) {
		parent->set_brnlen(0.0);
		parent->set_expect_subs_site(0.0);      
	  }
	  else {
		treein>>rel;
		treein>>rel;
		treein>>brn;
		parent->set_expect_subs_site(brn);
	  }
      parent->set_parent(0);
      parent->set_name("Root");

 }
  curr_tree->set_as_parent_child(parent, sibling1, 0);
  curr_tree->set_as_parent_child(parent, sibling2, 1);
 
  nest_level--;
  return (parent);
}  //End Read_PAUP_Tree::make_subtree



Branch* Read_PAUP_Tree::make_base_unrooted(ifstream &treein)
{
  int i;
  Branch *sibling1, *sibling2, *sibling3, *parent;
  
  if (rel=='(')
    nest_level++;
  
  treein>>rel;
  
  if (rel=='(')
   sibling1=make_subtree(treein);                                     //Recursion
  else
    sibling1=get_tip(treein);

  treein>>rel;
  treein>>rel;
 
  if (rel=='(')
    sibling2=make_subtree(treein);                                 //Recursion
  else
    sibling2=get_tip(treein);

  treein>>rel;
  treein>>rel;
 
  if (rel=='(')
    sibling3=make_subtree(treein);                                 //Recursion
  else
    sibling3=get_tip(treein);

#if defined (DEBUG)
  cout<<"Unrooted base: Sibling 1: "<<sibling1->get_name()<<" Sibling 2: "<<sibling2->get_name()
      <<" Sibling 3: "<<sibling3->get_name()<<endl;
#endif 

  curr_tree->set_as_siblings(sibling1, sibling2);
  
  parent=initialize_branch();
  parent->set_tip(FALSE);
  parent->set_brnlen(0.0);    
  parent->set_expect_subs_site(0.0);

  curr_tree->set_as_parent_child(parent, sibling1, 0);
  curr_tree->set_as_parent_child(parent, sibling2, 1);
  
  sibling1=sibling3;
  sibling2=parent;

  curr_tree->set_as_siblings(sibling1, sibling2);
 
  parent=initialize_branch();
  parent->set_tip(FALSE);
  parent->set_brnlen(0.0);
  parent->set_expect_subs_site(0.0);
  parent->set_name("Root");

   curr_tree->set_as_parent_child(parent, sibling1, 1);
   curr_tree->set_as_parent_child(parent, sibling2, 0);

#if defined (DEBUG)
  cout<<" Parent: "<<parent->get_name()<<")";
#endif

  nest_level--;
  return (parent);
}  //End Read_PAUP_Tree::make_base_unrooted



Tree * Two_Taxa_Tree::create_two_taxa_tree(Exchange *cexchange, double dist)
{
  curr_tree=new Tree(cexchange);

  (*curr_tree)[0]->set_taxa_id(0);
  (*curr_tree)[1]->set_taxa_id(1);

  (*curr_tree)[0]->set_name("Taxa1");
  (*curr_tree)[1]->set_name("Taxa2");
  (*curr_tree)[2]->set_name("Root");

  (*curr_tree)[0]->initialized();
  (*curr_tree)[1]->initialized();
  (*curr_tree)[2]->initialized();


  curr_tree->set_as_parent_child((*curr_tree)[2], (*curr_tree)[0], 0);
  curr_tree->set_as_parent_child((*curr_tree)[2], (*curr_tree)[1], 1);
  curr_tree->set_as_siblings((*curr_tree)[0], (*curr_tree)[1]);


  (*curr_tree)[0]->set_expect_subs_site(dist);
  (*curr_tree)[1]->set_expect_subs_site(0.0);
  (*curr_tree)[2]->set_expect_subs_site(0.0);

  (*curr_tree)[0]->set_tip(TRUE);
  (*curr_tree)[1]->set_tip(TRUE);
  (*curr_tree)[2]->set_tip(FALSE);

  curr_tree->set_root((*curr_tree)[0]);
  return(curr_tree);
}




Tree * Three_Taxa_Tree::create_three_taxa_tree(Exchange *cexchange, Sequence_dataset *curr_data, double dist[3])
{
  int i;
  curr_tree=new Tree(cexchange);


  (*curr_tree)[0]->set_taxa_id(0);
  (*curr_tree)[1]->set_taxa_id(1);
  (*curr_tree)[2]->set_taxa_id(2);

  (*curr_tree)[0]->set_name((*curr_data)[0].Sequence_name());
  (*curr_tree)[1]->set_name((*curr_data)[1].Sequence_name());
  (*curr_tree)[2]->set_name((*curr_data)[2].Sequence_name());
  (*curr_tree)[3]->set_name("Null");
  (*curr_tree)[4]->set_name("Root");

  (*curr_tree)[0]->initialized();
  (*curr_tree)[1]->initialized();
  (*curr_tree)[2]->initialized();
  (*curr_tree)[3]->initialized();
  (*curr_tree)[4]->initialized();

 
  curr_tree->set_as_parent_child((*curr_tree)[3], (*curr_tree)[0], 0);
  curr_tree->set_as_parent_child((*curr_tree)[3], (*curr_tree)[1], 1);
  curr_tree->set_as_siblings((*curr_tree)[0], (*curr_tree)[1]);
  curr_tree->set_as_parent_child((*curr_tree)[4], (*curr_tree)[2], 0);
  curr_tree->set_as_parent_child((*curr_tree)[4], (*curr_tree)[3], 1);
  curr_tree->set_as_siblings((*curr_tree)[2], (*curr_tree)[3]);
  (*curr_tree)[4]->null_parent();
  (*curr_tree)[4]->null_sibling();

 
  (*curr_tree)[0]->set_expect_subs_site(dist[0]);
  (*curr_tree)[1]->set_expect_subs_site(dist[1]);
  (*curr_tree)[2]->set_expect_subs_site(dist[2]);
  (*curr_tree)[3]->set_expect_subs_site(0.0);
  (*curr_tree)[4]->set_expect_subs_site(0.0);
  
  for(i=0; i<cexchange->get_num_branches(); i++)
    (*curr_tree)[i]->set_brnlen(0.0);

  (*curr_tree)[0]->set_tip(TRUE);
  (*curr_tree)[1]->set_tip(TRUE);
  (*curr_tree)[2]->set_tip(TRUE);
  (*curr_tree)[3]->set_tip(FALSE);
  (*curr_tree)[4]->set_tip(FALSE);

  curr_tree->set_root((*curr_tree)[0]);
  return(curr_tree);
}


Read_PAUP_w_Settings_Tree::Read_PAUP_w_Settings_Tree()
{
	cerr<<"Error: cannot use default constructor of Read_PAUP_w_Settings_Tree\n";
}

Read_PAUP_w_Settings_Tree::Read_PAUP_w_Settings_Tree(Exchange *cexchange)
{
	int i;

	curr_exchange=cexchange;
	evol_pattern=new int [curr_exchange->get_num_nonsyn_params()];

}





Read_PAUP_w_Settings_Tree::~Read_PAUP_w_Settings_Tree()
{
	delete[] evol_pattern;
}


Branch* Read_PAUP_w_Settings_Tree::get_tip(ifstream &treein)
{
  int i, curid;
  char dump;
  double blen;
  Branch *new_tip;

  curid=rel-48;

  //Handles Taxa IDs greater than 9
  treein>>rel1;

  while(rel1>=48 && rel1<=57)
    {
      curid=curid*10+(rel1-48);
      treein>>rel1;
    }
  
  new_tip=(*curr_tree)[curid-1];
  new_tip->set_tip(TRUE);
  treein>>blen>>dump;
  for(i=0; i<curr_exchange->get_num_nonsyn_params()-1; i++)
	  treein>>evol_pattern[i]>>dump;
  i=curr_exchange->get_num_nonsyn_params()-1;
  treein>>evol_pattern[i];
  
  switch(curr_exchange->get_model()) {
  case MG_94_JC:
  case MG_94_K2P:
  case MG_94_HKY:
	new_tip->set_p_nonsyn_num(evol_pattern[0]-1);
	break;
  case C_00_JC:
  case C_00_K2P:
  case C_00_HKY:
	new_tip->set_p_nonsyn_num(evol_pattern[0]-1);
	new_tip->set_p_inter_group_num(evol_pattern[1]-1);
	break;
  case LCAP_JC:
  case LCAP_K2P:
  case LCAP_HKY:
	  for(i=0; i<curr_exchange->get_num_nonsyn_params(); i++)
		  new_tip->set_aa_prop_num(evol_pattern[i]-1, curr_exchange->get_prop_num_n(i));
	  break;
  case AAMULMAT_JC:
  case AAMULMAT_K2P:
  case AAMULMAT_HKY:		
	  for(i=0; i<curr_exchange->get_num_nonsyn_params(); i++)
		  new_tip->set_matrix_coeff_num(i, evol_pattern[i]-1);
	  break;


  }

  new_tip->set_expect_subs_site(blen);
  
  
  return(new_tip);

}  //End Read_PAUP_Tree::get_tip



Branch* Read_PAUP_w_Settings_Tree::get_interior(ifstream &treein)
{
  int i;
  char dump;
  double blen;
  Branch *curbranch;
  
  treein>>rel1;
  treein>>rel1;

  curbranch=initialize_branch();
  curbranch->set_tip(FALSE);
  treein>>blen>>dump;
  for(i=0; i<curr_exchange->get_num_nonsyn_params()-1; i++)
	  treein>>evol_pattern[i]>>dump;
  i=curr_exchange->get_num_nonsyn_params()-1;
  treein>>evol_pattern[i];
  
  switch(curr_exchange->get_model()) {
  case MG_94_JC:
  case MG_94_K2P:
  case MG_94_HKY:
	curbranch->set_p_nonsyn_num(evol_pattern[0]-1);
	break;
  case C_00_JC:
  case C_00_K2P:
  case C_00_HKY:
	curbranch->set_p_nonsyn_num(evol_pattern[0]-1);
	curbranch->set_p_inter_group_num(evol_pattern[1]-1);
	break;
  case LCAP_JC:
  case LCAP_K2P:
  case LCAP_HKY:
	  for(i=0; i<curr_exchange->get_num_nonsyn_params(); i++)
		  curbranch->set_aa_prop_num(evol_pattern[i]-1, curr_exchange->get_prop_num_n(i));
	  break;
  case AAMULMAT_JC:
  case AAMULMAT_K2P:
  case AAMULMAT_HKY:		
	  for(i=0; i<curr_exchange->get_num_nonsyn_params(); i++)
		  curbranch->set_matrix_coeff_num(i, evol_pattern[i]-1);
	  break;
  }
  
  curbranch->set_expect_subs_site(blen);
  
 
  return(curbranch);
}  //End Read_PAUP_Tree::get_interior



Read_PAUP_w_brlen_CI_Tree::Read_PAUP_w_brlen_CI_Tree()
{
	cerr<<"Error: cannot use default constructor of Read_PAUP_w_brlen_CI_Tree\n";
}

Read_PAUP_w_brlen_CI_Tree::Read_PAUP_w_brlen_CI_Tree(Exchange *cexchange)
{
	int i;

	curr_exchange=cexchange;
}





Read_PAUP_w_brlen_CI_Tree::~Read_PAUP_w_brlen_CI_Tree()
{
	
}


Branch* Read_PAUP_w_brlen_CI_Tree::get_tip(ifstream &treein)
{
  int i, curid;
  char dump;
  double blen;
  Branch *new_tip;

  curid=rel-48;

  //Handles Taxa IDs greater than 9
  treein>>rel1;

  while(rel1>=48 && rel1<=57)
    {
      curid=curid*10+(rel1-48);
      treein>>rel1;
    }
  
	cout<<"Reading tip: "<<curid<<endl;

  new_tip=(*curr_tree)[curid-1];
  new_tip->set_tip(TRUE);
  treein>>blen>>dump;


  treein>>con_int[0]>>dump;
  treein>>con_int[1];
  
  
  new_tip->set_brlen_ci(FALSE, con_int[0]);
  new_tip->set_brlen_ci(TRUE, con_int[1]);

  new_tip->set_expect_subs_site(blen);
  cout<<"Set: "<<blen<<"\t"<<con_int[0]<<"\t"<<con_int[1]<<endl<<flush;
  
  return(new_tip);

}  //End Read_PAUP_Tree::get_tip



Branch* Read_PAUP_w_brlen_CI_Tree::get_interior(ifstream &treein)
{
  int i;
  char dump;
  double blen;
  Branch *curbranch;
  
  treein>>rel1;
  treein>>rel1;

  curbranch=initialize_branch();
  curbranch->set_tip(FALSE);
  treein>>blen>>dump;
  treein>>con_int[0]>>dump;
  treein>>con_int[1];

	cout<<"Set interior len: "<<blen<<" u: "<<con_int[1]<<" l: "<<con_int[0]<<endl<<flush;

  curbranch->set_brlen_ci(FALSE, con_int[0]);
  curbranch->set_brlen_ci(TRUE, con_int[1]);
  curbranch->set_expect_subs_site(blen);
  
 
  return(curbranch);
}  //End Read_PAUP_Tree::get_interior


Branch* Read_PAUP_w_brlen_CI_Tree::make_subtree(ifstream &treein)
{
  int i;
  double brn;
  char dump;
  Branch *sibling1, *sibling2, *parent;
  
  if (rel=='(')
    nest_level++;
  
  treein>>rel;
  
  if (rel=='(')
    sibling1=make_subtree(treein);                                     //Recursion
  else
    sibling1=get_tip(treein);

  treein>>rel;
  treein>>rel;
 
  if (rel=='(')
    sibling2=make_subtree(treein);                                 //Recursion
  else
    sibling2=get_tip(treein);

#if defined (DEBUG)
  cout<<" Sibling 1: "<<sibling1->get_name()<<" Sibling 2: "<<sibling2->get_name()<<" Nest level: "<<nest_level<<" ";
#endif 

  curr_tree->set_as_siblings(sibling1, sibling2);
  
  if (nest_level>1)
    parent=get_interior(treein);
  else
    {
      parent=initialize_branch();
      parent->set_tip(FALSE);
	  if (curr_tree->rooted_tree() == FALSE) {
		parent->set_brnlen(0.0);
		parent->set_expect_subs_site(0.0);      
	  }
	  else {
		treein>>rel;
		treein>>rel;
		treein>>brn>>dump;
		parent->set_expect_subs_site(brn);
		treein>>con_int[0]>>dump;
		treein>>con_int[1];

		cout<<"Set Root len: "<<brn<<" u: "<<con_int[1]<<" l: "<<con_int[0]<<endl<<flush;

		parent->set_brlen_ci(FALSE, con_int[0]);
		parent->set_brlen_ci(TRUE, con_int[1]);
	  }
      parent->set_parent(0);
      parent->set_name("Root");

 }
  curr_tree->set_as_parent_child(parent, sibling1, 0);
  curr_tree->set_as_parent_child(parent, sibling2, 1);
 
  nest_level--;
  return (parent);
}  //End Read_PAUP_Tree::make_subtree

