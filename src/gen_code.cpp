#include <iostream>
#include <iomanip>
#include "gen_code.h"
#include <math.h>

using namespace::std;

//Genetic_code Public functions
Genetic_code::Genetic_code ()
{
}



Genetic_code::Genetic_code(Exchange *cexchange, BOOL standard_code, const char *code_file)
{
 
  aa_diffs=0;
  cys_cc_zero=FALSE;
  curr_exchange=cexchange;
  if(standard_code==FALSE) {
      get_external_code(code_file);
      cout<<"Non-standard genetic code in use\n";
    } 
 else {
     //cout<<"Using standard genetic code\n";
     non_stops=61;
     num_stops=3;
    make_code();
   }
  

 aa_diffs=0;
 setup_aa_diffs();

}  //End Genetic_code::Genetic_code


Genetic_code::Genetic_code(Exchange *cexchange)
{
    
    aa_diffs=0;
    cys_cc_zero=FALSE;
    curr_exchange=cexchange;
    
    //cout<<"Using standard genetic code\n";
    non_stops=61;
    num_stops=3;
    make_code();
   
    
    aa_diffs=0;
    setup_aa_diffs();
    
}  //End Genetic_code::Genetic_code



int Genetic_code::get_pos_n (int codon_num, int position)
{
  switch (position)
    {
    case 0:
       return((int)(codon_num/16));
       break;
    case 1: 
      return((int)((codon_num%16)/4));
      break;
    case 2:
      return((int)((codon_num%16)%4));
      break;
    default:
        return(-1);
        break;
    }
} //End Genetic_code::get_pos_n




int Genetic_code::get_codon_num(int *codon)
{
  int i, j;
  BOOL is_ambig=FALSE;

  for(i=0; i<3; i++) 
    is_ambig=base_is_ambig(codon[i]);
  
  if (is_ambig == FALSE)
    return((codon[0]*16)+(codon[1]*4)+codon[2]);
  else
    {
      cerr<<"ERROR: Ambiguity code in input DNA sequence.  Replacing "
	  <<num_to_base(codon[0])<<num_to_base(codon[1])<<num_to_base(codon[2])<<" with ";
      for(i=0; i<3; i++) {
	j=0;
	while(valid_ambig_spec(codon[i], j)==FALSE) j++;
	codon[i]=j;
      }
      cerr<<num_to_base(codon[0])<<num_to_base(codon[1])<<num_to_base(codon[2])<<"\n";
      return(get_codon_num(codon));
    }
}  //End Genetic_code::get_codon_num


Sequence_dataset * Genetic_code::make_codon_sequence (Sequence_dataset *nucleotide_seqs)
{
  int taxa, i, j, temp_codon[3], num_codons;
  Sequence_dataset *new_data;

  if (curr_exchange->get_num_sites()==(*nucleotide_seqs)[0].Sequence_size())
      num_codons=curr_exchange->get_num_codons();
     else
       num_codons=(*nucleotide_seqs)[0].Sequence_size()/3;
  new_data = new Sequence_dataset (curr_exchange->get_num_taxa(), num_codons);



  for(taxa=0; taxa<curr_exchange->get_num_taxa(); taxa++)
    {
      (*new_data)[taxa].Assign_name((*nucleotide_seqs)[taxa].Sequence_name());
      for(i=0; i<num_codons; i++)
	{
	  for(j=0; j<3; j++)
	    temp_codon[j]=(*nucleotide_seqs)[taxa][i*3+j];
	  if (temp_codon[0]<=3 && temp_codon[1]<=3 && temp_codon[2]<=3)
	    (*new_data)[taxa].Assign_site(i, get_codon_num(temp_codon));
	  else
	    (*new_data)[taxa].Assign_site(i, 64);
	}
    } 
  return(new_data);
}

Sequence_dataset * Genetic_code::make_codon_sequence_arbitrary_size (Sequence_dataset *nucleotide_seqs)
{
    int taxa, i, j, temp_codon[3], *num_codons;
    Sequence_dataset *new_data;
    
    num_codons=new int [curr_exchange->get_num_taxa()];
    for(taxa=0; taxa<curr_exchange->get_num_taxa(); taxa++)
        num_codons[taxa]=(*nucleotide_seqs)[taxa].Sequence_size()/3;
    
    new_data = new Sequence_dataset (curr_exchange->get_num_taxa(), num_codons, CODON);
    
   // cout<<"Created dataset with arbitrary sequence sizes\n";
    
    for(taxa=0; taxa<curr_exchange->get_num_taxa(); taxa++)
    {
        (*new_data)[taxa].Assign_name((*nucleotide_seqs)[taxa].Sequence_name());
        for(i=0; i<num_codons[taxa]; i++)
        {
            for(j=0; j<3; j++)
                temp_codon[j]=(*nucleotide_seqs)[taxa][i*3+j];
            if (temp_codon[0]<=3 && temp_codon[1]<=3 && temp_codon[2]<=3)
                (*new_data)[taxa].Assign_site(i, get_codon_num(temp_codon));
            else
                (*new_data)[taxa].Assign_site(i, 64);
        }
    }
    delete[] num_codons;
    return(new_data);
}


BOOL Genetic_code::is_non_synon(int *codon1, int *codon2)
  //This function lets the likeluhood function decide if a non-syon 
  //codon subsitution has occured.  The probablity of a non-syn subsition 
  //is raised to the power of this function.  If the sub is syn, the function returns 0 
  //and no probablity is added
  
{
  if (gen_code[codon1[0]][codon1[1]][codon1[2]]==
      gen_code[codon2[0]][codon2[1]][codon2[2]])
    return(FALSE);
  else
    return(TRUE);
}  //end Genetic_code::is_non_synon





BOOL Genetic_code::is_stop (int codon_num)
{
  if (gen_code[get_pos_n(codon_num, 0)]
      [get_pos_n(codon_num, 1)][get_pos_n(codon_num, 2)]=='*')
    return(TRUE);
  else
    return(FALSE);
}  //End Genetic_code::is_stop



int Genetic_code::one_diff(int codon_num, int position_num, int iteration)
{
  int factor;

  switch (position_num)
    {
    case 0:
      factor=16;
      break;
    case 1: 
      factor=4;
      break;
    case 2:
      factor=1;
      break;
    }

  
  switch (Genetic_code::get_pos_n(codon_num, position_num))
    {
    case 0:
      switch (iteration) {
        case 0:
          return(codon_num+factor);
          break;
        case 1:
          return(codon_num+2*factor);
          break;
        case 2:
          return(codon_num+3*factor);
          break;
        default:
          return(-1);
          break;
	}
    case 1: 
      switch (iteration) {
        case 0:
          return(codon_num-factor);
          break;
        case 1:
          return(codon_num+factor);
          break;
        case 2:
          return(codon_num+2*factor);
          break;
        default:
          return(-1);
          break;
	}
      
    case 2:
      switch (iteration) {
        case 0:
          return(codon_num-(2*factor));
          break;
        case 1:
          return(codon_num-factor);
          break;
        case 2:
          return(codon_num+factor);
          break;
        default:
          return(-1);
          break;
	}
      
    case 3:
      switch (iteration) {
        case 0:
          return(codon_num-(3*factor));
          break;
        case 1:
          return(codon_num-(2*factor));
          break;
        case 2:
          return(codon_num-factor);
          break;
        default:
          return(-1);
          break;
        }
    default:
        return(-1);
        break;
    }
  
}  //End Genetic_code::one_diff




int Genetic_code::multiple_subs(int *codon1, int *codon2)
{
  int i, count=0;

  for(i=0; i<3; i++)
    if(codon1[i]!=codon2[i])
      count++;

  return(count);
} //End Genetic_code::multiple_subs




void Genetic_code::diff_pos(int *codon1, int *codon2, BOOL is_diff[3])
{
  int i, diff1,diff2;
  for(i=0; i<3; i++)
    {
      if(codon1[i]!=codon2[i])
	is_diff[i]=TRUE;
      else
	is_diff[i]=FALSE;
    }
}  //End Genetic_code::diff_pos





double Genetic_code::synon_paths(int num_subs, int *codon1, int *codon2)
{
  int i, j, inter_codon[3], inter_codon2[3], first_c, second_c, third_c, diff1, diff2, passes=0;
  double paths=0;
  BOOL diffs[3];

  switch(num_subs)
    {
    case 1:
      if (is_non_synon(codon1, codon2)!=TRUE)
	paths+=1;
      break;
    

    case 2:
      diff_pos(codon1, codon2, diffs);

      diff1=0; 
      while(diffs[diff1]==FALSE)
	diff1++; 

      diff2=diff1+1; 
      while(diffs[diff2]==FALSE)
	diff2++;      

      for(i=0; i<diff1; i++)
	inter_codon[i]=codon1[i];
      
      inter_codon[diff1]=codon2[diff1];
      
      for(i=diff1+1; i<3; i++)
	inter_codon[i]=codon1[i];
	      
      if(is_stop(get_codon_num(inter_codon))==FALSE)
	{
	  passes++;
	  if (is_non_synon(codon1, inter_codon)!=TRUE)
	    paths+=1;
	  if (is_non_synon(inter_codon, codon2)!=TRUE)
	    paths+=1;
	}     
  
      for(i=0; i<diff2; i++)
	inter_codon[i]=codon1[i];
      
      inter_codon[diff2]=codon2[diff2];
      
      for(i=diff2+1; i<3; i++)
	inter_codon[i]=codon1[i];
	      
      if(is_stop(get_codon_num(inter_codon))==FALSE)
	{
	  passes++;
	  if (is_non_synon(codon1, inter_codon)!=TRUE)
	    paths+=1;
	  if (is_non_synon(inter_codon, codon2)!=TRUE)
	    paths+=1;
	}    

      paths=paths/passes;
      break;
      

    case 3:
      for(first_c=0; first_c<3; first_c++)
	for(second_c=0; second_c<3; second_c++)
	  for(third_c=0; third_c<3; third_c++)
	    {
	      if (first_c!=second_c && first_c!=third_c && second_c!=third_c)
		{
		  for(i=0; i<3; i++)
		    inter_codon[i]=codon1[i];
		  
		  inter_codon[first_c]=codon2[first_c];

		  for(i=0; i<3; i++)
		    inter_codon2[i]=inter_codon[i];
		  
		  inter_codon2[second_c]=codon2[second_c];

		  if (is_stop(get_codon_num(inter_codon))==FALSE &&
		      is_stop(get_codon_num(inter_codon2))==FALSE)
		    {
		      passes++;
		      if (is_non_synon(codon1, inter_codon)!=TRUE)
			paths+=1;
		      if (is_non_synon(inter_codon, inter_codon2)!=TRUE)
			paths+=1;
		      if (is_non_synon(inter_codon2, codon2)!=TRUE)
			paths+=1;
		    }
		}
	    }
      paths=paths/passes;
      break;
    }
  return(paths);
  
} //End Genetic_code::synon_paths;




double Genetic_code::non_synon_paths(int num_subs, int *codon1, int *codon2)
{
  int i, diff1, diff2, first_c, second_c, third_c, inter_codon[3], inter_codon2[3], passes=0;
  double paths=0.0;
  BOOL diffs[3];
  
  switch(num_subs)
    {
    case 1:
      if (is_non_synon(codon1, codon2)==TRUE)
	paths+=1;
      break;
    

    case 2:
      diff_pos(codon1, codon2, diffs);

      diff1=0; 
      while(diffs[diff1]==FALSE)
	diff1++; 

      diff2=diff1+1; 
      while(diffs[diff2]==FALSE)
	diff2++;      

      for(i=0; i<diff1; i++)
	inter_codon[i]=codon1[i];
      
      inter_codon[diff1]=codon2[diff1];
      
      for(i=diff1+1; i<3; i++)
	inter_codon[i]=codon1[i];
	      
      if(is_stop(get_codon_num(inter_codon))==FALSE)
	{
	  passes++;
	  if (is_non_synon(codon1, inter_codon)==TRUE)
	    paths+=1;
	  if (is_non_synon(inter_codon, codon2)==TRUE)
	    paths+=1;
	}     
  
      for(i=0; i<diff2; i++)
	inter_codon[i]=codon1[i];
      
      inter_codon[diff2]=codon2[diff2];
      
      for(i=diff2+1; i<3; i++)
	inter_codon[i]=codon1[i];
	      
      if(is_stop(get_codon_num(inter_codon))==FALSE)
	{
	  passes++;
	  if (is_non_synon(codon1, inter_codon)==TRUE)
	    paths+=1;
	  if (is_non_synon(inter_codon, codon2)==TRUE)
	    paths+=1;
	}    

      paths=paths/passes;
      break;
      

    case 3:
      for(first_c=0; first_c<3; first_c++)
	for(second_c=0; second_c<3; second_c++)
	  for(third_c=0; third_c<3; third_c++)
	    {
	      if (first_c!=second_c && first_c!=third_c && second_c!=third_c)
		{
		  for(i=0; i<3; i++)
		    inter_codon[i]=codon1[i];
		  
		  inter_codon[first_c]=codon2[first_c];

		  for(i=0; i<3; i++)
		    inter_codon2[i]=inter_codon[i];
		  
		  inter_codon2[second_c]=codon2[second_c];

		  if (is_stop(get_codon_num(inter_codon))==FALSE &&
		      is_stop(get_codon_num(inter_codon2))==FALSE)
		    {
		      passes++;
		      if (is_non_synon(codon1, inter_codon)==TRUE)
			paths+=1;
		      if (is_non_synon(inter_codon, inter_codon2)==TRUE)
			paths+=1;
		      if (is_non_synon(inter_codon2, codon2)==TRUE)
			paths+=1;
		    }
		}
	    }
      paths=paths/passes;
      break;
    }
  return(paths);
 
} //End Genetic_code::non_synon_paths




void Genetic_code::count_by_degen(int *codon1, int *codon2, int diff_pos,
				  double &trs_paths, double &trv_paths)
{
  switch (Li_93_degen(codon1, diff_pos))
    {
    case 0:
      if (is_transition(codon1[diff_pos], codon2[diff_pos])==TRUE)
	++trs_paths;
      else
	++trv_paths; 
      break;

    case 2:
      if (is_non_synon(codon1, codon2)==TRUE)
	++trv_paths;
      if  (is_non_synon(codon1, codon2)==FALSE)
	++trs_paths;
      break;

    case 4:
      if (is_transition(codon1[diff_pos], codon2[diff_pos])==TRUE)
	++trs_paths;
      else
	++trv_paths; 
      break;
    }
  
} //End Genetic_code::count_by_degen

 


void Genetic_code::Comeron_levels(int *codon1, int *codon2, int diff_pos, 
				  double trs_paths[5], double trv_paths[5], double divisor, int degen)
{
  if (is_transition(codon1[diff_pos], codon2[diff_pos])==TRUE)
    trs_paths[degen]+=1/divisor;
  else
    trv_paths[degen]+=1/divisor;  
} //End Genetic_code::count_by_degen

 


int Genetic_code::Li_93_degen(int codon[3], int pos)
{
  int i, num_so_far=0, new_codon[3];

  for(i=0; i<3; i++)
    new_codon[i]=codon[i];

  for(i=0; i<codon[pos]; i++)
    {
      new_codon[pos]=i;
      if (is_non_synon(codon, new_codon)==FALSE)
	num_so_far++;
    }

  for(i=codon[pos]+1; i<4; i++)
    {
      new_codon[pos]=i;
      if (is_non_synon(codon, new_codon)==FALSE)
	num_so_far++;
    }
   

  switch (num_so_far)
    {
    case 0:
      return(0);
      break;
    case 1:
      return(2);
      break;
    case 2:
      return(2);
      break;
    case 3:
      return(4);
        break;
    default:
        return(-1);
        break;
      
    }
  
}  //End Li_93_degen




int Genetic_code::Comeron_95_degen(int codon[3], int pos)
{
  int i, s_so_far=0, v_so_far=0, new_codon[3], retval;

  for(i=0; i<3; i++)
    new_codon[i]=codon[i];

    for(i=0; i<codon[pos]; i++) {
        new_codon[pos]=i;

        if (is_non_synon(codon, new_codon)==FALSE) {
            if (is_transition(codon[pos], new_codon[pos])==TRUE)
              s_so_far++;
            else
              v_so_far++;
        }
    }

  for(i=codon[pos]+1; i<4; i++) {
      new_codon[pos]=i;
      if (is_non_synon(codon, new_codon)==FALSE) {
        if (is_transition(codon[pos], new_codon[pos])==TRUE)
          s_so_far++;
        else
          v_so_far++;
      }
    }
   

  if (v_so_far+s_so_far==0)
    retval=0;
  else if (v_so_far==2 && s_so_far==1)
    retval=4;
  else if (s_so_far==1 && (v_so_far==1 || v_so_far==0))
    retval=2;
  else if ((v_so_far==2 || v_so_far==1) && s_so_far==0)
    retval=1;
  
  return(retval);
  
}  //End Comeron_95_degen




double Genetic_code::possible_synon_subs(int *codon)
{
  int i,j, num_so_far=0, new_codon[3];

  for(j=0; j<3; j++)
    new_codon[j]=codon[j];

  for(i=0; i<3; i++)
    {
      for(j=0; j<codon[i]; j++)
	{
	  new_codon[i]=j;
	  if (is_non_synon(codon, new_codon)==FALSE)
	      num_so_far++;
	 
	}

      for(j=codon[i]+1; j<4; j++)
	{
	  new_codon[i]=j;
	  if (is_non_synon(codon, new_codon)==FALSE)
	    num_so_far++;   
	}

      new_codon[i]=codon[i];
     
    }
  return((1.0*num_so_far)/3.0);
} //End Genetic_code::possible_synon_subs




double Genetic_code::possible_non_synon_subs(int *codon)
{
  int i,j, num_so_far=0, new_codon[3];

  for(j=0; j<3; j++)
    new_codon[j]=codon[j];

  for(i=0; i<3; i++)
    {
      for(j=0; j<codon[i]; j++)
	{
	  new_codon[i]=j;
	  if (is_non_synon(codon, new_codon)==TRUE && 
	      is_stop(get_codon_num(new_codon))==FALSE)
	    num_so_far++;
	}

      for(j=codon[i]+1; j<4; j++)
	{
	  new_codon[i]=j;
	  if (is_non_synon(codon, new_codon)==TRUE && 
	      is_stop(get_codon_num(new_codon))==FALSE)
	    num_so_far++;
	}
      new_codon[i]=codon[i];

    }
  // cout<<"Codon: "<<codon[0]<<codon[1]<<codon[2]<<" Poss synon subs:" <<(1.0*num_so_far)/3<<endl;
  return((1.0*num_so_far)/3.0);
  
} //End Genetic_code::possible_non_synon_subs




char Genetic_code::return_AA(int *codon)
{
  int num;
  num=get_codon_num(codon);
  return(gen_code[get_pos_n(num, 0)][get_pos_n(num, 1)][get_pos_n(num, 2)]);
}

char Genetic_code::return_AA(int codon_num)
{
  return(gen_code[get_pos_n(codon_num, 0)][get_pos_n(codon_num, 1)][get_pos_n(codon_num, 2)]);
}


double Genetic_code::aa_pair_diff(int codon1[3], int codon2[3], AA_PROPERTIES prop)
{
  if (is_non_synon(codon1, codon2)==FALSE)
    return(0);
  else     
    return(aa_diffs[readchar_to_aa(return_AA(codon1))][readchar_to_aa(return_AA(codon2))]
	   [curr_exchange->get_prop_index_num(prop)]);
} //End Genetic_code::aa_pair_diff



double Genetic_code::aa_pair_diff(int aa1, int aa2, AA_PROPERTIES prop)
{
   return(aa_diffs[aa1][aa2]
	   [curr_exchange->get_prop_index_num(prop)]);


}


double Genetic_code::aa_pair_diff(int codon1[3], int codon2[3], int prop)
{
  if (is_non_synon(codon1, codon2)==FALSE)
    return(0);
  else     
    return(aa_diffs[readchar_to_aa(return_AA(codon1))][readchar_to_aa(return_AA(codon2))]
	   [prop]);
} //End Genetic_code::aa_pair_diff



double Genetic_code::aa_pair_diff(int aa1, int aa2, int prop)
{
   return(aa_diffs[aa1][aa2][prop]);


}

int Genetic_code::get_stop_num(int stop)
{
  if (stop<num_stops)
    return(stop_pos[stop]);
  else
    {
      cerr<<"Invalid stop num request\n";
      return(-1);
    }
}


void Genetic_code::use_zero_cys_cc()           
{
  cys_cc_zero=TRUE;
  setup_aa_diffs();
}


Genetic_code::~Genetic_code()
{
	int i, j;
	
	if (aa_diffs != 0) {
		for(i=0; i<20; i++) {
			for(j=0; j<20; j++)
				delete[] aa_diffs[i][j];
			delete[] aa_diffs[i];
		}
		delete[] aa_diffs;
	}
}


//End Genetic_code public_functions





//Genetic_code private functions
void Genetic_code::make_code()
  //The standard genetic code is hard coded into this function so that no extraneous files are
  //needed for the executable

{
  int i;

  gen_code[0][0][3]=gen_code[0][0][1]='N';
  gen_code[0][0][0]= gen_code[0][0][2]='K';
  
  for (i=0; i<4; i++)
    gen_code[0][1][i]='T';
  
  gen_code[0][2][1]=gen_code[0][2][3]='S';
  gen_code[0][2][0]=gen_code[0][2][2]='R';
  
  gen_code[0][3][0]=gen_code[0][3][1]=gen_code[0][3][3]='I';
  
  gen_code[0][3][2]='M';
  
  gen_code[1][0][3]=gen_code[1][0][1]='H';
  gen_code[1][0][0]= gen_code[1][0][2]='Q';
  
  for (i=0; i<4; i++)
    gen_code[1][1][i]='P';

  for (i=0; i<4; i++)
    gen_code[1][2][i]='R';

  for (i=0; i<4; i++)
    gen_code[1][3][i]='L';

  gen_code[2][0][3]=gen_code[2][0][1]='D';
  gen_code[2][0][0]= gen_code[2][0][2]='E';

  for (i=0; i<4; i++)
    gen_code[2][1][i]='A';

  for (i=0; i<4; i++)
    gen_code[2][2][i]='G';

  for (i=0; i<4; i++)
    gen_code[2][3][i]='V';

  gen_code[3][0][3]=gen_code[3][0][1]='Y';
  gen_code[3][0][0]= gen_code[3][0][2]='*';
  

  for (i=0; i<4; i++)
    gen_code[3][1][i]='S';

  gen_code[3][2][3]=gen_code[3][2][1]='C';
  gen_code[3][2][0]='*'; 
  gen_code[3][2][2]='W';

  gen_code[3][3][3]=gen_code[3][3][1]='F';
  gen_code[3][3][0]= gen_code[3][3][2]='L';
  
 
  locate_stops();


} //End Genetic_code::make_code




void Genetic_code::get_external_code(const char *filename)
{
  //Code stud for eventual ablity to read external genetic codes
  cerr<<"Module for external codes is not available\n";
}  //End Genetic_code::get_external_code



void Genetic_code::setup_aa_diffs()
{
  int i, j;
  
  if (aa_diffs == 0) {
  
	  aa_diffs = new double ** [20];
  
  
	  for(i=0; i<20; i++) {
		  aa_diffs[i]=new double * [20];
		  for(j=0; j<20; j++)
			  aa_diffs[i][j]=new double [curr_exchange->get_num_aa_props()];
	 }
  }


  for(i=0; i<20; i++)
    {
      for(j=0; j<20; j++)
	{
	  if(i==j)
	    {
	      aa_diffs[i][j][curr_exchange->get_prop_index_num(CHEM_COMP)]=0;
	      aa_diffs[i][j][curr_exchange->get_prop_index_num(POLARITY)]=0;
	      aa_diffs[i][j][curr_exchange->get_prop_index_num(VOLUME)]=0;
	      aa_diffs[i][j][curr_exchange->get_prop_index_num(ISO_ELEC)]=0;
		  aa_diffs[i][j][curr_exchange->get_prop_index_num(HYDROPATHY)]=0;
	    }
	  else
	    {
	      aa_diffs[i][j][curr_exchange->get_prop_index_num(CHEM_COMP)]=
		fabs(individual_aa_props(num_to_aa(i), CHEM_COMP)-
					     individual_aa_props(num_to_aa(j), CHEM_COMP));
	      aa_diffs[i][j][curr_exchange->get_prop_index_num(POLARITY)]=
		fabs(individual_aa_props(num_to_aa(i), POLARITY)-
					    individual_aa_props(num_to_aa(j), POLARITY));
	      aa_diffs[i][j][curr_exchange->get_prop_index_num(VOLUME)]=
		fabs(individual_aa_props(num_to_aa(i), VOLUME)-
					  individual_aa_props(num_to_aa(j), VOLUME));
	      aa_diffs[i][j][curr_exchange->get_prop_index_num(ISO_ELEC)]=
		fabs(individual_aa_props(num_to_aa(i), ISO_ELEC)-
					    individual_aa_props(num_to_aa(j), ISO_ELEC));
		  aa_diffs[i][j][curr_exchange->get_prop_index_num(HYDROPATHY)]=
		fabs(individual_aa_props(num_to_aa(i), HYDROPATHY)-
					    individual_aa_props(num_to_aa(j), HYDROPATHY));
	    }
	}
    }

}  //End Genetic_code::setup_aa_diffs



double Genetic_code::individual_aa_props(char aa, AA_PROPERTIES prop)
{
    double ret_val=0.0;
  switch (aa)
    {
    case 'A':
      switch (prop)
		{
		case CHEM_COMP:
			ret_val=0;
			break;
		case POLARITY:
			ret_val=8.1;
			break;
		case VOLUME:
			ret_val=3.1;
			break;
		case ISO_ELEC:
			ret_val=6;
			break;
		case HYDROPATHY:
			ret_val=1.8;
			break;
        default:
            break;
		}
		break;
    case 'C':
	  switch (prop)
		{
		case CHEM_COMP:
		  if( cys_cc_zero ==FALSE)
		    ret_val=2.75;
		  else
		    ret_val=0.5089;
		  break;
		case POLARITY:
			ret_val=5.5;
			break;
		case VOLUME:
			ret_val=5.5;
			break;
		case ISO_ELEC:
			ret_val=5.07;
			break;
		case HYDROPATHY:
			ret_val=2.5;
			break;
        default:
            break;
		}
		break;
	case 'D':
      switch (prop)
		{
		case CHEM_COMP:
			ret_val=1.38;
			break;
		case POLARITY:
			ret_val=13;
			break;
		case VOLUME:
			ret_val=5.4;
			break;
		case ISO_ELEC:
			ret_val=2.77;
			break;
		case HYDROPATHY:
			ret_val=-3.5;
			break;
        default:
            break;
		}
      break;
    case 'E':
      switch (prop)
		{
		case CHEM_COMP:
			ret_val=0.92;
			break;
		case POLARITY:
			ret_val=12.3;
			break;
		case VOLUME:
			ret_val=8.3;
			break;
		case ISO_ELEC:
			ret_val=3.22;
			break;
		case HYDROPATHY:
			ret_val=-3.5;
			break;
        default:
            break;
		}
      break;
    case 'F':
      switch (prop)
		{
		case CHEM_COMP:
            ret_val=0;
			break;
		case POLARITY:
			ret_val=5.2;
			break;
		case VOLUME:
			ret_val=13.2;
			break;
		case ISO_ELEC:
			ret_val=5.48;
			break;
		case HYDROPATHY:
			ret_val=2.8;
			break;
        default:
            break;
		}
      break;
    case 'G':
      switch (prop)
		{
		case CHEM_COMP:
			ret_val=0.74;
			break;
		case POLARITY:
			ret_val=9;
			break;
		case VOLUME:
			ret_val=0.3;
			break;
		case ISO_ELEC:
			ret_val=5.97;
			break;
		case HYDROPATHY:
			ret_val=-0.4;
			break;
        default:
            break;
		}
      break;
    case 'H':
      switch (prop) {
        case CHEM_COMP:
          ret_val=0.58;
          break;
        case POLARITY:
          ret_val=10.4;
          break;
        case VOLUME:
          ret_val=9.6;
          break;
        case ISO_ELEC:
          ret_val=7.59;
          break;
          case HYDROPATHY:
            ret_val=-3.2;
            break;
        default:
          break;
	}
      break;
    case 'I':
      switch (prop) {
        case CHEM_COMP:
          ret_val=0;
          break;
        case POLARITY:
          ret_val=5.2;
          break;
        case VOLUME:
          ret_val=11.1;
          break;
        case ISO_ELEC:
          ret_val=6.02;
          break;
          case HYDROPATHY:
            ret_val=4.5;
            break;
        default:
          break;
	}
      break;
    case 'K':
      switch (prop){
        case CHEM_COMP:
          ret_val=0.33;
          break;
        case POLARITY:
          ret_val=11.3;
          break;
        case VOLUME:
          ret_val=11.9;
          break;
        case ISO_ELEC:
          ret_val=9.74;
          break;
          case HYDROPATHY:
            ret_val=-3.9;
            break;
        default:
          break;
	}
      break;
    case 'L':
      switch (prop) {
        case CHEM_COMP:
          ret_val=0;
          break;
        case POLARITY:
          ret_val=4.9;
          break;
        case VOLUME:
          ret_val=11.1;
          break;
        case ISO_ELEC:
          ret_val=5.98;
          break;
          case HYDROPATHY:
            ret_val=3.8;
            break;
        default:
          break;
	}
      break;
    case 'M':
      switch (prop) {
        case CHEM_COMP:
          ret_val=0;
          break;
        case POLARITY:
          ret_val=5.7;
          break;
        case VOLUME:
          ret_val=10.5;
          break;
        case ISO_ELEC:
          ret_val=5.74;
          break;
          case HYDROPATHY:
            ret_val=1.9;
            break;
        default:
          break;
	}
      break;
    case 'N':
      switch (prop) {
        case CHEM_COMP:
          ret_val=1.33;
          break;
        case POLARITY:
          ret_val=11.6;
          break;
        case VOLUME:
          ret_val=5.6;
          break;
        case ISO_ELEC:
          ret_val=5.41;
          break;
          case HYDROPATHY:
            ret_val=-3.5;
            break;
          default:
          break;
	}
      break;
    case 'P':
      switch (prop){
        case CHEM_COMP:
          ret_val=0.39;
          break;
        case POLARITY:
          ret_val=8;
          break;
        case VOLUME:
          ret_val=3.25;
          break;
        case ISO_ELEC:
          ret_val=6.3;
          break;
        case HYDROPATHY:
            ret_val=-1.6;
            break;
        default:
          break;
	}
      break;
    case 'Q':
      switch (prop){
        case CHEM_COMP:
          ret_val=0.89;
          break;
        case POLARITY:
          ret_val=10.5;
          break;
        case VOLUME:
          ret_val=8.5;
          break;
        case ISO_ELEC:
          ret_val=5.65;
          break;
          case HYDROPATHY:
            ret_val=-3.5;
            break;
        default:
          break;
	}
      break;
    case 'R':
      switch (prop){
        case CHEM_COMP:
          ret_val=0.65;
          break;
        case POLARITY:
          ret_val=10.5;
          break;
        case VOLUME:
          ret_val=12.4;
          break;
        case ISO_ELEC:
          ret_val=10.76;
          break;
          case HYDROPATHY:
            ret_val=-4.5;
            break;
        default:
          break;
	}
      break;
    case 'S':
      switch (prop){
        case CHEM_COMP:
          ret_val=1.42;
          break;
        case POLARITY:
          ret_val=9.2;
          break;
        case VOLUME:
          ret_val=3.2;
          break;
        case ISO_ELEC:
          ret_val=5.68;
          break;
          case HYDROPATHY:
            ret_val=-0.8;
            break;
        default:
          break;
	}
      break;
    case 'T':
      switch (prop){
        case CHEM_COMP:
          ret_val=0.71;
          break;
        case POLARITY:
          ret_val=8.6;
          break;
        case VOLUME:
          ret_val=6.1;
          break;
        case ISO_ELEC:
          ret_val=6.16;
          break;
        case HYDROPATHY:
            ret_val=-0.7;
            break;
        default:
          break;
	}
      break;
    case 'V':
      switch (prop){
        case CHEM_COMP:
          ret_val=0;
          break;
        case POLARITY:
          ret_val=5.9;
          break;
        case VOLUME:
          ret_val=8.4;
          break;
        case ISO_ELEC:
          ret_val=5.96;
          break;
          case HYDROPATHY:
            ret_val=4.2;
            break;
          default:
            break;
	}
      break;
    case 'W':
      switch (prop){
        case CHEM_COMP:
          ret_val=0.13;
          break;
        case POLARITY:
          ret_val=5.4;
          break;
        case VOLUME:
          ret_val=17.0;
          break;
        case ISO_ELEC:
          ret_val=5.89;
          break;
        case HYDROPATHY:
            ret_val=-0.9;
            break;
        default:
          break;
	}
      break;
    case 'Y':
      switch (prop) {
        case CHEM_COMP:
          ret_val=0.2;
          break;
        case POLARITY:
          ret_val=6.2;
          break;
        case VOLUME:
          ret_val=13.6;
          break;
        case ISO_ELEC:
          ret_val=5.66;
          break;
          case HYDROPATHY:
            ret_val=-1.3;
            break;
        default:
          break;
	}
      break;

    default:
      ret_val=0;
      break;

    }
    return(ret_val);
}



void Genetic_code::locate_stops()
{
  int i, sp=0;
  for (i=0; i<64; i++) {
    if (is_stop(i) == TRUE) 
      {
	stop_pos[sp]=i;
	sp++;
      }
  }
}








//Genetic_code private functions
Vert_mito_genetic_code::Vert_mito_genetic_code(Exchange *cexchange) :
  Genetic_code()
{
    curr_exchange=cexchange;

    make_code();
    non_stops=62;
    num_stops=2;

    aa_diffs=0;
    setup_aa_diffs();

}  //End Genetic_code::Genetic_code


void Vert_mito_genetic_code::make_code()
  //The standard genetic code is hard coded into this function so that no extraneous files are
  //needed for the executable

{
  int i;

  gen_code[0][0][3]=gen_code[0][0][1]='N';
  gen_code[0][0][0]= gen_code[0][0][2]='K';
  
  for (i=0; i<4; i++)
    gen_code[0][1][i]='T';
  
  gen_code[0][2][1]=gen_code[0][2][3]='S';
  gen_code[0][2][0]=gen_code[0][2][2]='*';
  
  gen_code[0][3][0]='M';
  gen_code[0][3][1]=gen_code[0][3][3]='I';
  
  gen_code[0][3][2]='M';
  
  gen_code[1][0][3]=gen_code[1][0][1]='H';
  gen_code[1][0][0]= gen_code[1][0][2]='Q';
  
  for (i=0; i<4; i++)
    gen_code[1][1][i]='P';

  for (i=0; i<4; i++)
    gen_code[1][2][i]='R';

  for (i=0; i<4; i++)
    gen_code[1][3][i]='L';

  gen_code[2][0][3]=gen_code[2][0][1]='D';
  gen_code[2][0][0]= gen_code[2][0][2]='E';

  for (i=0; i<4; i++)
    gen_code[2][1][i]='A';

  for (i=0; i<4; i++)
    gen_code[2][2][i]='G';

  for (i=0; i<4; i++)
    gen_code[2][3][i]='V';

  gen_code[3][0][3]=gen_code[3][0][1]='Y';
  gen_code[3][0][0]= gen_code[3][0][2]='*';

  for (i=0; i<4; i++)
    gen_code[3][1][i]='S';

  gen_code[3][2][3]=gen_code[3][2][1]='C';
  gen_code[3][2][0]='W'; 
  gen_code[3][2][2]='W';

  gen_code[3][3][3]=gen_code[3][3][1]='F';
  gen_code[3][3][0]= gen_code[3][3][2]='L';

  locate_stops();
  
} //End Vert_mito_genetic_code::make_code


Yeast_mito_genetic_code::Yeast_mito_genetic_code(Exchange *cexchange) :
  Genetic_code()
{
    curr_exchange=cexchange;

     make_code();
     non_stops=60;
     num_stops=4;

    aa_diffs=0;
    setup_aa_diffs();

}  //End Genetic_code::Genetic_code


void Yeast_mito_genetic_code::make_code()
  //The standard genetic code is hard coded into this function so that no extraneous files are
  //needed for the executable

{
  int i;

  gen_code[0][0][3]=gen_code[0][0][1]='N';
  gen_code[0][0][0]= gen_code[0][0][2]='K';
  
  for (i=0; i<4; i++)
    gen_code[0][1][i]='T';
  
  gen_code[0][2][1]=gen_code[0][2][3]='S';
  gen_code[0][2][0]=gen_code[0][2][2]='R';
  
  gen_code[0][3][0]='M';
  gen_code[0][3][1]=gen_code[0][3][3]='I';
  
  gen_code[0][3][2]='M';
  
  gen_code[1][0][3]=gen_code[1][0][1]='H';
  gen_code[1][0][0]= gen_code[1][0][2]='Q';
  
  for (i=0; i<4; i++)
    gen_code[1][1][i]='P';

  //Absent==#
  gen_code[1][2][0]=gen_code[1][2][1]='*';
  for (i=2; i<4; i++)
    gen_code[1][2][i]='R';

  for (i=0; i<4; i++)
    gen_code[1][3][i]='T';

  gen_code[2][0][3]=gen_code[2][0][1]='D';
  gen_code[2][0][0]= gen_code[2][0][2]='E';

  for (i=0; i<4; i++)
    gen_code[2][1][i]='A';

  for (i=0; i<4; i++)
    gen_code[2][2][i]='G';

  for (i=0; i<4; i++)
    gen_code[2][3][i]='V';

  gen_code[3][0][3]=gen_code[3][0][1]='Y';
  gen_code[3][0][0]= gen_code[3][0][2]='*';

  for (i=0; i<4; i++)
    gen_code[3][1][i]='S';

  gen_code[3][2][3]=gen_code[3][2][1]='C';
  gen_code[3][2][0]='W'; 
  gen_code[3][2][2]='W';

  gen_code[3][3][3]=gen_code[3][3][1]='F';
  gen_code[3][3][0]= gen_code[3][3][2]='L';

  locate_stops();
  
} //End Yeast_mito_genetic_code::make_code



Mold_mito_genetic_code::Mold_mito_genetic_code(Exchange *cexchange) :
  Genetic_code()
{
    curr_exchange=cexchange;

    non_stops=62;
    num_stops=2;
    make_code();

    aa_diffs=0;
    setup_aa_diffs();

}  //End Genetic_code::Genetic_code


void Mold_mito_genetic_code::make_code()
  //The standard genetic code is hard coded into this function so that no extraneous files are
  //needed for the executable

{
  int i;

  gen_code[0][0][3]=gen_code[0][0][1]='N';
  gen_code[0][0][0]= gen_code[0][0][2]='K';
  
  for (i=0; i<4; i++)
    gen_code[0][1][i]='T';
  
  gen_code[0][2][1]=gen_code[0][2][3]='S';
  gen_code[0][2][0]=gen_code[0][2][2]='R';
  
  gen_code[0][3][0]=gen_code[0][3][1]=gen_code[0][3][3]='I';
  
  gen_code[0][3][2]='M';
  
  gen_code[1][0][3]=gen_code[1][0][1]='H';
  gen_code[1][0][0]= gen_code[1][0][2]='Q';
  
  for (i=0; i<4; i++)
    gen_code[1][1][i]='P';

  for (i=0; i<4; i++)
    gen_code[1][2][i]='R';

  for (i=0; i<4; i++)
    gen_code[1][3][i]='L';

  gen_code[2][0][3]=gen_code[2][0][1]='D';
  gen_code[2][0][0]= gen_code[2][0][2]='E';

  for (i=0; i<4; i++)
    gen_code[2][1][i]='A';

  for (i=0; i<4; i++)
    gen_code[2][2][i]='G';

  for (i=0; i<4; i++)
    gen_code[2][3][i]='V';

  gen_code[3][0][3]=gen_code[3][0][1]='Y';
  gen_code[3][0][0]= gen_code[3][0][2]='*';

  for (i=0; i<4; i++)
    gen_code[3][1][i]='S';

  gen_code[3][2][3]=gen_code[3][2][1]='C';
  gen_code[3][2][0]='W'; 
  gen_code[3][2][2]='W';

  gen_code[3][3][3]=gen_code[3][3][1]='F';
  gen_code[3][3][0]= gen_code[3][3][2]='L';

  locate_stops();
  
} //End Mold_mito_genetic_code::make_code



Invert_mito_genetic_code::Invert_mito_genetic_code(Exchange *cexchange) :
  Genetic_code()
{
    curr_exchange=cexchange;

    non_stops=62;
    num_stops=2;
    make_code();

    aa_diffs=0;
    setup_aa_diffs();

}  //End Genetic_code::Genetic_code



void Invert_mito_genetic_code::make_code()
  //The standard genetic code is hard coded into this function so that no extraneous files are
  //needed for the executable

{
  int i;

  gen_code[0][0][3]=gen_code[0][0][1]='N';
  gen_code[0][0][0]= gen_code[0][0][2]='K';
  
  for (i=0; i<4; i++)
    gen_code[0][1][i]='T';
  
  gen_code[0][2][1]=gen_code[0][2][3]='S';
  gen_code[0][2][0]=gen_code[0][2][2]='S';
  
  gen_code[0][3][0]='M';
  gen_code[0][3][1]=gen_code[0][3][3]='I';
  
  gen_code[0][3][2]='M';
  
  gen_code[1][0][3]=gen_code[1][0][1]='H';
  gen_code[1][0][0]= gen_code[1][0][2]='Q';
  
  for (i=0; i<4; i++)
    gen_code[1][1][i]='P';

  for (i=0; i<4; i++)
    gen_code[1][2][i]='R';

  for (i=0; i<4; i++)
    gen_code[1][3][i]='L';

  gen_code[2][0][3]=gen_code[2][0][1]='D';
  gen_code[2][0][0]= gen_code[2][0][2]='E';

  for (i=0; i<4; i++)
    gen_code[2][1][i]='A';

  for (i=0; i<4; i++)
    gen_code[2][2][i]='G';

  for (i=0; i<4; i++)
    gen_code[2][3][i]='V';

  gen_code[3][0][3]=gen_code[3][0][1]='Y';
  gen_code[3][0][0]= gen_code[3][0][2]='*';

  for (i=0; i<4; i++)
    gen_code[3][1][i]='S';

  gen_code[3][2][3]=gen_code[3][2][1]='C';
  gen_code[3][2][0]='W'; 
  gen_code[3][2][2]='W';

  gen_code[3][3][3]=gen_code[3][3][1]='F';
  gen_code[3][3][0]= gen_code[3][3][2]='L';

  locate_stops();
  
} //End Invert_mito_genetic_code::make_code



Ciliate_nuc_genetic_code::Ciliate_nuc_genetic_code(Exchange *cexchange) :
  Genetic_code()
{
    curr_exchange=cexchange;

    make_code();
    num_stops=1;
    non_stops=63;

    aa_diffs=0;
    setup_aa_diffs();

}  //End Genetic_code::Genetic_code


void Ciliate_nuc_genetic_code::make_code()
  //The standard genetic code is hard coded into this function so that no extraneous files are
  //needed for the executable

{
  int i;

  gen_code[0][0][3]=gen_code[0][0][1]='N';
  gen_code[0][0][0]= gen_code[0][0][2]='K';
  
  for (i=0; i<4; i++)
    gen_code[0][1][i]='T';
  
  gen_code[0][2][1]=gen_code[0][2][3]='S';
  gen_code[0][2][0]=gen_code[0][2][2]='R';
  
  gen_code[0][3][0]=gen_code[0][3][1]=gen_code[0][3][3]='I';
  
  gen_code[0][3][2]='M';
  
  gen_code[1][0][3]=gen_code[1][0][1]='H';
  gen_code[1][0][0]= gen_code[1][0][2]='Q';
  
  for (i=0; i<4; i++)
    gen_code[1][1][i]='P';

  for (i=0; i<4; i++)
    gen_code[1][2][i]='R';

  for (i=0; i<4; i++)
    gen_code[1][3][i]='L';

  gen_code[2][0][3]=gen_code[2][0][1]='D';
  gen_code[2][0][0]= gen_code[2][0][2]='E';

  for (i=0; i<4; i++)
    gen_code[2][1][i]='A';

  for (i=0; i<4; i++)
    gen_code[2][2][i]='G';

  for (i=0; i<4; i++)
    gen_code[2][3][i]='V';

  gen_code[3][0][3]=gen_code[3][0][1]='Y';
  gen_code[3][0][0]= gen_code[3][0][2]='Q';

  for (i=0; i<4; i++)
    gen_code[3][1][i]='S';

  gen_code[3][2][3]=gen_code[3][2][1]='C';
  gen_code[3][2][0]='*'; 
  gen_code[3][2][2]='W';

  gen_code[3][3][3]=gen_code[3][3][1]='F';
  gen_code[3][3][0]= gen_code[3][3][2]='L';

  locate_stops();
  
} //End Ciliate_nuc_genetic_code::make_code


Echino_mito_genetic_code::Echino_mito_genetic_code(Exchange *cexchange) :
  Genetic_code()
{
    curr_exchange=cexchange;

    make_code();
    num_stops=2;
    non_stops=62;

    aa_diffs=0;
    setup_aa_diffs();

}  //End Genetic_code::Genetic_code

void Echino_mito_genetic_code::make_code()
  //The standard genetic code is hard coded into this function so that no extraneous files are
  //needed for the executable

{
  int i;

  gen_code[0][0][3]=gen_code[0][0][1]='N';
  gen_code[0][0][0]= 'N';
  gen_code[0][0][2]='K';
  
  for (i=0; i<4; i++)
    gen_code[0][1][i]='T';
  
  gen_code[0][2][1]=gen_code[0][2][3]='S';
  gen_code[0][2][0]=gen_code[0][2][2]='S';
  
  gen_code[0][3][0]=gen_code[0][3][1]=gen_code[0][3][3]='I';
  
  gen_code[0][3][2]='M';
  
  gen_code[1][0][3]=gen_code[1][0][1]='H';
  gen_code[1][0][0]= gen_code[1][0][2]='Q';
  
  for (i=0; i<4; i++)
    gen_code[1][1][i]='P';

  for (i=0; i<4; i++)
    gen_code[1][2][i]='R';

  for (i=0; i<4; i++)
    gen_code[1][3][i]='L';

  gen_code[2][0][3]=gen_code[2][0][1]='D';
  gen_code[2][0][0]= gen_code[2][0][2]='E';

  for (i=0; i<4; i++)
    gen_code[2][1][i]='A';

  for (i=0; i<4; i++)
    gen_code[2][2][i]='G';

  for (i=0; i<4; i++)
    gen_code[2][3][i]='V';

  gen_code[3][0][3]=gen_code[3][0][1]='Y';
  gen_code[3][0][0]= gen_code[3][0][2]='*';

  for (i=0; i<4; i++)
    gen_code[3][1][i]='S';

  gen_code[3][2][3]=gen_code[3][2][1]='C';
  gen_code[3][2][0]='W'; 
  gen_code[3][2][2]='W';

  gen_code[3][3][3]=gen_code[3][3][1]='F';
  gen_code[3][3][0]= gen_code[3][3][2]='L';

  locate_stops();
  
} //End Echino_mito_genetic_code::make_code



Mycoplasma_genetic_code::Mycoplasma_genetic_code(Exchange *cexchange) :
  Genetic_code()
{
    curr_exchange=cexchange;

     make_code();
     num_stops=2;
     non_stops=62;

    aa_diffs=0;
    setup_aa_diffs();

}  //End Mycoplasma_genetic_code::Mycoplasma_genetic_code

void Mycoplasma_genetic_code::make_code()
  //The correct genetic code is hard coded into this function so that no extraneous files are
  //needed for the executable

{
  int i;

  gen_code[0][0][3]=gen_code[0][0][1]='N';
  gen_code[0][0][0]= gen_code[0][0][2]='K';
  
  for (i=0; i<4; i++)
    gen_code[0][1][i]='T';
  
  gen_code[0][2][1]=gen_code[0][2][3]='S';
  gen_code[0][2][0]=gen_code[0][2][2]='R';
  
  gen_code[0][3][0]=gen_code[0][3][1]=gen_code[0][3][3]='I';
  
  gen_code[0][3][2]='M';
  
  gen_code[1][0][3]=gen_code[1][0][1]='H';
  gen_code[1][0][0]= gen_code[1][0][2]='Q';
  
  for (i=0; i<4; i++)
    gen_code[1][1][i]='P';

  for (i=0; i<4; i++)
    gen_code[1][2][i]='R';

  for (i=0; i<4; i++)
    gen_code[1][3][i]='L';

  gen_code[2][0][3]=gen_code[2][0][1]='D';
  gen_code[2][0][0]= gen_code[2][0][2]='E';

  for (i=0; i<4; i++)
    gen_code[2][1][i]='A';

  for (i=0; i<4; i++)
    gen_code[2][2][i]='G';

  for (i=0; i<4; i++)
    gen_code[2][3][i]='V';

  gen_code[3][0][3]=gen_code[3][0][1]='Y';
  gen_code[3][0][0]= gen_code[3][0][2]='*';
  

  for (i=0; i<4; i++)
    gen_code[3][1][i]='S';

  gen_code[3][2][3]=gen_code[3][2][1]='C';
  gen_code[3][2][0]='W'; 
  gen_code[3][2][2]='W';

  gen_code[3][3][3]=gen_code[3][3][1]='F';
  gen_code[3][3][0]= gen_code[3][3][2]='L';
  
 
  locate_stops();
} //End Mycoplasma_genetic_code::make_code





Genetic_code * create_genetic_code(Exchange *curr_exchange)
{
  Genetic_code *code;

  switch (curr_exchange->get_genetic_code()) {
  case UNIVERSAL:
    code = new Univ_genetic_code(curr_exchange);
    break;
  case VERT_MITO:
    code = new Vert_mito_genetic_code(curr_exchange);
    break;
  case YEAST_MITO:
    code = new Yeast_mito_genetic_code(curr_exchange);
    break;
  case MOLD_MITO:
    code = new Mold_mito_genetic_code(curr_exchange);
    break;
    case INVERT_MITO:
      code = new Invert_mito_genetic_code(curr_exchange);
    break;
  case CILIATE_NUC:
    code = new Ciliate_nuc_genetic_code(curr_exchange);
    break;
  case ECHINO_MITO:
    code = new Echino_mito_genetic_code(curr_exchange);
    break;
  case MYCOPLASMA:
	  code = new Mycoplasma_genetic_code(curr_exchange);
	  break;
  }

  return(code);
 
}
