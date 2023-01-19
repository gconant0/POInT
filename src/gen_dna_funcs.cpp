//Copyright 1999-2002 Gavin Conant


#include <iostream>
#include <math.h>
#include <string.h>
#include <fstream>
#include "gen_dna_funcs.h"


using namespace::std;
#define RECIP_SQRT_TWO 0.707106781

Molecule_Sequence::Molecule_Sequence ()
{
  cerr<<"Call to default constructor of Molecule_Sequnce class\n";
    sequence=0;
    size=0;
}


Molecule_Sequence::Molecule_Sequence (int sequence_length)
{
  size=sequence_length;
  sequence=new int [size];
  type=NUCLEIC;
}  //End Molecule_Sequence::Molecule_Sequence (int num_seqs, int sequence_length)


Molecule_Sequence::Molecule_Sequence (int sequence_length, DATATYPE tp)
{	
  size=sequence_length;
  sequence=new int [size];
  type=tp;
}

int Molecule_Sequence::operator[](int element)
{
  if (element<size)
    return(sequence[element]);
  else
    {
      cerr<<"Invalid Sequence location\n";
      return(0);
    }
}


void Molecule_Sequence::Assign_site(int site_num, int assignment)
{
   if (site_num<size)
     sequence[site_num]=assignment;
  else
      cerr<<"Invalid Sequence location\n";
}


int Molecule_Sequence::operator== (Molecule_Sequence &compareto)
{
  int i, retval=1;
 
  if ((size != compareto.Sequence_size()) || (type != compareto.get_datatype()))
    return(0);
  else 
    {
      for(i=0; i<size; i++)
	if (sequence[i] != compareto[i])
	    retval=0;
	
    }
  return(retval);

}



Molecule_Sequence & Molecule_Sequence::operator= (Molecule_Sequence & assign_from)
{
  int i;
  
 
  if (this == &assign_from)
      return(*this);

  else
  {
    if (size == assign_from.Sequence_size())
    {
      strcpy(name, assign_from.Sequence_name());
     
      for (i=0; i<size; i++)
         sequence[i]=assign_from[i];
    }
    else
       cerr<<"Incompatible Sequence length\n";
  }
  return(*this);
}


BOOL Molecule_Sequence::compare_elements(int element, Molecule_Sequence * other_seq, int other_element)
{
	if (sequence[element] == (*other_seq)[other_element])
		return(TRUE);
	else
		return(FALSE);
}


Sequence_dataset::Sequence_dataset()
{
  cerr<<"Call to default constructor of Sequence_dataset class\n";
    dummy=0;
    actual_sequences=0;
}



Sequence_dataset::Sequence_dataset(int num_taxa, int sequence_length)
{
  int i;
  num_sequences=num_taxa;

  dummy=new Molecule_Sequence(0);

  actual_sequences=new Molecule_Sequence * [num_taxa];
  for(i=0; i<num_sequences; i++)
    actual_sequences[i]=new Molecule_Sequence(sequence_length); 
}



Sequence_dataset::Sequence_dataset(int num_taxa, int *lenghts)
{
   int i;
  num_sequences=num_taxa;
  
  dummy=new Molecule_Sequence(0);
  
  actual_sequences=new Molecule_Sequence * [num_taxa];
  
  for(i=0; i<num_sequences; i++)
    actual_sequences[i]=new Molecule_Sequence(lenghts[i]); 
}




Sequence_dataset::Sequence_dataset(int num_taxa, int sequence_length, DATATYPE tp)
{
  int i;
  num_sequences=num_taxa;

  type=tp;

  if (type != NUCLEIC)
	dummy=new Molecule_Sequence(0, tp);
  else
	  dummy=new DNA_Ambig_Molecule_Sequence(0);

  actual_sequences=new Molecule_Sequence * [num_taxa];
  for(i=0; i<num_sequences; i++) {
	if(tp != NUCLEIC)
	  actual_sequences[i]=new Molecule_Sequence(sequence_length, tp); 
	else
		actual_sequences[i]=new DNA_Ambig_Molecule_Sequence(sequence_length);
  }
}



Sequence_dataset::Sequence_dataset(int num_taxa, int *lenghts, DATATYPE tp)
{
   int i;
  num_sequences=num_taxa;
  
  type=tp;

  if (type != NUCLEIC)
	dummy=new Molecule_Sequence(0, tp);
  else
	  dummy=new DNA_Ambig_Molecule_Sequence(0);
  
  actual_sequences=new Molecule_Sequence * [num_taxa];
  
  for(i=0; i<num_sequences; i++) {
	if (tp != NUCLEIC)
	  actual_sequences[i]=new Molecule_Sequence(lenghts[i], tp);
	else
		actual_sequences[i] = new DNA_Ambig_Molecule_Sequence(lenghts[i]);
  }
}



Sequence_dataset  & Sequence_dataset::operator=(Sequence_dataset  & assign_from)
{
  int i;
  
 if (this == &assign_from)
 return(*this);

  else
  {
    if (num_sequences == assign_from.Num_sequences())
    {
      for (i=0; i<num_sequences; i++)
      {   
		delete actual_sequences[i];
		if (assign_from.get_datatype() != NUCLEIC)
			actual_sequences[i]= new Molecule_Sequence(assign_from[i].Sequence_size(), assign_from.get_datatype());
		else
			actual_sequences[i] = new DNA_Ambig_Molecule_Sequence(assign_from[i].Sequence_size());
		*actual_sequences[i] = assign_from[i];
      }
    
    }
    else
       cerr<<"Incompatible Sequence numbers\n";
  }
 return(*this);
}


Molecule_Sequence& Sequence_dataset::operator[](int element)
{
 if (element<num_sequences)
    return(*actual_sequences[element]);
  else
    {
      cerr<<"Invalid Sequence number\n";
      return(*dummy);
    }  
}




Sequence_dataset::~Sequence_dataset()
{
  int i;
  for (i=0; i<num_sequences; i++)
    delete actual_sequences[i];

  delete[] actual_sequences;
  delete dummy;
}

DNA_Ambig_Molecule_Sequence::DNA_Ambig_Molecule_Sequence () : Molecule_Sequence()
{
	type=NUCLEIC;
}


DNA_Ambig_Molecule_Sequence::DNA_Ambig_Molecule_Sequence (int sequence_length) : 
Molecule_Sequence(sequence_length)
{
	type=NUCLEIC;
}


BOOL DNA_Ambig_Molecule_Sequence::compare_elements(int element, Molecule_Sequence * other_seq, int other_element)
{
	if (other_seq->get_datatype() != NUCLEIC)
		return(FALSE);
	else 
		return(bases_equal(sequence[element], (*other_seq)[other_element]));
}


Amino_acid_group::Amino_acid_group()
{
  cerr<<"Call to default constructor of class Amino_acid_group\n";
}


Amino_acid_group::Amino_acid_group(int ngroups)
{
  int i;

  num_groups=ngroups;
  for (i=0; i<20; i++)
	groups[i]=0;
 
}
 
Amino_acid_group::Amino_acid_group(int ngroups, const char *filename)
{
  int i, aa, group;
  char aa_symbol;
  ifstream groupsin;

  num_groups=ngroups;
  groupsin.open(filename);

  for (i=0; i<20; i++)
	groups[i]=0;
  
  for (i=0; i<20; i++)
    {
      groupsin>>aa_symbol>>group;
      aa=readchar_to_aa(aa_symbol);
      groups[aa]=group-1;
    }
}


Amino_acid_group & Amino_acid_group::operator=(Amino_acid_group  & assign_from)
{
  int i;
  if (num_groups==assign_from.get_num_groups())
  {
	  for (i=0; i<20; i++)
		groups[i]=assign_from.get_group(i);
	return((*this));
  }

  else
	  cerr<<"Assignment between group descriptions of unequal size\n";
  return(*this);
}
 
void Amino_acid_group::assign_to_group(int aa, int group)
{
  if (group<num_groups)
    groups[aa]=group;
  else
    cerr<<"Group assignment out of range\n";
}

Btree_node::Btree_node()
{
	children[0]=children[1]=0;
	parent=0;
	internal_node_num=-1;
	tip_num=-1;
}

Constrain_Param_Lookup::Constrain_Param_Lookup()
{
	cerr<<"Error: Invalid call to default constructor of Constrain_Param_Lookup\n";

}
	
Constrain_Param_Lookup::Constrain_Param_Lookup(int num_p)
{
	int i;

	num_categories =num_p;
	internal_node_count=0;	
	search_tree=new Btree_node[2*num_categories-1];

	//Holds the new independently variable parameters
	cat_proportions = new double [num_categories-1];

	//We assume we will never have more that 2^32 different parameters!
	pow2_masks[0]=1;
	for(i=1; i<32; i++)
		pow2_masks[i]=2*pow2_masks[i-1];
	
	internal_node_pos=num_categories;

	for(i=0; i<num_categories; i++) {
		search_tree[i].children[0]=0;
		search_tree[i].children[1]=0;
		search_tree[i].tip_num=i;
	}

	root=build_tree(0, num_categories-1);
}


double Constrain_Param_Lookup::get_param_value(int param_num)
{
	return(cat_proportions[param_num]);
}


double Constrain_Param_Lookup::get_prob_value(int category_num)
{
	return(traverse_tree(category_num, root));
}


void Constrain_Param_Lookup::set_param_proportion (int param_num, double param_val)
{
	//Note that values passed to this function should be between 0 and 1
	cat_proportions[param_num]=param_val;
}


Constrain_Param_Lookup::~Constrain_Param_Lookup()
{
	delete[] cat_proportions;
	delete[] search_tree;
}


Btree_node * Constrain_Param_Lookup::build_tree(int min_id, int max_id) 
{
	int diff, pow2_id, new_max1, new_min2;
	Btree_node *child1, *child2, *new_node;;

	diff=max_id-min_id+1;

	pow2_id=0;

	while(pow2_masks[pow2_id] < diff) pow2_id++;

	new_max1 = min_id + pow2_masks[pow2_id-1]-1;
	new_min2 = new_max1+1;

	if (new_max1 > min_id) 
		child1=build_tree(min_id, new_max1);
	else
		child1=&search_tree[min_id];

	if(new_min2 < max_id)
		child2=build_tree(new_min2, max_id);
	else
		child2=&search_tree[max_id];

	new_node=get_next_free_node();

	child1->parent=child2->parent=new_node;
	new_node->children[0]=child1;
	new_node->children[1]=child2;
	new_node->left_max=new_max1;

	cat_proportions[new_node->internal_node_num]=(double)(new_max1-min_id+1)/(double)diff;
	return(new_node);

}

Btree_node * Constrain_Param_Lookup::get_next_free_node()
{
	while(search_tree[internal_node_pos].children[0] != 0)	internal_node_pos++;
	search_tree[internal_node_pos].internal_node_num=internal_node_count++;
	return(&search_tree[internal_node_pos]);
}


double Constrain_Param_Lookup::traverse_tree(int index, Btree_node *curr_node)
{
	if (curr_node->tip_num==index)
		return(1.0);
	else {
		if (index > curr_node->left_max)
			return((1.0-cat_proportions[curr_node->internal_node_num])*traverse_tree(index, curr_node->children[1]));
		else
			return(cat_proportions[curr_node->internal_node_num]*traverse_tree(index, curr_node->children[0]));
	}
}
	

BinaryString::BinaryString ()
{
	cerr<<"Error: call to BinaryString default constructor: no size given\n";

}



BinaryString::BinaryString (int l)
{
	int i;

	len=l;
	int_val_correct=FALSE;
	the_string=new int [len];
	for(i=0; i<len; i++)
		the_string[i]=0;
}


BinaryString& BinaryString::operator= (BinaryString &assign_from)
{
	int i;

	int_val_correct=FALSE;
	if (len == assign_from.get_len()) {
		for(i=0; i<len; i++)
			the_string[i]=assign_from.element_n(i);
	}
	else
		cerr<<"Error: Assignment to BinaryString from unequal size string\n";
		
	
	return(*this);
}

int BinaryString::operator[] (int n)
{
	if ((n>=0) && (n<len)){
		return(element_n(n));
	}
	else {
		return(-1);
	}
}

int BinaryString::element_n(int n)
{
	return(the_string[n]);
}
	

void BinaryString::set_element_n(int n, int val)
{
	if (((val ==0) || (val==1)) && ((n>=0) && (n<len)))
		the_string[n]=val;
	int_val_correct=FALSE;
}


BinaryString::~BinaryString()
{
	delete[] the_string;
}

int BinaryString::string_int_value()
{
	int i;

	if (int_val_correct != TRUE) {
		intval=0;
		if (len <= 32) {
			for(i=0; i<len; i++) 
				intval = intval & (the_string[i]<<i);
			int_val_correct=TRUE;
		
		}
		else {
			cerr<<"Error: Binary string is too long to represent as 32-bit integer\n";
			intval=-1;
		}
	}
	return(intval);

}


extern int readchar_to_base(char inbase)
  //Converts standard nucleic acid character representations
  //to my internal integer representation.  This allows direct lookups
  //in arrays like transition probability matrices
{
  switch (inbase)
    {
    case 'a':
    case 'A':
		return(0);
    case 'c':
    case 'C':
		return(1);
    case 'g':
    case 'G':
		return(2);
    case 't':
    case 'T':
		return(3);
    case 'y':
    case 'Y':
		return(4);
    case 'r':
    case 'R':
		return(5);
    case 'n':
    case 'N':
		return(6);      
    case '-':
		return(7);
	case 's':
	case 'S':
		return(8);
	case 'w':
	case 'W':
		return(9);
	case 'm':
	case 'M':
		return(10);
	case 'k':
	case 'K':
		return(11);
	case 'v':
	case 'V':
		return(12);
	case 'h':
	case 'H':
		return(13);
	case 'd':
	case 'D':
		return(14);
	case 'b':
	case 'B':
		return(15);
    default:
      return(16);
    }

} //End readchar_to_base


int readchar_to_aa(char inaa) 
  //Converts standard amino acid character representations
  //to my internal integer representation.  This allows direct lookups
  //in arrays like transition probability matrices
{
  switch(inaa)
    {
    case 'a':
    case 'A':
      return(0);
    case 'c':
    case 'C':
      return(1);
    case 'd':
    case 'D':
      return(2);
    case 'e':
    case 'E':
      return(3);
    case 'f':
    case 'F':
      return(4);
    case 'g':
    case 'G':
      return(5);
    case 'h':
    case 'H':
      return(6);
    case 'i':
    case 'I':
      return(7);
    case 'k':
    case 'K':
      return(8);
    case 'l':
    case 'L':
      return(9);
    case 'm':
    case 'M':
      return(10);
    case 'n':
    case 'N':
      return(11);
    case 'p':
    case 'P':
      return(12);
    case 'q':
    case 'Q':
      return(13);
    case 'r':
    case 'R':
      return(14);
    case 's':
    case 'S':
      return(15);
    case 't':
    case 'T':
      return(16);
    case 'v':
    case 'V':
      return(17);
    case 'w':
    case 'W':
      return(18);
    case 'y':
    case 'Y':
      return(19);
    case 'x':
    case 'X':
      return(20);
    case '-':
      return(21);
    case '$':
	case '*':
      return(22);
      break;
    case '#':
      return(23);
      break;
    default:
      return(24);
    }
  //Return for xlC compile--should never be reached
  return(-1);
}

 
extern int loss_state_to_dupl(DUPL_LOSS_STATES loss_state)
{
	//Translates dupl model named states to integers for reference in data and trans. prob
	//matricies

	switch (loss_state) {
	case BOTH_PRESENT:
		return(3);
		break;
	case COPY1:
		return(1);
		break;
	case COPY2:
		return(2);
		break;
	case BOTH_PRESENT_FIXED:
		return(6);
		break;
	case BOTH_1_BIAS:
		return(7);
		break;
	case BOTH_2_BIAS:
		return(8);
		break;
	case COPY1_OR_BOTH:
		return(9);
		break;
	case COPY2_OR_BOTH:
		return(10);
		break;
	case BOTH_FIXED_SUBF:
		return(11);
		break;
	case COPY1_BIAS:
		return(12);
		break;
	case COPY2_BIAS:
		return(13);
		break;
	case GENERIC_SINGLE_COPY:
		return(4);
		break;
	case MISSING:
		return(5);
		break;
	case LOST:
		return(0);
		break;
	}

}


extern DUPL_LOSS_STATES dupl_to_loss_state (int dupl)
{
	//Translates integers to dupl model named states for reference in data and trans. prob
	//matricies

	switch (dupl) {
        case 3:
            return(BOTH_PRESENT);
            break;
        case 1:
            return(COPY1);
            break;
        case 2:
            return(COPY2);
            break;
        case 6:
            return(BOTH_PRESENT_FIXED);
            break;
        case 7:
            return(BOTH_1_BIAS);
            break;
        case 8:
            return(BOTH_2_BIAS);
            break;
        case 9:
            return(COPY1_OR_BOTH);
            break;
        case 10:
            return(COPY2_OR_BOTH);
            break;
        case 11:
            return(BOTH_FIXED_SUBF);
            break;
        case 12:
            return(COPY1_BIAS);
            break;
        case 13:
            return (COPY2_BIAS);
            break;
        case 4:
            return(GENERIC_SINGLE_COPY);
            break;
        case 5:
            return(MISSING);
            break;
        case 0:
            return(LOST);
            break;
        default:
            return(BOTH_PRESENT);
            break;
	}


}

extern char num_to_dupl_data(int indupl)
{
	switch(dupl_to_loss_state(indupl)) {
	case LOST:
		return('N');
		break;
	case COPY1:
		return('1');
		break;
	case COPY2:
		return('2');
		break;
	case BOTH_PRESENT:
		return('B');
		break;
	case GENERIC_SINGLE_COPY:
		return('#');
		break;
	case MISSING:
		return('-');
		break;
	case BOTH_PRESENT_FIXED:
		return('F');
		break;
	case BOTH_1_BIAS:
		return('X');
		break;
	case BOTH_2_BIAS:
		return('Y');
		break;
	case BOTH_FIXED_SUBF:
		return('S');
		break;
	case COPY1_BIAS:
		return('I');
		break;
	case COPY2_BIAS:
		return('J');
		break;
	case COPY1_OR_BOTH:
		return('C');
		break;
	case COPY2_OR_BOTH:
		return('D');
		break;
	}
}


extern int readchar_to_dupl(char indupl)
//The intermediate loss_state_to_dupl function is silly, but it allows us to use names in
//other places and be sure we are referring to the same state
{
	switch(indupl) {
	case 'n':
	case 'N':
		return(loss_state_to_dupl(LOST));
		break;
	case '1':
		return(loss_state_to_dupl(COPY1));
		break;
	case '2':
		return(loss_state_to_dupl(COPY2));
		break;
	case 'B':
	case 'b':
		return(loss_state_to_dupl(BOTH_PRESENT));
		break;
	case '#':
		return(loss_state_to_dupl(GENERIC_SINGLE_COPY));
		break;
	case '-':
		return(loss_state_to_dupl(MISSING));
		break;
	case 'c':
	case 'C':
		return(loss_state_to_dupl(COPY1_OR_BOTH));
		break;
	case 'd':
	case 'D':
		return(loss_state_to_dupl(COPY2_OR_BOTH));
		break;
	case 'F':
	case 'f':
		cerr<<"ERROR: Hidden state code used in input: F.  Continuing...\n";
		return(loss_state_to_dupl(BOTH_PRESENT_FIXED));
		break;
	case 'x':
	case 'X':
		cerr<<"ERROR: Hidden state code used in input: X.  Continuing...\n";
		return(loss_state_to_dupl(BOTH_1_BIAS));
		break;
	case 'y':
	case 'Y':
		cerr<<"ERROR: Hidden state code used in input: Y.  Continuing...\n";
		return(loss_state_to_dupl(BOTH_2_BIAS));
		break;
	case 's':
	case 'S':
		cerr<<"ERROR: Hidden state code used in input: S.  Continuing...\n";
		return(loss_state_to_dupl(BOTH_FIXED_SUBF));
		break;
	case 'i':
	case 'I':
		cerr<<"ERROR: Hidden state code used in input: I.  Continuing...\n";
		return(loss_state_to_dupl(COPY1_BIAS));
		break;
	case 'j':
	case 'J':
		cerr<<"ERROR: Hidden state code used in input: J.  Continuing...\n";
		return(loss_state_to_dupl(COPY2_BIAS));
		break;
	default:
		return(loss_state_to_dupl(MISSING));
		break;
	}
}


DUPL_LOSS_STATES get_observable_dupl_state(DUPL_LOSS_STATES inval)
//Takes the internal duplication states and maps them to the observable states
{
	if (observable_dupl_state(inval) == TRUE)
		return(inval);
	else {
		switch (inval) {
			case BOTH_PRESENT_FIXED:
			case BOTH_1_BIAS:
			case BOTH_2_BIAS:
			case BOTH_FIXED_SUBF:
				return(BOTH_PRESENT);
				break;
			case COPY1_BIAS:
				return(COPY1);
				break;
			case COPY2_BIAS:
				return(COPY2);
				break;
            default:
                return(BOTH_PRESENT);
                break;
		}
	
	}
}


extern char num_to_base(int inbase)
  //Converse of readchar_to_base--converts from my integer representation
  //to the standard codes
{
  switch (inbase)
    {
    case 0:
      return('A');
    case 1:
      return('C');
    case 2:
      return('G');
    case 3:
      return('T');
    case 4:
      return('Y');
    case 5:
      return('R');
    case 6:
      return('N');
    case 7:
      return('-');
	case 8:
		return('S');
	case 9:
		return('W');
	case 10:
		return('M');
	case 11:
		return('K');
	case 12:
		return('V');
	case 13:
		return('H');
	case 14:
		return('D');
	case 15:
		return('B');
    default:
      return(' ');
    }

} //End num_to_base



char num_to_aa(int inaa)
  //Converse of readchar_to_aa--converts from my integer representation
  //to the standard codes
{
 switch(inaa)
    {
    case 0:
      return('A');
    case 1:
      return('C');
    case 2:
      return('D');
    case 3:
      return('E');
    case 4:
      return('F');
    case 5:
      return('G');
    case 6:
      return('H');
    case 7:
      return('I');
    case 8:
      return('K');
    case 9:
      return('L');
    case 10:
      return('M');
    case 11:
      return('N');
    case 12:
      return('P');
    case 13:
      return('Q');
    case 14:
      return('R');
    case 15:
      return('S');
    case 16:
      return('T');
    case 17:
      return('V');
    case 18:
      return('W');
    case 19:
      return('Y');
    case 20:
      return('X');
    case 21:
      return('-');
    case 22:
      return('*');
    case 23:
      return('#');
    default:
      return(' ');
    }

}


extern BOOL is_base(char inbase)
{
  switch (inbase)
    {
    case 'a':
    case 'A':
      return(TRUE);
    case 'c':
    case 'C':
      return(TRUE);
    case 'g':
    case 'G':
      return(TRUE);
    case 't':
    case 'T':
      return(TRUE);
    case 'y':
    case 'Y':
      return(TRUE);
    case 'r':
    case 'R':
      return(TRUE);
    case 'n':
    case 'N':
      return(TRUE);
	case 'M':
	case 'm':
		return(TRUE);
	case 's':
	case 'S':
		return(TRUE);
	case 'K':
	case 'k':
		return(TRUE);
	case 'w':
	case 'W':
		return(TRUE);
	case 'v':
	case 'V':
		return(TRUE);
	case 'h':
	case 'H':
		return(TRUE);
	case 'd':
	case 'D':
		return(TRUE);
	case 'b':
	case 'B':
		return(TRUE);
    case '-':
      return(TRUE);
    default:
      return(FALSE);
    }
 
}  //End is_base



extern BOOL is_aa(char inaa)
{
  if ((inaa>=65 && inaa<=90) || (inaa>=97 && inaa<=122) || inaa=='-')
    {
      switch (inaa)
	{
	case 'b':
	case 'B':
	case 'j':
	case 'J':
	case 'o':
	case 'O':
	case 'u':
	case 'U':
	case 'z':
	case 'Z':
	  return (FALSE);
	  break;
	default:
	  return(TRUE);
	  break;
	}
     
    }
 
  return(FALSE);
}



extern BOOL is_dupl_data(char indupl)
{
	switch (indupl) {
	case 'N':
	case 'n':
	case '1':
	case '2':
	case 'b':
	case 'B':
	case 'c':
	case 'C':
	case 'd':
	case 'D':
	case '-':
	case '#':
	case 'X':
	case 'x':
	case 'Y':
	case 'y':
	case 'F':
	case 'f':
	case 's':
	case 'S':
	case 'i':
	case 'I':
		return(TRUE);
		break;
	default:
		return(FALSE);

	}
}


extern BOOL base_is_ambig(int base)
{
	switch(num_to_base(base))
	{
	case 'Y':
	case 'R':
	case 'N':
	case 'S':
	case 'W':
	case 'M':
	case 'K':
	case 'V':
	case 'H':
	case 'D':
	case 'B':
		return(TRUE);
	default:
		return(FALSE);
	}

}

extern BOOL valid_ambig_spec(int ambig_base, int test_base)
{
  if (base_is_ambig(ambig_base)==TRUE) {
    switch (num_to_base(ambig_base)) {
    case 'Y':
      if((num_to_base(test_base) == 'T') || (num_to_base(test_base) == 'C'))
		return(TRUE);
      else
		return(FALSE);
      break;
    case 'R':
      if((num_to_base(test_base) == 'T') || (num_to_base(test_base) == 'C'))
		return(FALSE);
      else
		return(TRUE);
      break;
	case 'S':
		if((num_to_base(test_base) == 'G') || (num_to_base(test_base) == 'C'))
			return(TRUE);
		else
			return(FALSE);
		break;
	case 'W':
		if((num_to_base(test_base) == 'T') || (num_to_base(test_base) == 'A'))
			return(TRUE);
		else
			return(FALSE);
		break;
	case 'M':
		 if((num_to_base(test_base) == 'A') || (num_to_base(test_base) == 'C'))
			return(TRUE);
		else
			return(FALSE);
		break;
	case 'K':
		 if((num_to_base(test_base) == 'T') || (num_to_base(test_base) == 'G'))
			return(TRUE);
		else
			return(FALSE);
		break;
	case 'V':
		if((num_to_base(test_base) == 'A') || (num_to_base(test_base) == 'C') || (num_to_base(test_base) == 'G'))
			return(TRUE);
		else
			return(FALSE);
		break;
	case 'H':
		if((num_to_base(test_base) == 'A') || (num_to_base(test_base) == 'C') || (num_to_base(test_base) == 'T'))
			return(TRUE);
		else
			return(FALSE);
		break;
	case 'D':
		if((num_to_base(test_base) == 'A') || (num_to_base(test_base) == 'T') || (num_to_base(test_base) == 'G'))
			return(TRUE);
		else
			return(FALSE);
		break;
	case 'B':
		if((num_to_base(test_base) == 'T') || (num_to_base(test_base) == 'C') || (num_to_base(test_base) == 'G'))
			return(TRUE);
		else
			return(FALSE);
		break;
    case 'N':
      return(TRUE);
      break;
    default:
        return(FALSE);
        break;
    }
  }
  else
    {
      if (ambig_base==test_base)
	return(TRUE);
      else
	return(FALSE);
    }
}

extern int get_num_ambig_states(int inbase)
{
	switch (num_to_base(inbase) ){
	case 'N':
		return(4);
	case 'Y':
	case 'R':
	case 'S':
	case 'W':
	case 'K':
	case 'M':
		return(2);
	case 'V':
	case 'H':
	case 'D':
	case 'B':
		return(3);
	case 'A':
	case 'C':
	case 'G':
	case 'T':
	case '-':
	default:
		return(1);

	}



}


extern int get_ambig_state_n(int inbase, int ambig_num)
{
	switch (num_to_base(inbase)) {
	case 'Y':
		if (ambig_num == 0)
			return(readchar_to_base('C'));
		else
			return(readchar_to_base('T'));
	case 'R':
		if (ambig_num == 0)
			return(readchar_to_base('A'));
		else
			return(readchar_to_base('G'));
	case 'S':
		if (ambig_num == 0)
			return(readchar_to_base('C'));
		else
			return(readchar_to_base('G'));
	case 'W':
		if (ambig_num == 0)
			return(readchar_to_base('A'));
		else
			return(readchar_to_base('T'));
	case 'K':
		if (ambig_num == 0)
			return(readchar_to_base('T'));
		else
			return(readchar_to_base('G'));
	case 'M':
		if (ambig_num == 0)
			return(readchar_to_base('A'));
		else
			return(readchar_to_base('C'));
	case 'V':
		if (ambig_num == 0)
			return(readchar_to_base('A'));
		else if (ambig_num == 1)
			return(readchar_to_base('G'));
		else
			return(readchar_to_base('C'));
	case 'H':
		if (ambig_num == 0)
			return(readchar_to_base('A'));
		else if (ambig_num == 1)
			return(readchar_to_base('T'));
		else
			return(readchar_to_base('C'));
	case 'D':
		if (ambig_num == 0)
			return(readchar_to_base('A'));
		else if (ambig_num == 1)
			return(readchar_to_base('G'));
		else
			return(readchar_to_base('T'));
	case 'B':
		if (ambig_num == 0)
			return(readchar_to_base('T'));
		else if (ambig_num == 1)
			return(readchar_to_base('G'));
		else
			return(readchar_to_base('C'));
	case 'N':
		if (ambig_num == 0)
			return(readchar_to_base('A'));
		else if (ambig_num == 1)
			return(readchar_to_base('C'));
		else if (ambig_num == 2)
			return(readchar_to_base('G'));
		else
			return(readchar_to_base('T'));
	case 'A':
	case 'C':
	case 'G':
	case 'T':
	case '-':
	default:
		return(inbase);
	}

}


extern BOOL bases_equal(int inbase1, int inbase2) 
{
	int i, j;
	BOOL retval=FALSE;

	for(i=0; i<get_num_ambig_states(inbase1); i++)
		for(j=0; j<get_num_ambig_states(inbase2); j++)
			if (get_ambig_state_n(inbase1, i) ==
				get_ambig_state_n(inbase2, j))
				retval=TRUE;

	return(retval);
}


extern BOOL observable_dupl_state(DUPL_LOSS_STATES the_state)
{
	switch (the_state) {
	case COPY1:
	case COPY2:
	case GENERIC_SINGLE_COPY:
	case MISSING:
	case LOST:
	case BOTH_PRESENT:
	case COPY1_OR_BOTH:
	case COPY2_OR_BOTH:
		return(TRUE);
		break;
	case BOTH_PRESENT_FIXED:
	case BOTH_1_BIAS:
	case BOTH_2_BIAS:
	case BOTH_FIXED_SUBF:
	case COPY1_BIAS:
	case COPY2_BIAS:
		return(FALSE);
		break;
	}
}


extern SNP_STATE snp_to_snpstate(int state)
{
	switch(state) {
	case 0:
		return(TYPE_A);
		break;
	case 1:
		return(TYPE_B);
		break;
	case 2:
		return(HETEROZYGOUS);
		break;
	case 3:
		return(SNP_ABSENT);
		break;
	default:
		return(DATA_MISSING);
		break;
	
	};
}

extern int snpstate_to_snp(SNP_STATE inval)
{
	switch(inval) {
	case TYPE_A:
		return(0);
		break;
	case TYPE_B:
		return(1);
		break;
	case HETEROZYGOUS:
		return(2);
		break;
	case SNP_ABSENT:
		return(3);
		break;
	case DATA_MISSING:
		return(4);
		break;
	};
}

extern char num_to_snp_state(int instate)
{
	switch(snp_to_snpstate(instate)) {
	case TYPE_A:
		return('0');
		break;
	case TYPE_B:
		return('1');
		break;
	case HETEROZYGOUS:
		return('B');
		break;
	case SNP_ABSENT:
		return('-');
		break;
	case DATA_MISSING:
		return('X');
		break;
	};
}


extern int readsnp_to_snpstate(char instate)
{
	switch(instate) {
        case '0':
            return(snpstate_to_snp(TYPE_A));
            break;
        case '1':
            return(snpstate_to_snp(TYPE_B));
            break;
        case 'B':
        case 'b':
            return(snpstate_to_snp(HETEROZYGOUS));
            break;
        case '-':
            return(snpstate_to_snp(SNP_ABSENT));
            break;
        case 'x':
        case 'X':
            return(snpstate_to_snp(DATA_MISSING));
            break;
        default:
            return(snpstate_to_snp(TYPE_A));
            break;
	}
}


extern BOOL is_snpstate(char instate)
{
	switch (instate) {
	case '0':
	case '1':
	case 'b':
	case 'B':
	case '-':
	case 'x':
	case 'X':
		return(TRUE);
		break;
	default:
		return(FALSE);
		break;
	}
}

extern double observed_basefreqs (Sequence_dataset *curr_data, int freq)
  //Calculates the observed basefrequencies (excluding gaps) of a sequence or
  //sequence alignment
{
  int i, j, count=0, total_count=0;

 
  if (freq>=0 && freq<=3)
    {
      for (i=0; i<curr_data->Num_sequences(); i++)
	
	for (j=0; j<(*curr_data)[i].Sequence_size(); j++)
	  {
	    if (((*curr_data)[i][j]>=0) &&((*curr_data)[i][j]<4)) {
	      if ((*curr_data)[i][j]==freq)
		count++;
	      total_count++;
	    }
	  }
      
      return(((double)count)/total_count);
    }

  else
    {
      cerr<<"Invalid frequency number in observed_basefreqs\n";
      return(0);
    } 
}  //End observed_basefreqs



extern double observed_codon_basefreqs(Sequence_dataset *curr_data, int codon_pos, int freq)
  //Calculates the base frequencies at each codon position for a sequence or sequence
  //alignment.  Requires that the sequence length be a multiple of three
{
  int i, j, count=0, total_count=0;
  BOOL fail=FALSE;

  if (!(freq>=0 && freq<=3))
    fail=TRUE;
  if(codon_pos<0 || codon_pos>2)
    fail=TRUE;

  if (fail==FALSE)
    for (i=0; i<curr_data->Num_sequences(); i++) {
      if ((*curr_data)[i].Sequence_size()%3==0)
	for (j=0; j<(*curr_data)[i].Sequence_size()/3; j++)
	  { 
	    if (((*curr_data)[i][3*j+codon_pos]>=0) &&((*curr_data)[i][3*j+codon_pos]<4)) {
	      if ((*curr_data)[i][3*j+codon_pos]==freq)
		count++;
	      total_count++;
	    }
	  }
      else
	fail=TRUE;
  
    }

  if (fail==FALSE)
    return(((double)count)/total_count);
  else
    {
      cerr<<"Invalid data for observed_codon_basefreqs\n";
      return(0);
    } 

}


extern Sequence_dataset * remove_gaps (Sequence_dataset *orig_data, BOOL protein_seq)
  //Given a sequence alignment, returns an alignment where all gap positions have been removed.
  //Warning--allocation memory to store the new sequence alignment
{
  int i, j, ntaxa, new_nchars, nchars;
  char (*convert_to_char)(int);
  BOOL found_gap, first=TRUE;
  Non_gap_sites *non_gap_list, *start;
  Sequence_dataset *new_dataset;

  ntaxa=orig_data->Num_sequences();
  nchars=(*orig_data)[0].Sequence_size();

 
  if (protein_seq == TRUE)
    convert_to_char=&num_to_aa;
  else
    convert_to_char=&num_to_base;

  new_nchars=nchars;
    
  non_gap_list=new Non_gap_sites;
  start=non_gap_list;
  start->next=0;

  for(j=0; j<nchars; j++)
    {
      found_gap=FALSE;
      i=0;
      while ((found_gap==FALSE) && (i<ntaxa))
	{
	  if (convert_to_char((*orig_data)[i][j])=='-')
	    found_gap=TRUE;
	  i++;
	}
      if (found_gap==TRUE)
	new_nchars--;
      
      else
	{
	  initialize_next_list_element (non_gap_list, first);
	  non_gap_list->site_num=j;
	  first=FALSE;
	}
    }
  
  new_dataset=new Sequence_dataset(ntaxa, new_nchars);
  
  for (i=0; i<ntaxa; i++)
      (*new_dataset)[i].Assign_name((*orig_data)[i].Sequence_name());
	
    non_gap_list=start;

    for(j=0; j<new_nchars; j++)
      {
	for (i=0; i<ntaxa; i++)
	  (*new_dataset)[i].Assign_site(j, (*orig_data)[i][non_gap_list->site_num]);
	non_gap_list=non_gap_list->next;
      }      
    
    non_gap_list=start;
    while(start !=0) {
      start=non_gap_list->next;
      delete non_gap_list;
      non_gap_list=start;
    }

    return(new_dataset);
}


char to_upper(char inchar)
{
  if (inchar>=97 && inchar<=122)
    return(inchar-32);
  else
    return(inchar);
}  //End to_upper




char to_lower(char inchar)
{
 if (inchar>=65 && inchar<=90)
    return(inchar+32);
  else
    return(inchar);

}  //End to_lower




int word_match(std::string my_string, std::string my_word)
{
    std::transform(my_string.begin(), my_string.end(), my_string.begin(), ::tolower);
    std::transform(my_word.begin(), my_word.end(), my_word.begin(), ::tolower);
     if (my_string.find(my_word) == std::string::npos)
        return(0);
    else
        return(1);
   
}

int word_match(char string1[], int stringlen, char word[], int wordlen)
{
    int i=0, returnval=0;
    char *parstring;
    
    while (i<stringlen && returnval==0)
    {
        parstring=&string1[i];
        returnval=exact_match(parstring, stringlen-i, word, wordlen);
        i++;
    }
    
    return(returnval);
}





int loc_word_match (char string1[], int stringlen, char word[], int wordlen)
{
 int i=0, returnval=-1;
  char *parstring;

  while (i<stringlen && returnval==-1)
    {
      parstring=&string1[i];

   if(exact_match(parstring, stringlen-i, word, wordlen)==1)
     returnval=i;
      i++;
    }

  return(returnval);


}


int exact_match (char string1[], int stringlen, char word[], int wordlen)
{
  char  *parstring, *parword;
  if (to_ucase(string1[0])==to_ucase(word[0]))
    {
      if (wordlen>1)
	{
	      parword=word+1;
	      parstring=string1+1;
	      return(exact_match(parstring, stringlen-1, parword, wordlen-1));
	}
      else
	return(1);
    }
  else
    return(0);
}


int exact_match (Molecule_Sequence *seq, int pos1,  Molecule_Sequence *motif, int pos2)
{
	//New function (3/10/06 which allows ambig specs in DNA sequences if 
	//class DNA_Ambig_Molecule_Sequence is used
  if (seq->compare_elements(pos1, motif, pos2) ==TRUE )
    {
      if (motif->Sequence_size()-pos2>1)
	{
	  return(exact_match(seq, pos1+1, motif, pos2+1));
	}
      else
	return(1);
    }
  else
    return(0);
}


int string_to_int(const char *instring)
{
  int i=0, value;
  BOOL neg=FALSE;

  while(instring[i]<48 || instring[i]>57) {
    if ((instring[i]) == '-' &&( instring[i+1]>=48 && instring[i+1]<=57))
      neg=TRUE;
    i++;
  }  
  value=instring[i]-48;
  i++;

  while (instring[i]>=48 && instring[i]<=57)
       {
         value=value*10+(instring[i]-48);
         i++;
       }

  if (neg==TRUE)
    value *=-1;

  return(value);
} //End string_to_int




double string_to_float(const char *instring)
{
  int i=0; 
  BOOL neg=FALSE;
  double value=0, tenths=0.1;

  while(instring[i]<48 || instring[i]>57)
    i++;
  
  if (instring[i-1]=='-')
    neg=TRUE;
  value=instring[i]-48;
  i++;

  while (instring[i]>=48 && instring[i]<=57)
       {
         value=value*10+(instring[i]-48);
         i++;
       }
  if (instring[i]=='.')
    {
      i++;
      while (instring[i]>=48 && instring[i]<=57)
       {
         value=value+((instring[i]-48)*tenths);
	 tenths*=0.1;
         i++;
       }
    }

  if (neg==TRUE)
    value=-1.0*value;

 return(value);

}


extern void double_to_string (char *the_string, int string_len, int precision, double val)
{
    int i, pre_dec, exp_size=0, neg=0, less_than_one, pos=0, exp_val,exp_digits;
    double temp=1;
    BOOL long_enough=TRUE, is_zero=FALSE, exp_form=FALSE;

    temp=fabs(val);

    if (val<0) {
        neg=1;
    }
    else {
        if (val==0.0)
            is_zero=TRUE;
    }
 
    if (log10(temp)>=0)
        less_than_one=0;
    else
        less_than_one=1;

    if (is_zero==FALSE) {
        exp_size=less_than_one + trunc(log10(temp));
        exp_val = std::abs((int)trunc(log10(temp)));
        exp_digits =(int)log10(exp_val)+1;
    }
     
 
    if (string_len<(2+exp_size+precision+neg+less_than_one))
        long_enough=FALSE;
 

    if (long_enough==TRUE) {
          if (is_zero==TRUE)
              strcpy(the_string, "0.0");
          else {
              if (((int)log10(temp)+1+2*less_than_one+(neg+1+precision))>string_len)
                  exp_form=TRUE;
          
              if (temp >1e9) exp_form=TRUE;
              
              if (exp_form==FALSE) {
                  if (neg==1) {
                      the_string[pos]='-';
                      pos++;
                  }

                  if (less_than_one==1) {
                      the_string[pos]='0';
                      pos++;
                  }
                  else {
                      pre_dec=((int)log10(temp))+1;
                
                      for (i=pre_dec; i>0; i--) {
                          the_string[pos]=( (int) (temp/((int)pow(10,i-1))) )+48;
                          temp-=( (int) (temp/((int)pow(10,i-1))) )*((int)pow(10,i-1));
                          pos++;
                      }
                  }
                  the_string[pos]='.';
                  pos++;
                  
                  temp=10.0*(temp-(int)temp);
                  
                  for (i=-1; i>=-precision; i--) {
                      the_string[pos]=((int)temp)+48;
                      temp=10.0*(temp-((int)temp));
                      pos++;
                  }
                  
                  the_string[pos]='\0';
              }
              else {
                 // cout<<"Storing "<<val<<" in expontial form. Prec. is "<<precision<<"\n";
                  if (neg==1) {
                      the_string[pos]='-';
                      pos++;
                  }

                  pre_dec=(int)log10(temp)+1;
                  the_string[pos]=(int)(temp/pow(10,pre_dec-1))+48;
                  pos++;
                  the_string[pos]='.';
                  pos++;
                  temp-=((double)((int)(temp/pow(10,pre_dec-1))))*(pow(10,pre_dec-1));
                  
                 // cout<<"Entering loop: temp is "<<temp<<endl;
                  
                  for (i=0; i<precision; i++) {
                      the_string[pos]=(int)(temp/pow(10,pre_dec-(i+2)))+48;
                      temp-=((double)((int)(temp/pow(10,pre_dec-(i+2)))))*pow(10,pre_dec-(i+2));
                     // cout<<i<<": "<<temp<<endl;
                      pos++;
                  }
                  
                  the_string[pos]='e';
                  pos++;
                  if (less_than_one==1)
                      the_string[pos]='-';
                  else
                    the_string[pos]='+';
                  pos++;
                  
                  
                  temp=exp_val;
                  //cout<<"Entering loop: temp is "<<temp<<endl;
                  for(i=0; i<exp_digits; i++) {
                      the_string[pos]=(int)(temp/pow(10,exp_digits-i-1))+48;
                      temp-=((double)((int)(temp/(pow(10,exp_digits-i-1)))))*pow(10,exp_digits-i-1);
                      //cout<<i<<": "<<temp<<endl;
                       pos++;
                  }
                  the_string[pos]='\0';
                  //cout<<"Result was "<<the_string<<" with "<<exp_val<<" and "<<pre_dec<<" and size "<<exp_size<<" digits: "<<exp_digits<<endl;
              }
          }
    }
    else
        cerr<<"String too short to hold requested value\n";
}



void int_to_string (char *the_string, int string_len, int val)
{
  int len=0, digits, i=0;
  
  
  if (val!=0)
    {
      digits=(int)(log10(fabs(val)))+1;
     
      if (val<0)
	len=digits+1;
      else 
	len=digits;
      
       
      if (len>string_len-1)
	cerr<<"String too short to hold requested value\n";
      else
	{
	  while (digits>0)
	    {
	      digits--;
	      
	      the_string[i]=(int)(val/(pow(10,digits)))+48;
	      i++;
	      val-=((int)(val/(pow(10,digits))))*(pow(10,digits));
	    }
	  the_string[i]='\0';
	}
    }
  else
    { 
    if (string_len>2)
      {
	the_string[0]=48;
	the_string[1]='\0';
      }
    else
      cerr<<"String too short to hold requested value\n"; 
    }
}


void to_ucase(char *instring)
{
	int i;

	for(i=0; i<strlen(instring); i++)
	{
		if(instring[i]>=97 && instring[i]<=122)
			instring[i]=instring[i]-32;
	}
	

}


char to_ucase(char inchar)
{
  if(inchar>=97 && inchar<=122)
    return(inchar-32);
  else
    return(inchar);
}


BOOL is_transition (int base_1, int base_2)
{
  switch (base_1)
    {
    case(0):
      if (base_2==2)
	return(TRUE);
      else
	return(FALSE);
      break;
    case(1):
      if (base_2==3)
	return(TRUE);
      else
	return(FALSE);
      break;
    case(2):
      if (base_2==0)
	return(TRUE);
      else
	return(FALSE);
      break;
    case(3):
      if (base_2==1)
	return(TRUE);
      else
	return(FALSE);
      break;
    default:
      return(FALSE);
      break;
    }

}



double integer_power(double probablity, int affects)
{
  if (affects!=0)
    return(probablity);
  else
    return(1);
}



double logadd(double val1, double val2)
  //Takes two values that are the log of 2 numbers and returns
  //the log of the sum of those two numbers (this function avoids
  //The precision problems that can result in adding two very small
  //numbers)
{
	
	if(val1 == MIN_FLOAT ) return(val2);
	else if (val2== MIN_FLOAT) return(val1);
	else {
	  if ((val1>val2) && (val1-710>val2))
		  return(val1);
	  else if ((val2>val1) && (val2-710>val1))
		  return(val2);
	  else
	    return((double)(val2+log(1+exp((long double)(val1-val2)))));
  }
}  //End logadd



double logsubtract(double val1, double val2)
{

  if (val2==1.0)
    return(val1);
  else
    return((double)(val1+log(1-exp((long double)(val2-val1)))));
}  //End logsubtract




int delta (int a, int b)
  //Simple delta function used when calculating transprobs
{

 if (a==b)
     return(1);
 else
     return(0); 

}  //End delta




extern DATAFORMAT guess_dataformat (const char *filename, int len)
  //Tries to identify the format of a sequence alignment by the file extension:
  //nex for PAUP/NEXUS, phy for PHYLIP, fas for FASTA, and pir for PIR
{
  int i, choice;
  char getline[800];
  BOOL identified=FALSE;
  DATAFORMAT retval;
  
  i=0;
  
  while (filename[i]!='.' && i<len)
    i++;
  i++;

  switch (filename[i])
    {
    case 'n':
    case 'N':
      if ((filename[i+1]=='e' || filename[i+1]=='E') && (filename[i+2]=='x' || filename[i+2]=='X'))
	{
	  retval=NEXUS;
	  identified=TRUE;
	}
      break;
    case 'f':
    case 'F':
      if ((filename[i+1]=='a' || filename[i+1]=='A') && (filename[i+2]=='s' || filename[i+2]=='S'))
	{
	  retval=FASTA;
	  identified=TRUE;
	}
      break;
    case 'p':
    case 'P':
      switch (filename[i+1])
	{
	case 'i':
	case 'I':
	  if ((filename[i+2]=='r' || filename[i+2]=='R'))
	    {
	      retval=PIR;
	      identified=TRUE;
	    }
	  break;
	  
	case 'h':
	case 'H':
	  if ((filename[i+2]=='y' || filename[i+2]=='Y'))
	    {
	      retval=PHYLIP;
	      identified=TRUE;
	    }
	  break;      

	}
      break;
    }

  while (identified!=TRUE)
    {
      cout<<"Please enter the number of the format of the file "<<filename<<endl;
      cout<<"Choices are:  1) PIR/CODA\n"
	  <<"              2) NEXUS   \n"
	  <<"              3) PHYLIP (Interleaved only)\n"
	  <<"              4) FASTA\n";
      cin>>getline;
      choice=string_to_int(getline);
      
      switch (choice)
	{
	case 1:
	  retval=PIR;
	  identified=TRUE;
	  break;  
	case 2:
	  retval=NEXUS;
	  identified=TRUE;
	  break;  
	case 3:
	  retval=PHYLIP;
	  identified=TRUE;
	  break;  
	case 4:
	  retval=FASTA;
	  identified=TRUE;
	  break;      
	  
	}

    }

  return(retval);

}



void initialize_next_list_element (Non_gap_sites *&the_list, BOOL first)
{

  if (first==TRUE)
    the_list->last=0;
  else
    {
      the_list->next=new Non_gap_sites;
      the_list->next->last=the_list;
      the_list->next->next=0;
      the_list=the_list->next;
    }
}


double calc_pearson_correl (int num_elements, double *val1, double *val2) 
{
  int i;
  double cov=0, var1=0, var2=0, avg1=0, avg2=0;


  for(i=0; i<num_elements; i++) {
    avg1+=val1[i];
    avg2+=val2[i];
  }
  avg1/=num_elements;
  avg2/=num_elements;
  

  for(i=0; i<num_elements; i++) {
    var1+=pow((val1[i]-avg1),2);
    var2+=pow((val2[i]-avg2),2);
    cov+=(val1[i]-avg1)*(val2[i]-avg2);
  }

      
 
  return(cov/sqrt(var1*var2));

}


int rank_cmp(const void * a, const void* b)
{
  if (((Entry*)a)->value < ((Entry*)b)->value)
    return(-1);
  else if (((Entry*)a)->value==((Entry*)b)->value)
    return(0);
  else
    return(1);
}



double calc_spearman_correl (int num_elements, double *val1, double *val2) 
{
  int i, last_tie1, last_tie2, num_ties1=0, num_ties2=0;
  double *ranks1, *ranks2, tie_sum1=0, tie_sum2=0, n_sum=0, rank_diff, avg1, avg2;
  BOOL *dones;
  Entry *sorted_vals1, *sorted_vals2;
  Ties *ties1, *ties2;

  sorted_vals1=new Entry [num_elements];
  sorted_vals2=new Entry [num_elements];
  ties1=new Ties [num_elements/2];
  ties2=new Ties [num_elements/2];
  dones=new BOOL[num_elements];
  
  ranks1=new double [num_elements];
  ranks2=new double [num_elements];
 
  for(i=0; i<num_elements; i++) {
	sorted_vals1[i].value=val1[i];
	sorted_vals2[i].value=val2[i];
	sorted_vals1[i].position=i;
	sorted_vals2[i].position=i;
  }
     
  qsort(sorted_vals1, num_elements, sizeof(Entry), rank_cmp);
  qsort(sorted_vals2, num_elements, sizeof(Entry), rank_cmp);
  
  ranks1[sorted_vals1[0].position]=1;
  ranks2[sorted_vals2[0].position]=1;

  for(i=1; i<num_elements; i++) {
    if (sorted_vals1[i].value==sorted_vals1[i-1].value)
      ranks1[sorted_vals1[i].position]=ranks1[sorted_vals1[i-1].position];
    else
      ranks1[sorted_vals1[i].position]=i+1;
	
    if (sorted_vals2[i].value==sorted_vals2[i-1].value)
      ranks2[sorted_vals2[i].position]=ranks2[sorted_vals2[i-1].position];
    else
      ranks2[sorted_vals2[i].position]=i+1;
    
  }
       

	
  last_tie1=0;
  last_tie2=0;
  for(i=1; i<num_elements; i++) {
	
    if (sorted_vals1[i].value == sorted_vals1[i-1].value)
      {
	if (last_tie1 == 0) {
	  ties1[num_ties1].value=ranks1[sorted_vals1[i].position];
	  ties1[num_ties1].num=2;
	  num_ties1++;
	  last_tie1=1;
	}
	else {
	  ties1[num_ties1-1].num++;
	}	  
      }
    else {
      last_tie1=0;
    }
	
    if (sorted_vals2[i].value == sorted_vals2[i-1].value)
      {
	if (last_tie2 == 0) {
	  ties2[num_ties2].value=ranks2[sorted_vals2[i].position];
	  ties2[num_ties2].num=2;
	  num_ties2++;
	  last_tie2=1;
	}
	else {
	  ties2[num_ties2-1].num++;
	}	  
      }
    else {
      last_tie2=0;
    }
  }
  
    
  last_tie1=0;
  last_tie2=0;      
  for(i=0; i<num_elements; i++) {
    while((ties1[last_tie1].value < ranks1[sorted_vals1[i].position]) && (last_tie1<num_ties1)) {last_tie1++;}

    if (!(last_tie1==num_ties1))
      if ( ties1[last_tie1].value ==  ranks1[sorted_vals1[i].position])
	ranks1[sorted_vals1[i].position]+=(double)(ties1[last_tie1].num-1)/2.0;


    while(ties2[last_tie2].value < ranks2[sorted_vals2[i].position] && (last_tie2<num_ties2)) {last_tie2++;}
    if (!(last_tie2==num_ties2))
      if ( ties2[last_tie2].value ==  ranks2[sorted_vals2[i].position])
	ranks2[sorted_vals2[i].position]+=(double)(ties2[last_tie2].num-1)/2.0;
  }
  
 
  tie_sum1=0;
  for(i=0; i<num_ties1; i++) {
    tie_sum1 += (ties1[i].num*ties1[i].num*ties1[i].num)-
      ties1[i].num;
   // cout<<"Tie value 1: "<<ties1[i].value<<" Num: "<<ties1[i].num<<endl;
  }
  
  tie_sum2=0;
  for(i=0; i<num_ties2; i++) {
    tie_sum2 += (ties2[i].num*ties2[i].num*ties2[i].num)-
      ties2[i].num;
    //cout<<"Tie value 2: "<<ties2[i].value<<" Num: "<<ties2[i].num<<endl;
  }
  

  n_sum=(num_elements*num_elements*num_elements)-num_elements;
     
  rank_diff=0;

      
  for(i=0; i<num_elements; i++) 
    rank_diff+=(ranks1[i]-ranks2[i])*(ranks1[i]-ranks2[i]);

  delete[] sorted_vals1;
  delete[] sorted_vals2;
  delete[] ties1;
  delete[] ties2;
  delete[] dones;
  
  delete[] ranks1;
  delete[] ranks2;
 
 
  return((1.0-((6.0/n_sum)*(rank_diff+(1.0/12.0)*tie_sum1+(1.0/12.0)*tie_sum2)))/
	 (sqrt(1-tie_sum1/n_sum)*sqrt(1-tie_sum2/n_sum)));
	
 

}



extern double calc_shannon_entropy(int *counts, int num_items)
{
	int i, sum;
	double val, logval, total;

	total=0.0;
	sum=0;

	for(i=0; i<num_items; i++) {
		sum += counts[i];
	}

	for(i=0; i<num_items; i++) {
		if (counts[i] != 0) {
			val=(double)counts[i]/(double)sum;
			logval=-log(val)/LN_2;
			total+= val*logval;
		}
	}
	return(total);

}

double calc_mutual_information (double *indep_freqs, double **joint_freqs, int num_items)
{
	int i, j;
	double retval=0;
	
	for(i=0; i<num_items; i++) {
		for (j=i+1; j<num_items; j++) {
			if (joint_freqs[i][j] > 0.0)
				retval+=joint_freqs[i][j]* log (joint_freqs[i][j]/(indep_freqs[i]*indep_freqs[j]));
			//cout<<i<<"\t"<<joint_freqs[i][j]<<"/ "<<indep_freqs[i]<<"* "<<indep_freqs[j]<<" ="<<log (joint_freqs[i][j]/(indep_freqs[i]*indep_freqs[j]))<<" SF: "<<retval<<endl;
		}
	}
	return(retval);
	
}


double calc_variance(double *vals, int num_items)
{
    int i;
    double xsqr_sum, x_sum, ret_val;
    //Calculates the variance of the list in one pass--use the (n-1) correction to the scaling of the sums.
    
    xsqr_sum=0;
    x_sum=0;
    
    for(i=0; i<num_items; i++) {
        xsqr_sum += (vals[i]) * (vals[i]);
        x_sum += vals[i];
    }
    
    x_sum = x_sum/(double)num_items;
    
    
    ret_val=(xsqr_sum/((double)num_items-1.0) - (double)num_items*(x_sum*x_sum)/((double)num_items-1.0));
    return(ret_val);
}


extern void enum_permutations(int n, int &num_perm, int *ids, int **&perms)
{
	int i, j, k, x, num_temp, **temp_perms, *temp_ids;

	
	num_perm=recurse_factorial(n);

	perms = new int* [num_perm];
	for(i=0; i<num_perm; i++)
		perms[i]=new int [n];
	
	if (n > 2) {
		num_temp =n-1;

		temp_ids=new int[num_temp];		
		x=0;

		for(i=0; i<n; i++) {
			for(j=0; j<i; j++)
				temp_ids[j]=ids[j];
			for(j=i+1; j<n; j++)
				temp_ids[j-1]=ids[j];
			enum_permutations(n-1, num_temp, temp_ids, temp_perms);

			for(j=0; j<num_temp; j++) {
				for(k=0; k<n-1; k++)
					perms[x][k+1]=temp_perms[j][k];
				perms[x++][0]=ids[i];
			}


			for(j=0; j<num_temp; j++)
				delete[] temp_perms[j];
			delete[] temp_perms;
		}	
		delete[] temp_ids;
	}
	else 
	{
		perms[0][0]=ids[0];
		perms[0][1]=ids[1];
		perms[1][0]=ids[1];
		perms[1][1]=ids[0];
	}


}


extern int recurse_factorial(int i)
{
	if (i>1)
		return(i*recurse_factorial(i-1));
	else
		return(1);
}


extern long double float_factorial(int i)
{
	int j;
	long double val=1.0;

	if (i> 1)
	{
		for(j=2; j<=i; j++) 
			val = val * (long double)j;
		return(val);
	}
	else
		return(1.0);
}



extern long double log_factorial(int i)
{
	int j;
	long double val=0.0, last_val=-1.0;
	
	if (i> 1)
	{
		j=i;
		while((j>0) && (last_val < val)) { 
			last_val=val;
			val = val + log(j);
			j--;
		}
		return(val);
	}
	else
		return(0.0);
}




extern int get_dna_comp(int inbase)
{
	if (is_base(num_to_base(inbase))==TRUE) {
		switch(num_to_base(inbase)) {
		case 'C':
			return(readchar_to_base('G'));
			break;
		case 'G':
			return(readchar_to_base('C'));
			break;
		case 'A':
			return(readchar_to_base('T'));
			break;
		case 'T':
			return(readchar_to_base('A'));
			break;
		case 'Y':
			return(readchar_to_base('R'));
			break;
		case 'R':
			return(readchar_to_base('Y'));
			break;
		case '-':
			return(readchar_to_base('-'));
			break;
		default:
			return(readchar_to_base('N'));
			break;
		}
	}
	else
		return(readchar_to_base('N'));

}


extern Sequence_dataset * complement_dna (Sequence_dataset *orig_data)
{
	int i, j, k, ntaxa, nchars, *lens;
	Sequence_dataset *new_dataset;

	ntaxa=orig_data->Num_sequences();
	lens=new int [ntaxa];
	for(i=0; i<ntaxa; i++)
		lens[i]=(*orig_data)[i].Sequence_size();

	new_dataset=new Sequence_dataset(ntaxa, lens);

	for(i=0; i<ntaxa; i++) {
		(*new_dataset)[i].Assign_name((*orig_data)[i].Sequence_name());
		k=0;
		for(j=lens[i]-1; j>=0; j--) {
			(*new_dataset)[i].Assign_site(k++, get_dna_comp((*orig_data)[i][j]));	
		}


	}

	
	delete[] lens;
	return(new_dataset);
}


int n_choose_k(int n, int k)
//Note that this function assumes that n choose k can be represented as the 32-bit int
{
	int i, retval=1, large_denom, small_denom;


	if (k > (n-k)) {
		large_denom=k;
		small_denom=n-k;
	}
	else {
		large_denom=n-k;
		small_denom=k;
	}
	for(i=n; i>large_denom; i--)
		retval*=i;

	retval=retval/recurse_factorial(small_denom);
	return(retval);
}


long double float_n_choose_k(int n, int k)
{
	int i; 
	long double fretval=1.0, large_denom, small_denom;


	if (k > (n-k)) {
		large_denom=k;
		small_denom=n-k;
	}
	else {
		large_denom=n-k;
		small_denom=k;
	}
	for(i=n; i>large_denom; i--) {
		fretval=fretval*(long double)i;
	}

	fretval=fretval/float_factorial(small_denom);
	return(fretval);

}


extern long double table_prob(int vals[2][2])
{
	//Note--assumes that vals are a 2X2 matrix

	int sums[2];

	sums[0]=vals[0][0]+vals[1][0];
	sums[1]=vals[0][1]+vals[1][1];

	return((float_n_choose_k(sums[0], vals[0][0])*float_n_choose_k(sums[1], vals[0][1]))/
		float_n_choose_k(sums[0]+sums[1],vals[0][0]+vals[0][1]));
}


double calc_fishers_exact_p(int vals[2][2], BOOL left)
{
	int i, c_sums[2], r_sums[2], small_group, large_group, new_vals[2][2], max, min;
	long double real_pval, new_pval, sum_pval=0, prop_diff, new_prop_diff;

	real_pval=table_prob(vals);

	c_sums[0]=vals[0][0]+vals[1][0];
	c_sums[1]=vals[0][1]+vals[1][1];

	r_sums[0]=vals[0][0]+vals[0][1];
	r_sums[1]=vals[1][0]+vals[1][1];

	prop_diff=((double)vals[0][0]/(double)r_sums[0]-(double)vals[1][0]/(double)r_sums[1]);

	if (r_sums[0] > r_sums[1]) {
		small_group=r_sums[1];
		large_group=r_sums[0];
		
		//Group defs swap in calculation below
		prop_diff *= -1.0;
	}
	else {
		small_group=r_sums[0];
		large_group=r_sums[1];
	}

	if (small_group < c_sums[0])
		max=small_group;
	else
		max=c_sums[0];

	if (c_sums[0] > large_group)
		min = c_sums[0]-large_group;
	else
		min=0;

	for(i=min; i<=max; i++) {
		new_vals[0][0]=i;
		new_vals[0][1]=small_group-i;
		new_vals[1][0]=c_sums[0]-i;
		new_vals[1][1]=large_group-new_vals[1][0];

		new_pval=table_prob(new_vals);


		new_prop_diff=(((double)new_vals[0][0]/(double)(new_vals[0][0]+new_vals[0][1]))
			-((double)new_vals[1][0]/(double)(new_vals[1][0]+new_vals[1][1])));

		if (left == TRUE) {
			if (new_prop_diff >= prop_diff)
				sum_pval+=new_pval;
		}
		else
			if (new_prop_diff <= prop_diff)
				sum_pval+=new_pval;
	}


	return((double)sum_pval);
}


extern Sequence_dataset * make_writable_dupl_dataset(Sequence_dataset *orig_data)
{
	int i, j, k, ntaxa, nchars, *lens;
	Sequence_dataset *new_dataset;

	ntaxa=orig_data->Num_sequences();
	lens=new int [ntaxa];
	for(i=0; i<ntaxa; i++)
		lens[i]=(*orig_data)[i].Sequence_size();

	new_dataset=new Sequence_dataset(ntaxa, lens);

	for(i=0; i<ntaxa; i++) {
		(*new_dataset)[i].Assign_name((*orig_data)[i].Sequence_name());
		for(j=0; j<(*orig_data)[i].Sequence_size(); j++) {
			(*new_dataset)[i].Assign_site(j, 
					loss_state_to_dupl(get_observable_dupl_state(dupl_to_loss_state((*orig_data)[i][j]))));

		}


	}

	
	delete[] lens;
	return(new_dataset);

}
