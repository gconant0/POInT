//Copyright 1999-2002 Gavin Conant

#ifndef ___READ_SEQ_H___
#define ___READ_SEQ_H___

#include <fstream>
#include <iostream>
using namespace std;

#include "gen_dna_funcs.h"



struct sequence {
  int site;
  sequence *next_site;
};


BOOL test_interleave(char *filename);


class Read_Sequence 
{
public:
	//Functions
	Read_Sequence();
	virtual Sequence_dataset * get_dataset(int &num_taxa, int &num_sites, 
			const char *filename, BOOL is_aa_seq)=0;
	virtual Sequence_dataset * get_dataset(int &num_taxa, int &num_sites, 
			const char *filename, DATATYPE type)=0;
	virtual Sequence_dataset  * get_dataset(int &num_taxa, int &num_sites, DATATYPE type)=0;
	int get_num_taxa();
	int get_num_chars();
	BOOL protein_data();
	DATATYPE get_type()  {return(ctype);};
    virtual ~Read_Sequence() {};

protected:
	int size, ntaxa, *sizes;
	char last_char, inbase;
	BOOL std_in, finished;
	istream *var_stream;
	ifstream file_stream;
	BOOL aa_data;
	DATATYPE ctype;
	sequence *in_seq, *in_seq_start, *temp;
	List_Molecule_Sequence *build_taxa, *taxa_start, *temp_taxa;
	int (*convert_char)(char);
	BOOL(*valid_input)(char);

	//Functions
	void get_sequence_array_size_known(Molecule_Sequence *new_sequence);
	void get_sequence_array(Molecule_Sequence *&new_sequence);
	void set_datatype(DATATYPE cdata);
	void get_next_list_element(List_Molecule_Sequence *&thelist);
	void get_next_sequence_element(sequence *&thelist);
	BOOL read_finished();
	char next_char();
};



class Read_PIR : public Read_Sequence 
{
public:
    Read_PIR(): Read_Sequence() {};
    ~Read_PIR();
	Sequence_dataset * get_dataset(int &num_taxa, int &num_sites, 
			const char *filename, BOOL is_aa_seq);
	Sequence_dataset * get_dataset(int &num_taxa, int &num_sites, 
			const char *filename, DATATYPE type);
	Sequence_dataset * get_dataset(int &num_taxa, int &num_sites, DATATYPE type);

};



class Read_Nexus : public Read_Sequence 
{
public:
    Read_Nexus(): Read_Sequence() {};
    ~Read_Nexus();
	Sequence_dataset * get_dataset(int &num_taxa, int &num_sites, 
			const char *filename, BOOL is_aa_seq);
	Sequence_dataset * get_dataset(int &num_taxa, int &num_sites, 
			const char *filename, DATATYPE type);
	Sequence_dataset * get_dataset(int &num_taxa, int &num_sites, DATATYPE type);
};



class Read_Phylip_interleave : public Read_Sequence 
{
public:
    Read_Phylip_interleave(): Read_Sequence() {};
    ~Read_Phylip_interleave();
	Sequence_dataset * get_dataset(int &num_taxa, int &num_sites, 
			const char *filename, BOOL is_aa_seq);
	Sequence_dataset * get_dataset(int &num_taxa, int &num_sites, 
			const char *filename, DATATYPE type);
	Sequence_dataset * get_dataset(int &num_taxa, int &num_sites, DATATYPE type);
};



class Read_Phylip_noninterleave : public Read_Sequence 
{
public:
    Read_Phylip_noninterleave(): Read_Sequence() {};
    ~Read_Phylip_noninterleave();
	Sequence_dataset * get_dataset(int &num_taxa, int &num_sites, 
				 const char *filename, BOOL is_aa_seq);
	Sequence_dataset * get_dataset(int &num_taxa, int &num_sites, 
				 const char *filename, DATATYPE type);
	Sequence_dataset * get_dataset(int &num_taxa, int &num_sites, DATATYPE type);
};




class Read_FASTA : public Read_Sequence
{
public:
    Read_FASTA() :Read_Sequence() {finished=FALSE; started=FALSE;}
	Sequence_dataset * get_dataset(int &num_taxa, int &num_sites, 
		   const char *filename, BOOL is_aa_seq);
	Sequence_dataset * get_dataset(int &num_taxa, int &num_sites, 
		   const char *filename, DATATYPE type);
	Sequence_dataset * get_dataset(int &num_taxa, int &num_sites, DATATYPE type);
    Sequence_dataset * get_data_subset(int size, string filename, int &num_taxa, DATATYPE type);
    BOOL finished_reading()  {return(finished);};
	~Read_FASTA();
protected:
    BOOL finished, started;
};




#endif




