//Copyright 1999-2002 Gavin Conant

#ifndef ___WRITE_SEQ_H___
#define ___WRITE_SEQ_H___


#include <fstream>
#include <iostream>
#include "gen_dna_funcs.h"

using namespace::std;

class Write_Sequence 
{
 public:
  //Functions
  Write_Sequence();
  Write_Sequence(const char *output_file, DATATYPE cdata);
	Write_Sequence(const char *output_file, DATATYPE cdata, BOOL append);
  virtual void write_dataset(int ntaxa, int nchars, Sequence_dataset *dataset);
  virtual void write_dataset(int ntaxa, Sequence_dataset *dataset);
  virtual void write_sequence(int nchars, Molecule_Sequence &seq_data)=0;
  BOOL fail();

 protected:
  DATATYPE curr_data;
  ofstream outfile;
  char (*convert_int)(int);
  BOOL failed;

  //Functions
  void set_output_type(DATATYPE cdata);
};



class Write_PIR : public Write_Sequence 
{
 public:
  Write_PIR();
  Write_PIR(const char *output_file, DATATYPE cdata);
  void write_sequence(int nchars, Molecule_Sequence &seq_data);
};



class Write_Nexus : public Write_Sequence 
{
 public: 
  Write_Nexus();
  Write_Nexus(const char *output_file, DATATYPE cdata);
  Write_Nexus(const char *output_file, DATATYPE cdata, char** data_comments, int c_lines);
  void write_dataset(int ntaxa, int nchars, Sequence_dataset *dataset);
  void write_dataset(int ntaxa, Sequence_dataset *dataset);
  void write_sequence(int nchars, Molecule_Sequence &seq_data);
  ~Write_Nexus();
 protected:
  int comment_lines;
  char **comments;
  BOOL has_comments;

};



class Write_Phylip_interleave : public Write_Sequence 
{
 public:
  Write_Phylip_interleave();
  Write_Phylip_interleave(const char *output_file, DATATYPE cdata);
  void write_dataset(int ntaxa, int nchars, Sequence_dataset *dataset); 
  void write_dataset(int ntaxa, Sequence_dataset *dataset);
  void write_sequence(int nchars, Molecule_Sequence &seq_data);
   
};


class Write_Phylip_noninterleave : public Write_Sequence 
{
 public:
  Write_Phylip_noninterleave();
  Write_Phylip_noninterleave(const char *output_file, DATATYPE cdata);
  void write_dataset(int ntaxa, int nchars, Sequence_dataset *dataset); 
  void write_dataset(int ntaxa, Sequence_dataset *dataset);
  void write_sequence(int nchars, Molecule_Sequence &seq_data);

};


class Write_FASTA : public Write_Sequence 
{
 public:
  Write_FASTA();
  Write_FASTA(const char *output_file, DATATYPE cdata);
  Write_FASTA(const char *output_file, DATATYPE cdata, BOOL append);
  void write_sequence(int nchars, Molecule_Sequence &seq_data);
};

#endif









