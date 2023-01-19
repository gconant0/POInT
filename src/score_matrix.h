#ifndef ___SCORE_MATRIX_H___
#define ___SCORE_MATRIX_H___

#include <fstream>
#include "gen_dna_funcs.h"

using namespace::std;

//A set of classes to contain/obtain amino acid/nucleotide substitution matrices (such as BLOSUM)
//from various sources

class Matrix_line 
{
public:
	Matrix_line();
	Matrix_line(int sz);
	double operator[](int element);
	void set_value (int element, double val);
	~Matrix_line();

protected:
	int size;
	double *line;

};

class Score_matrix 
{
public:
	Score_matrix();
	Score_matrix(int mat_size);
	Matrix_line& operator[](int element);
	int get_matrix_size()                {return(matrix_size);};
	char* get_matrix_name()    {return(matrix_name);};
	void set_entry(int i, int j, double val);
	~Score_matrix();

protected:
	int matrix_size;
	char matrix_name[40];
	Matrix_line **matrix;
};


class Nuc_matrix : public Score_matrix
{
public:	
	Nuc_matrix ();
	Nuc_matrix (double match, double mismatch);
};

class Nuc_arb_matrix : public Nuc_matrix 
{
public:
	Nuc_arb_matrix ();
	Nuc_arb_matrix (double **the_matrix);
};


class AA_matrix : public Score_matrix
{
public:
	AA_matrix();
	
protected:

};


class Simple_AA_matrix : public AA_matrix
{
public:
	Simple_AA_matrix();
	Simple_AA_matrix(double match, double mismatch);
};


class Arb_AA_matrix : public Simple_AA_matrix
{
public:
	Arb_AA_matrix ();
	Arb_AA_matrix (double **the_matrix);

};


class File_AA_matrix : public Simple_AA_matrix
{
public:
	File_AA_matrix();
	File_AA_matrix(char *filename);
protected:
	char code[24];
};


class BLOSUM_62_matrix : public AA_matrix
{
public:
	BLOSUM_62_matrix();

};



#endif
