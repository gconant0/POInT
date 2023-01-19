#include "score_matrix.h"
#include "gen_dna_funcs.h"
#include <iostream>
#include "string.h"

using namespace std;

Matrix_line::Matrix_line()
{
	size=20;
	line=new double[size];
}


Matrix_line::Matrix_line(int sz)
{
	size=sz;
	line=new double[size];
}


double Matrix_line::operator[](int element)
{
	if (element >= size) return(line[size-1]);
	else return(line[element]);
}


void Matrix_line::set_value (int element, double val)
{
	line[element]=val;

}


Matrix_line::~Matrix_line()
{
	delete[] line;
}	

Score_matrix::Score_matrix() 
{
	int i, j;

	strcpy(matrix_name, "Unnamed AA matrix");
	matrix_size=20;

	matrix=new Matrix_line* [matrix_size];
	for(i=0; i<matrix_size; i++)
		matrix[i]=new Matrix_line(matrix_size);
}


Score_matrix::Score_matrix(int mat_size)
{
	int i, j;

	matrix_size=mat_size;

	strcpy(matrix_name, "Unnamed matrix");

	matrix= new Matrix_line* [matrix_size];


	for(i=0; i<matrix_size; i++)
		matrix[i]=new Matrix_line (matrix_size);

}


Matrix_line& Score_matrix::operator[](int element)
{
	if (element >= matrix_size) return((*matrix[matrix_size-1]));
	else   return((*matrix[element])); 
}


void Score_matrix::set_entry(int i, int j, double val)
{
	matrix[i]->set_value(j, val);
}



Score_matrix::~Score_matrix()
{
	int i;

	for(i=0; i<matrix_size; i++)
		delete matrix[i];
	delete[] matrix;

}


Nuc_matrix::Nuc_matrix() : Score_matrix (5)
{
	int i, j;

  
	strcpy(matrix_name, "Default Nuc. matrix: 5/-4");

	for(i=0; i<matrix_size; i++) {
		for(j=0; j<matrix_size; j++)
		{
			if ((i==(matrix_size-1)) || (j == (matrix_size-1))) 
				matrix[i]->set_value(j, 5.0);
			else {
				if (i==j)
					matrix[i]->set_value(j, 5.0);
				else
					matrix[i]->set_value(j, -4.0);
			}
		}
	}
}


Nuc_matrix::Nuc_matrix(double match, double mismatch) : Score_matrix (5)
{
	int i, j;
	char param_string[10];

   	strcpy(matrix_name, "Nuc. matrix: ");
	double_to_string(param_string, 9, 3, match);
	strcat(matrix_name, param_string);
	strcat(matrix_name, ", ");
	double_to_string(param_string, 9, 3, mismatch);
	strcat(matrix_name, param_string);

	for(i=0; i<matrix_size; i++) {
       for(j=0; j<matrix_size; j++)
		{
			if ((i==(matrix_size-1)) || (j == (matrix_size-1))) 
				matrix[i]->set_value(j, 5.0);
			else {
				if (i==j)
					matrix[i]->set_value(j, match);
				else
					matrix[i]->set_value(j, mismatch);
			}
		}
	}
}



Nuc_arb_matrix::Nuc_arb_matrix () : Nuc_matrix ()
{
	cerr<<"Error: Call to arbitrary nucleotide matrix without input matrix\n";
}



Nuc_arb_matrix::Nuc_arb_matrix (double **the_matrix) : Nuc_matrix ()
{
	int i, j;

	strcpy(matrix_name, "Nuc. matrix w/ input matrix");

	for(i=0; i<matrix_size; i++)
		for(j=0; j<matrix_size; j++)
			matrix[i]->set_value(j, the_matrix[i][j]);

}


AA_matrix::AA_matrix() : Score_matrix(20)
{

}

Simple_AA_matrix::Simple_AA_matrix() : AA_matrix ()
{
	int i, j;

	strcpy(matrix_name, "Default AA matrix: 1/-1");

	for(i=0; i<matrix_size; i++)
		for (j=0; j<matrix_size; j++)
		{
			if (i==j)
				matrix[i]->set_value(j, 1);
			else
				matrix[i]->set_value(j, -1);

		}
}


Simple_AA_matrix::Simple_AA_matrix(double match, double mismatch) : AA_matrix ()
{
	int i, j;
	char param_string[10];

	strcpy(matrix_name, "AA matrix: ");
	double_to_string(param_string, 9, 3, match);
	strcat(matrix_name, param_string);
	strcat(matrix_name, ", ");
	double_to_string(param_string, 9, 3, mismatch);
	strcat(matrix_name, param_string);

	
	for(i=0; i<matrix_size; i++)
		for (j=0; j<matrix_size; j++)
		{
			if (i==j)
				matrix[i]->set_value(j, match);
			else
				matrix[i]->set_value(j, mismatch);

		}

}


Arb_AA_matrix::Arb_AA_matrix () : Simple_AA_matrix ()
{
	cerr<<"Error: Call to arbitrary amino acid matrix without input matrix\n";
}



Arb_AA_matrix::Arb_AA_matrix (double **the_matrix) : Simple_AA_matrix ()
{
	int i, j;

	strcpy(matrix_name, "AA matrix w/ Arb. elements ");
	

	for(i=0; i<matrix_size; i++)
		for(j=0; j<matrix_size; j++)
			matrix[i]->set_value(j, the_matrix[i][j]);

}

File_AA_matrix::File_AA_matrix() : Simple_AA_matrix()
{
	cerr<<"Error: Call to file-based AA matrix without filename\n";
}


File_AA_matrix::File_AA_matrix(char *filename) : Simple_AA_matrix() 
{
	char inchar, line[300];
	int i, j;
	double innum;
	ifstream infile;

	infile.open(filename);
	
	if (infile.fail()) {
		cerr<<"Can't find file "<<filename<<"--Using default matrix: Match =0, Mis-match=-1\n";
	}
	else {
		strcpy(matrix_name, "AA matrix from file ");
		strcat(matrix_name, filename);


		for (i=0; i<6; i++) {
			infile.getline(line, 299);
		}

		infile.get(inchar);
		while(inchar == 32)
			infile.get(inchar);

		for(i=0; i<24; i++)
		{
			code[i]=inchar;
			infile.get(inchar);
			infile.get(inchar);    
			infile.get(inchar);
		}

		for(i=0; i<24; i++) {
			for(j=0; j<24; j++)
			{
				infile>>innum;
			
				if ((readchar_to_aa(code[i])<20) && (readchar_to_aa(code[j])<20)) {
				
					matrix[readchar_to_aa(code[i])]->set_value(readchar_to_aa(code[j]), innum);
					
				}
			}
			
			infile>>inchar;
		}
	}


	

}

BLOSUM_62_matrix::BLOSUM_62_matrix() : AA_matrix ()
{
int i,j ;
    
	strcpy(matrix_name, "Built-in BLOSUM62 matrix");

matrix[0]->set_value(0, 4);
matrix[0]->set_value(1, 0);
matrix[0]->set_value(2, -2);
matrix[0]->set_value(3, -1);
matrix[0]->set_value(4, -2);
matrix[0]->set_value(5, 0);
matrix[0]->set_value(6, -2);
matrix[0]->set_value(7, -1);
matrix[0]->set_value(8, -1);
matrix[0]->set_value(9, -1);
matrix[0]->set_value(10, -1);
matrix[0]->set_value(11, -2);
matrix[0]->set_value(12, -1);
matrix[0]->set_value(13, -1);
matrix[0]->set_value(14, -1);
matrix[0]->set_value(15, 1);
matrix[0]->set_value(16, 0);
matrix[0]->set_value(17, 0);
matrix[0]->set_value(18, -3);
matrix[0]->set_value(19, -2);
matrix[1]->set_value(0, 0);
matrix[1]->set_value(1, 9);
matrix[1]->set_value(2, -3);
matrix[1]->set_value(3, -4);
matrix[1]->set_value(4, -2);
matrix[1]->set_value(5, -3);
matrix[1]->set_value(6, -3);
matrix[1]->set_value(7, -1);
matrix[1]->set_value(8, -3);
matrix[1]->set_value(9, -1);
matrix[1]->set_value(10, -1);
matrix[1]->set_value(11, -3);
matrix[1]->set_value(12, -3);
matrix[1]->set_value(13, -3);
matrix[1]->set_value(14, -3);
matrix[1]->set_value(15, -1);
matrix[1]->set_value(16, -1);
matrix[1]->set_value(17, -1);
matrix[1]->set_value(18, -2);
matrix[1]->set_value(19, -2);
matrix[2]->set_value(0, -2);
matrix[2]->set_value(1, -3);
matrix[2]->set_value(2, 6);
matrix[2]->set_value(3, 2);
matrix[2]->set_value(4, -3);
matrix[2]->set_value(5, -1);
matrix[2]->set_value(6, -1);
matrix[2]->set_value(7, -3);
matrix[2]->set_value(8, -1);
matrix[2]->set_value(9, -4);
matrix[2]->set_value(10, -3);
matrix[2]->set_value(11, 1);
matrix[2]->set_value(12, -1);
matrix[2]->set_value(13, 0);
matrix[2]->set_value(14, -2);
matrix[2]->set_value(15, 0);
matrix[2]->set_value(16, -1);
matrix[2]->set_value(17, -3);
matrix[2]->set_value(18, -4);
matrix[2]->set_value(19, -3);
matrix[3]->set_value(0, -1);
matrix[3]->set_value(1, -4);
matrix[3]->set_value(2, 2);
matrix[3]->set_value(3, 5);
matrix[3]->set_value(4, -3);
matrix[3]->set_value(5, -2);
matrix[3]->set_value(6, 0);
matrix[3]->set_value(7, -3);
matrix[3]->set_value(8, 1);
matrix[3]->set_value(9, -3);
matrix[3]->set_value(10, -2);
matrix[3]->set_value(11, 0);
matrix[3]->set_value(12, -1);
matrix[3]->set_value(13, 2);
matrix[3]->set_value(14, 0);
matrix[3]->set_value(15, 0);
matrix[3]->set_value(16, -1);
matrix[3]->set_value(17, -2);
matrix[3]->set_value(18, -3);
matrix[3]->set_value(19, -2);
matrix[4]->set_value(0, -2);
matrix[4]->set_value(1, -2);
matrix[4]->set_value(2, -3);
matrix[4]->set_value(3, -3);
matrix[4]->set_value(4, 6);
matrix[4]->set_value(5, -3);
matrix[4]->set_value(6, -1);
matrix[4]->set_value(7, 0);
matrix[4]->set_value(8, -3);
matrix[4]->set_value(9, 0);
matrix[4]->set_value(10, 0);
matrix[4]->set_value(11, -3);
matrix[4]->set_value(12, -4);
matrix[4]->set_value(13, -3);
matrix[4]->set_value(14, -3);
matrix[4]->set_value(15, -2);
matrix[4]->set_value(16, -2);
matrix[4]->set_value(17, -1);
matrix[4]->set_value(18, 1);
matrix[4]->set_value(19, 3);
matrix[5]->set_value(0, 0);
matrix[5]->set_value(1, -3);
matrix[5]->set_value(2, -1);
matrix[5]->set_value(3, -2);
matrix[5]->set_value(4, -3);
matrix[5]->set_value(5, 6);
matrix[5]->set_value(6, -2);
matrix[5]->set_value(7, -4);
matrix[5]->set_value(8, -2);
matrix[5]->set_value(9, -4);
matrix[5]->set_value(10, -3);
matrix[5]->set_value(11, 0);
matrix[5]->set_value(12, -2);
matrix[5]->set_value(13, -2);
matrix[5]->set_value(14, -2);
matrix[5]->set_value(15, 0);
matrix[5]->set_value(16, -2);
matrix[5]->set_value(17, -3);
matrix[5]->set_value(18, -2);
matrix[5]->set_value(19, -3);
matrix[6]->set_value(0, -2);
matrix[6]->set_value(1, -3);
matrix[6]->set_value(2, -1);
matrix[6]->set_value(3, 0);
matrix[6]->set_value(4, -1);
matrix[6]->set_value(5, -2);
matrix[6]->set_value(6, 8);
matrix[6]->set_value(7, -3);
matrix[6]->set_value(8, -1);
matrix[6]->set_value(9, -3);
matrix[6]->set_value(10, -2);
matrix[6]->set_value(11, 1);
matrix[6]->set_value(12, -2);
matrix[6]->set_value(13, 0);
matrix[6]->set_value(14, 0);
matrix[6]->set_value(15, -1);
matrix[6]->set_value(16, -2);
matrix[6]->set_value(17, -3);
matrix[6]->set_value(18, -2);
matrix[6]->set_value(19, 2);
matrix[7]->set_value(0, -1);
matrix[7]->set_value(1, -1);
matrix[7]->set_value(2, -3);
matrix[7]->set_value(3, -3);
matrix[7]->set_value(4, 0);
matrix[7]->set_value(5, -4);
matrix[7]->set_value(6, -3);
matrix[7]->set_value(7, 4);
matrix[7]->set_value(8, -3);
matrix[7]->set_value(9, 2);
matrix[7]->set_value(10, 1);
matrix[7]->set_value(11, -3);
matrix[7]->set_value(12, -3);
matrix[7]->set_value(13, -3);
matrix[7]->set_value(14, -3);
matrix[7]->set_value(15, -2);
matrix[7]->set_value(16, -1);
matrix[7]->set_value(17, 3);
matrix[7]->set_value(18, -3);
matrix[7]->set_value(19, -1);
matrix[8]->set_value(0, -1);
matrix[8]->set_value(1, -3);
matrix[8]->set_value(2, -1);
matrix[8]->set_value(3, 1);
matrix[8]->set_value(4, -3);
matrix[8]->set_value(5, -2);
matrix[8]->set_value(6, -1);
matrix[8]->set_value(7, -3);
matrix[8]->set_value(8, 5);
matrix[8]->set_value(9, -2);
matrix[8]->set_value(10, -1);
matrix[8]->set_value(11, 0);
matrix[8]->set_value(12, -1);
matrix[8]->set_value(13, 1);
matrix[8]->set_value(14, 2);
matrix[8]->set_value(15, 0);
matrix[8]->set_value(16, -1);
matrix[8]->set_value(17, -2);
matrix[8]->set_value(18, -3);
matrix[8]->set_value(19, -2);
matrix[9]->set_value(0, -1);
matrix[9]->set_value(1, -1);
matrix[9]->set_value(2, -4);
matrix[9]->set_value(3, -3);
matrix[9]->set_value(4, 0);
matrix[9]->set_value(5, -4);
matrix[9]->set_value(6, -3);
matrix[9]->set_value(7, 2);
matrix[9]->set_value(8, -2);
matrix[9]->set_value(9, 4);
matrix[9]->set_value(10, 2);
matrix[9]->set_value(11, -3);
matrix[9]->set_value(12, -3);
matrix[9]->set_value(13, -2);
matrix[9]->set_value(14, -2);
matrix[9]->set_value(15, -2);
matrix[9]->set_value(16, -1);
matrix[9]->set_value(17, 1);
matrix[9]->set_value(18, -2);
matrix[9]->set_value(19, -1);
matrix[10]->set_value(0, -1);
matrix[10]->set_value(1, -1);
matrix[10]->set_value(2, -3);
matrix[10]->set_value(3, -2);
matrix[10]->set_value(4, 0);
matrix[10]->set_value(5, -3);
matrix[10]->set_value(6, -2);
matrix[10]->set_value(7, 1);
matrix[10]->set_value(8, -1);
matrix[10]->set_value(9, 2);
matrix[10]->set_value(10, 5);
matrix[10]->set_value(11, -2);
matrix[10]->set_value(12, -2);
matrix[10]->set_value(13, 0);
matrix[10]->set_value(14, -1);
matrix[10]->set_value(15, -1);
matrix[10]->set_value(16, -1);
matrix[10]->set_value(17, 1);
matrix[10]->set_value(18, -1);
matrix[10]->set_value(19, -1);
matrix[11]->set_value(0, -2);
matrix[11]->set_value(1, -3);
matrix[11]->set_value(2, 1);
matrix[11]->set_value(3, 0);
matrix[11]->set_value(4, -3);
matrix[11]->set_value(5, 0);
matrix[11]->set_value(6, 1);
matrix[11]->set_value(7, -3);
matrix[11]->set_value(8, 0);
matrix[11]->set_value(9, -3);
matrix[11]->set_value(10, -2);
matrix[11]->set_value(11, 6);
matrix[11]->set_value(12, -2);
matrix[11]->set_value(13, 0);
matrix[11]->set_value(14, 0);
matrix[11]->set_value(15, 1);
matrix[11]->set_value(16, 0);
matrix[11]->set_value(17, -3);
matrix[11]->set_value(18, -4);
matrix[11]->set_value(19, -2);
matrix[12]->set_value(0, -1);
matrix[12]->set_value(1, -3);
matrix[12]->set_value(2, -1);
matrix[12]->set_value(3, -1);
matrix[12]->set_value(4, -4);
matrix[12]->set_value(5, -2);
matrix[12]->set_value(6, -2);
matrix[12]->set_value(7, -3);
matrix[12]->set_value(8, -1);
matrix[12]->set_value(9, -3);
matrix[12]->set_value(10, -2);
matrix[12]->set_value(11, -2);
matrix[12]->set_value(12, 7);
matrix[12]->set_value(13, -1);
matrix[12]->set_value(14, -2);
matrix[12]->set_value(15, -1);
matrix[12]->set_value(16, -1);
matrix[12]->set_value(17, -2);
matrix[12]->set_value(18, -4);
matrix[12]->set_value(19, -3);
matrix[13]->set_value(0, -1);
matrix[13]->set_value(1, -3);
matrix[13]->set_value(2, 0);
matrix[13]->set_value(3, 2);
matrix[13]->set_value(4, -3);
matrix[13]->set_value(5, -2);
matrix[13]->set_value(6, 0);
matrix[13]->set_value(7, -3);
matrix[13]->set_value(8, 1);
matrix[13]->set_value(9, -2);
matrix[13]->set_value(10, 0);
matrix[13]->set_value(11, 0);
matrix[13]->set_value(12, -1);
matrix[13]->set_value(13, 5);
matrix[13]->set_value(14, 1);
matrix[13]->set_value(15, 0);
matrix[13]->set_value(16, -1);
matrix[13]->set_value(17, -2);
matrix[13]->set_value(18, -2);
matrix[13]->set_value(19, -1);
matrix[14]->set_value(0, -1);
matrix[14]->set_value(1, -3);
matrix[14]->set_value(2, -2);
matrix[14]->set_value(3, 0);
matrix[14]->set_value(4, -3);
matrix[14]->set_value(5, -2);
matrix[14]->set_value(6, 0);
matrix[14]->set_value(7, -3);
matrix[14]->set_value(8, 2);
matrix[14]->set_value(9, -2);
matrix[14]->set_value(10, -1);
matrix[14]->set_value(11, 0);
matrix[14]->set_value(12, -2);
matrix[14]->set_value(13, 1);
matrix[14]->set_value(14, 5);
matrix[14]->set_value(15, -1);
matrix[14]->set_value(16, -1);
matrix[14]->set_value(17, -3);
matrix[14]->set_value(18, -3);
matrix[14]->set_value(19, -2);
matrix[15]->set_value(0, 1);
matrix[15]->set_value(1, -1);
matrix[15]->set_value(2, 0);
matrix[15]->set_value(3, 0);
matrix[15]->set_value(4, -2);
matrix[15]->set_value(5, 0);
matrix[15]->set_value(6, -1);
matrix[15]->set_value(7, -2);
matrix[15]->set_value(8, 0);
matrix[15]->set_value(9, -2);
matrix[15]->set_value(10, -1);
matrix[15]->set_value(11, 1);
matrix[15]->set_value(12, -1);
matrix[15]->set_value(13, 0);
matrix[15]->set_value(14, -1);
matrix[15]->set_value(15, 4);
matrix[15]->set_value(16, 1);
matrix[15]->set_value(17, -2);
matrix[15]->set_value(18, -3);
matrix[15]->set_value(19, -2);
matrix[16]->set_value(0, 0);
matrix[16]->set_value(1, -1);
matrix[16]->set_value(2, -1);
matrix[16]->set_value(3, -1);
matrix[16]->set_value(4, -2);
matrix[16]->set_value(5, -2);
matrix[16]->set_value(6, -2);
matrix[16]->set_value(7, -1);
matrix[16]->set_value(8, -1);
matrix[16]->set_value(9, -1);
matrix[16]->set_value(10, -1);
matrix[16]->set_value(11, 0);
matrix[16]->set_value(12, -1);
matrix[16]->set_value(13, -1);
matrix[16]->set_value(14, -1);
matrix[16]->set_value(15, 1);
matrix[16]->set_value(16, 5);
matrix[16]->set_value(17, 0);
matrix[16]->set_value(18, -2);
matrix[16]->set_value(19, -2);
matrix[17]->set_value(0, 0);
matrix[17]->set_value(1, -1);
matrix[17]->set_value(2, -3);
matrix[17]->set_value(3, -2);
matrix[17]->set_value(4, -1);
matrix[17]->set_value(5, -3);
matrix[17]->set_value(6, -3);
matrix[17]->set_value(7, 3);
matrix[17]->set_value(8, -2);
matrix[17]->set_value(9, 1);
matrix[17]->set_value(10, 1);
matrix[17]->set_value(11, -3);
matrix[17]->set_value(12, -2);
matrix[17]->set_value(13, -2);
matrix[17]->set_value(14, -3);
matrix[17]->set_value(15, -2);
matrix[17]->set_value(16, 0);
matrix[17]->set_value(17, 4);
matrix[17]->set_value(18, -3);
matrix[17]->set_value(19, -1);
matrix[18]->set_value(0, -3);
matrix[18]->set_value(1, -2);
matrix[18]->set_value(2, -4);
matrix[18]->set_value(3, -3);
matrix[18]->set_value(4, 1);
matrix[18]->set_value(5, -2);
matrix[18]->set_value(6, -2);
matrix[18]->set_value(7, -3);
matrix[18]->set_value(8, -3);
matrix[18]->set_value(9, -2);
matrix[18]->set_value(10, -1);
matrix[18]->set_value(11, -4);
matrix[18]->set_value(12, -4);
matrix[18]->set_value(13, -2);
matrix[18]->set_value(14, -3);
matrix[18]->set_value(15, -3);
matrix[18]->set_value(16, -2);
matrix[18]->set_value(17, -3);
matrix[18]->set_value(18, 11);
matrix[18]->set_value(19, 2);
matrix[19]->set_value(0, -2);
matrix[19]->set_value(1, -2);
matrix[19]->set_value(2, -3);
matrix[19]->set_value(3, -2);
matrix[19]->set_value(4, 3);
matrix[19]->set_value(5, -3);
matrix[19]->set_value(6, 2);
matrix[19]->set_value(7, -1);
matrix[19]->set_value(8, -2);
matrix[19]->set_value(9, -1);
matrix[19]->set_value(10, -1);
matrix[19]->set_value(11, -2);
matrix[19]->set_value(12, -3);
matrix[19]->set_value(13, -1);
matrix[19]->set_value(14, -2);
matrix[19]->set_value(15, -2);
matrix[19]->set_value(16, -2);
matrix[19]->set_value(17, -1);
matrix[19]->set_value(18, 2);
matrix[19]->set_value(19, 7);


}
