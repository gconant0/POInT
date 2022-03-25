//Gavin Conant  10/28/1999
//Class for linear algebra routines needed to calculate
//transition probablities for codon model of evolution
//Uses LAPACK library, either in subdirectory /lin_alg
//or machine local
#include <iostream>
#include <math.h>
#include <iomanip>
#include "maxlike.h"

using namespace::std;

#ifndef ___LIN_ALG_H___
#include "lin_alg.h"
#endif

Generic_Linear_Algebra::Generic_Linear_Algebra(int matrix_s)
//Constructor for the Linear_Algebra class--initializes the 
//variables to be passed to the fortran LAPACK libraries
//These initializations are for 64X64 matrices
{
	matrix_size=new long int;
	
	*matrix_size=matrix_s;
	
	work_size=new long int;
	*work_size=6*(*matrix_size);
	
	leadingd= new long int;
	leadingdr= new long int;
	*leadingd=*leadingdr=(*matrix_size);
	
	leadingdl= new long int;
	lwork=new long int;
	*leadingdl=1;
	*lwork=1;
	
	infor=new long int;
	
	arg_len=new long int;
	*arg_len=1;
	
	arg_len2=new long int;
	*arg_len2=1;
	
	char_1=new char;
	char_2= new char;
	
	pivots= new long int [*matrix_size];
	
	workspace= new double [*work_size];
	pass_matrix= new double [*matrix_size*(*matrix_size)];
	pass_matrix_b= new double [*matrix_size*(*matrix_size)];
	pass_matrix_c= new double [*matrix_size*(*matrix_size)];
	store_eigen_v= new double [*matrix_size*(*matrix_size)];
	store_c_eigen_v= new double [*matrix_size*(*matrix_size)];
	store_inv_eigen_v= new double [*matrix_size*(*matrix_size)];
	store_eigen_values= new double [(*matrix_size)];
	
	pass_double_dummy=new double;
	pass_zero_dummy= new double;
}



void Generic_Linear_Algebra::find_eigen_matrix (double **q_mat)
{
	int i,j; 
	
	*infor=0;
	
	for (i=0; i<(*matrix_size); i++)
		for (j=0; j<(*matrix_size); j++)
			pass_matrix[*matrix_size*i+j]=q_mat[j+num_stops_before(j)][i+num_stops_before(i)];
	
	
	*char_1='N';
	*char_2='V';
	
	dgeev_ (char_1, char_2, matrix_size, pass_matrix, leadingd, store_eigen_values, store_c_eigen_v,
			pass_double_dummy,  leadingdl, store_eigen_v, leadingdr, workspace, work_size, infor,
			arg_len, arg_len2);
	
	
	for (i=0; i<(int)(*matrix_size*(*matrix_size)); i++)
		store_inv_eigen_v[i]=store_eigen_v[i];
	
	
	*infor=0;
	
	dgetrf_ (matrix_size, matrix_size, store_inv_eigen_v, leadingd, pivots, infor);
	
	
	*infor=0;
	
	dgetri_(matrix_size, store_inv_eigen_v, leadingd, pivots, workspace, work_size, infor);
	
	
}  //End Generic_Linear_Algebra::find_eigen_matrix





void Generic_Linear_Algebra::find_scaled_matrix(double **end_matrix, double ut, int deriv_level)
{
	
	int i,j;
	
	*infor=0;
	
	
	if (deriv_level == 0)
		for (i=0; i<(*matrix_size); i++)
			for (j=0; j<(*matrix_size); j++)
				pass_matrix_b[*matrix_size*i+j]=
				exp(store_eigen_values[i]*ut)*store_eigen_v[*matrix_size*i+j];
	else if (deriv_level == 1)
		for (i=0; i<(*matrix_size); i++)
			for (j=0; j<(*matrix_size); j++)
				pass_matrix_b[*matrix_size*i+j]=
				store_eigen_values[i]*exp(store_eigen_values[i]*ut)*store_eigen_v[*matrix_size*i+j];
	else
		for (i=0; i<(*matrix_size); i++)
			for (j=0; j<(*matrix_size); j++)
				pass_matrix_b[*matrix_size*i+j]=
				store_eigen_values[i]*store_eigen_values[i]*exp(store_eigen_values[i]*ut)*store_eigen_v[*matrix_size*i+j];
	
	
	for (i=0; i<(int)(*matrix_size*(*matrix_size)); i++)
	{
		pass_matrix_c[i]=store_inv_eigen_v[i];
		pass_matrix[i]=0.0;
	}
	
	*infor=0;
	*char_1='N';
	*char_2='N';
	*pass_double_dummy=1.0;
	*pass_zero_dummy=0.0;
	
	dgemm_ (char_1, char_2, matrix_size, matrix_size, matrix_size, pass_double_dummy, 
			pass_matrix_b,  leadingd, pass_matrix_c, leadingd, pass_zero_dummy, pass_matrix, 
			leadingd, arg_len, arg_len2);	
	
	
	
	
	
	for (i=0; i<*matrix_size; i++)
		for (j=0; j<*matrix_size; j++)
			end_matrix[j+num_stops_before(j)][i+num_stops_before(i)]=pass_matrix[i*(*matrix_size)+j];
	
	
}  //Generic_Linear_Algebra::find_scaled_matrix



int Generic_Linear_Algebra::num_stops_before (int lin_alg_index)
{
	return(0);
}




Linear_Algebra::Linear_Algebra(Genetic_code *ccode) :
Generic_Linear_Algebra(ccode->get_non_stops())
  //Constructor for the Linear_Algebra class--initializes the 
  //variables to be passed to the fortran LAPACK libraries
  //These initializations are for 64X64 matrices
{
  curr_code=ccode;
}




int Linear_Algebra::num_stops_before (int lin_alg_index)
{
  int i, ret_val=0;

 
  while ((ret_val<curr_code->get_num_stops()) && ((lin_alg_index+ret_val) >= curr_code->get_stop_num(ret_val)))
    ret_val++;
  
 
  return(ret_val);
}



Linear_Algebra::~Linear_Algebra()
{
	delete matrix_size;
	delete infor;
	delete work_size;
	delete leadingd;
	delete leadingdl;
	delete leadingdr;
	delete lwork;
	delete pivots;
	delete char_1;
	delete char_2;
	delete[] workspace;
	delete[] pass_matrix;
	delete[] store_eigen_v;
	delete[] store_inv_eigen_v;
	delete[] pass_matrix_c;
	delete[] pass_matrix_b;
	delete[] store_eigen_values;
	delete pass_double_dummy;
	delete pass_zero_dummy;
	
}






