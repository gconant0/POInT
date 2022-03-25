//Gavin Conant  10/28/1999
//Class for linear algebra routines needed to calculate
//transition probablities for codon model of evolution
//Uses LAPACK library, either in subdirectory /lin_alg
//or machine local
#include <iostream>


#ifndef ___LIN_ALG_H___
#define ___LIN_ALG_H___

#include "liblapack_orig.h"
#include "gen_code.h"

class Generic_Linear_Algebra {
	public:
	//The variables must be public so the external code
	//can access them
	long int *matrix_size, *infor, *work_size, 
    *leadingd, *leadingdl, *leadingdr, *lwork, *pivots;
	long int *arg_len, *arg_len2;
	char *char_1, *char_2;
	double *workspace, *pass_matrix, *pass_matrix_b, *pass_matrix_c,
    *store_eigen_v, *store_c_eigen_v, 
    *store_inv_eigen_v, *store_eigen_values, *pass_double_dummy,
    *pass_zero_dummy;
	
	//Functions
    Generic_Linear_Algebra()   {std::cerr<<"Call to default Linear_Algebra constructor\n";};
	Generic_Linear_Algebra(int matrix_s);
	void find_eigen_matrix (double **q_mat);
	virtual int num_stops_before (int lin_alg_index);
	void find_scaled_matrix(double **end_matrix, double ut, int deriv_level);
    virtual ~Generic_Linear_Algebra() {};
	
	
};

class Linear_Algebra : public Generic_Linear_Algebra
{
public:
   //Functions
    Linear_Algebra()   {std::cerr<<"Call to default Linear_Algebra constructor\n";};
  Linear_Algebra(Genetic_code *curr_code);
  int num_stops_before (int lin_alg_index);
    ~Linear_Algebra();

private:
 Genetic_code *curr_code;
  
};


#endif













