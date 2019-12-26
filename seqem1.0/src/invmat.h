#ifndef __INVMAT_H__
#define __INVMAT_H__

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

const int maxN = 200;

class Mat{
public:
	Mat();
	~Mat();
	void init(int dim, double **matrix);
	void swap_rows(int, int, int, double**);		// swaps rows
	double** find_inverse();						// finds inverse
	double det();									// finds determinant of matrix
	Mat & operator=(const Mat &rhs);
	Mat & operator*=(const Mat &rhs);				// matrix multiplication
	Mat & operator*=(const int rhs);				// scalar multiplication
	Mat & operator+=(const Mat &rhs);
	Mat & operator-=(const Mat &rhs);
	void print_matrix();
	void print_matrix_to_file(FILE *fptr);			// assumes file is open for writing/appending
	void print_aug();

private:
	int dim;										// dimension of matrix
	double **aug;									// augmented amtrix
	double *row;									// row for interchanging rows in case of leading zero
	double **matrix;								// actual matrix
	double **uppertriag;							// upper triangular form, useful for determinant
};

#endif
