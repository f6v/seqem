#include "invmat.h"

//*********************************

Mat::Mat()
{
	Mat::aug = NULL;
	Mat::row = NULL;
	Mat::matrix = NULL;
	Mat::uppertriag = NULL;
	Mat::dim = 0;
}

//*********************************

Mat::~Mat()
{
	int i;
	if(NULL != row){
		delete [] row;
	}
	if(NULL != aug){
		for(i = 0; i < dim; i++){
			delete [] aug[i];
		}
		delete [] aug;
	}
	if(NULL != matrix){
		for(i = 0; i < dim; i++){
			delete [] matrix[i];
		}
		delete [] matrix;
	}
	if(NULL != uppertriag){
		for(i = 0; i < dim; i++){
			delete [] uppertriag[i];
		}
		delete [] uppertriag;
	}
}

//*********************************

void Mat::init(int n, double **m)
{
	int i;
	if(Mat::dim == n){
		for(i = 0; i < dim; i++){
			memcpy((void*)matrix[i], (void*)m[i], dim*sizeof(double));
		}
	}
	else if(Mat::dim > 0){
		for(i = 0; i < dim; i++){
			delete [] matrix[i];
		}
		delete [] matrix;
		dim = 0;
		matrix = NULL;
	}
	if(Mat::dim == 0){
		dim = n;
		matrix = new double*[dim];
		for(i = 0; i < dim; i++){
			matrix[i] = new double[dim];
			memcpy((void*)matrix[i], (void*)m[i], dim*sizeof(double));
		}
	}
}

//*********************************

void Mat::swap_rows(int a, int b, int size, double **dim2array)
{
	memcpy((void*)row, (void*)dim2array[a], size*sizeof(double));
	memcpy((void*)dim2array[a], (void*)dim2array[b], size*sizeof(double));
	memcpy((void*)dim2array[b], (void*)row, size*sizeof(double));
}

//*********************************

double** Mat::find_inverse()
{
	double **inverse = NULL;
	int i,j,k,l;
	int n = dim;
	int non_zero_pivot;
	double pivot;
	int inverse_exists = 1;
	aug = new double*[n];
	for(i = 0; i < n; i++){
		aug[i] = new double[2*n];
		memset((void*)aug[i], 0, 2*n*sizeof(double));
		memcpy((void*)aug[i], (void*)matrix[i], n*sizeof(double));
		j = n+i;
		aug[i][j] = 1.0;
	}
	row = new double[2*n];
	for(j = 0; j < dim; j++){												// columns
		non_zero_pivot = 0;
		i = j;
		while(i < dim){														// find non-zero pivot element
			if(0 != aug[i][j]){
				non_zero_pivot = 1;
				break;
			}
			else{
				i++;
			}
		}
		if(0 == non_zero_pivot){
			inverse_exists = 0;
			break;															// we will return NULL pointer, matrix has no inverse
		}
		if(i != j){
			swap_rows(i, j, 2*dim, aug);									// get non-zero pivot element into [j][j] position
		}
		pivot = aug[j][j];
		aug[j][j] = 1;
		i = j;
		for(k = i+1; k < 2*dim; k++){										// pivot row divided by pivot element
			aug[i][k] /= pivot;		
		} // for(k = i+1; k < dim; k++)
		k = 0;																// now perform addition/subtr so that the j_th column has all zeros except the 1 in the j_th row
		while(k < dim){
			if(i != k){
				pivot = aug[k][j];
				for(l = j; l < 2*dim; l++){
					aug[k][l] -= aug[i][l] * pivot;
				} // for(l = j; l < 2*dim; l++)
			} // if(i != k)
			k++;
		} // while(k < dim)
	} // for(j = 0; j < dim; j++)
	if(inverse_exists){
		inverse = new double*[dim];
		for(i = 0; i < dim; i++){
			inverse[i] = new double[dim];
			for(j = 0; j < dim; j++){
				inverse[i][j] = aug[i][dim+j];
			}
		}
	}
	for(i = 0; i < dim; i++){
		delete [] aug[i];
	}
	delete [] aug;
	aug = NULL;
	return(inverse);
}

//*********************************

double Mat::det()
{
	double pivot, factor;
	int i,j,k;
	int non_zero_pivot;
	double result;
	uppertriag = new double*[dim];
	for(i = 0; i < dim; i++){
		uppertriag[i] = new double[dim];
		memcpy((void*)uppertriag[i], (void*)matrix[i], dim*sizeof(double));
	}
	// now try to convert this matrix into upper triangular form
	for(j = 0; j < dim; j++){
		non_zero_pivot = 0;
		i = j;
		while(i < dim){														// find non-zero pivot element
			if(0 != uppertriag[i][j]){
				non_zero_pivot = 1;
				break;
			}
			else{
				i++;
			}
		}
		if(0 == non_zero_pivot){
			result = 0;														// singular matirx, det(A) = 0
			break;															// we will return NULL pointer, matrix has no inverse
		}	
		if(i != j){
			swap_rows(i, j, dim, uppertriag);								// get non-zero pivot element into [j][j] position
		}
		pivot = uppertriag[j][j];
		for(i = j+1; i < dim; i++){											// for each row below
			if(0 != uppertriag[i][j]){										// if the leading entry != 0
				factor = uppertriag[i][j] / pivot;								
				for(k = j; k < dim; k++){
					uppertriag[i][k] -= factor * uppertriag[j][k];			// subtract "factor * pivot_row" from i_th row
				} // for(k = j; k < dim; k++)
			} // if(0 != temp[i][j])
		} // for(i = j+1; i < dim; i++)
	} // for(j = 0; j < dim; j++)
	result = 1;
	for(i = 0; i < dim; i++){												// determinant is product of main diagonal
		result *= uppertriag[i][i];
	}
	return(result);
}

//*********************************

void Mat::print_matrix()
{
	int i,j;
	for(i = 0; i < dim; i++){
		for(j = 0; j < dim; j++){
			printf("%8.4f ",matrix[i][j]);
		}
		printf("\n");
	}
}

//*********************************

void Mat::print_matrix_to_file(FILE *fptr)
{
	int i,j;
	for(i = 0; i < dim; i++){
		for(j = 0; j < dim; j++){
			fprintf(fptr,"%10.6f ",matrix[i][j]);
		}
		fprintf(fptr,"\n");
	}
}

//*********************************

void Mat::print_aug()
{
	int i,j;
	for(i = 0; i < dim; i++){
		for(j = 0; j < 2*dim; j++){
			printf("%8.4f ", aug[i][j]);
		}
		printf("\n");
	}
}

//*********************************

//*********************************

Mat & Mat::operator=(const Mat &rhs)
{
	int i;
	if(&rhs != this){
		if(rhs.dim == this->dim){													// check if they're the same size
			for(i = 0; i < this->dim; i++){
				memcpy((void*)this->matrix[i], (void*)rhs.matrix[i], dim*sizeof(double));
			} // for(i = 0; i < this->dim; i++)
		} // if(rhs.dim == this->dim)
		else if(this->dim > 0){														// need to first delete old
			for(i = 0; i < this->dim; i++){
				delete [] matrix[i];
			}
			delete [] matrix;
			this->dim = 0;
			this->matrix = NULL;
		}
		if(this->dim == 0){
			this->dim = rhs.dim;
			this->matrix = new double*[dim];
			for(i = 0; i < this->dim; i++){
				this->matrix[i] = new double [dim];
				memcpy((void*)this->matrix[i], rhs.matrix[i], dim*sizeof(double));
			}
		}
	}
	else{
		fprintf(stderr,"Mat::operator=(const Mat &rhs) Detected self assignment. Bailing out.\n");
		exit(0);
	}
	return(*this);
}

//*********************************

Mat & Mat::operator*=(const Mat &rhs)
{
	int i,j,k;
	double **res;
	if(rhs.dim == this->dim){
		res = new double*[this->dim];
		for(i = 0; i < this->dim; i++){
			res[i] = new double[this->dim];
		}
		for(i = 0; i < dim; i++){
			for(j = 0; j < dim; j++){
				res[i][j] = 0;
				for(k = 0; k < dim; k++){
					res[i][j] += this->matrix[i][k] * rhs.matrix[k][j]; 
				}
			}
		}
		for(i = 0; i < this->dim; i++){
			memcpy((void*)this->matrix[i], (void*)res[i], this->dim*sizeof(double));
			delete [] res[i];
		}
		delete [] res;
	}
	else{
		fprintf(stderr,"Mat::operator*=(const Mat &rhs) this->dim = %i\n", this->dim);
		exit(0);
	}
	return(*this);
}

//*********************************

Mat & Mat::operator*=(const int rhs)
{
	int i,j;
	if(this->dim > 0){
		for(i = 0; i < dim; i++){
			for(j = 0; j < dim; j++){
				this->matrix[i][j] *= rhs;
			}
		}
	}
	else{
		fprintf(stderr,"Mat::operator*=(const int rhs) this->dim = %i\n", this->dim);
		exit(0);
	}
	return(*this);
}

//*********************************

Mat & Mat::operator+=(const Mat &rhs)
{
	int i,j;
	if(rhs.dim == this->dim){
		for(i = 0; i < this->dim; i++){
			for(j = 0; j < this->dim; j++){
				this->matrix[i][j] += rhs.matrix[i][j];
			}
		}
	}
	else{
		fprintf(stderr,"Mat::operator+= Dimension mismatch. this->dim = %i, rhd.dim = %i\n", this->dim, rhs.dim);
		exit(0);
	}
	return(*this);
}

//*********************************

Mat & Mat::operator-=(const Mat &rhs)
{
	int i,j;
	if(rhs.dim == this->dim){
		for(i = 0; i < this->dim; i++){
			for(j = 0; j < this->dim; j++){
				this->matrix[i][j] -= rhs.matrix[i][j];
			}
		}	
	}
	else{
		fprintf(stderr,"Mat::operator+= Dimension mismatch. this->dim = %i, rhd.dim = %i\n", this->dim, rhs.dim);
		exit(0);
	}
	return(*this);
}

//*********************************
