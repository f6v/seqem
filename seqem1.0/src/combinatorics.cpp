#include "combinatorics.h"


//**************************************

int* prime_factors(int N, int &factors)
{
	if(N > max_N){
		fprintf(stderr,"combinatorics.cpp::prime_factors. Argument exceeds max_N = %i. Bye\n", max_N);
		exit(0);
	}
	if(N < 0){
		N *= -1;
	}
	int *factor_list = new int[max_factors];
	memset((void*)factor_list, 0, max_factors*sizeof(int));
	factors = 0;
	int index = 0;
	while(N > 1){
		if(0 == N % prime_archive[index]){
			N /= prime_archive[index];
			factor_list[factors] = prime_archive[index];
			factors++;
		}
		else{
			index++;
		}
	}
	return(factor_list);
}

//**************************************

unsigned long long n_choose_k(int n, int k, int &overflow)
{
	unsigned long long result = 1;							// 0 to 18,446,744,073,709,551,615
	unsigned long long temp;
	overflow = 0;
	if((n == k) || (k == 0)){
		return(result);
	}
	if((k > n) || (n < 1) || (k < 0)){
		fprintf(stderr,"unsigned long long n_choose_k. n = %i, k = %i. Bye\n", n, k);
		exit(0);
	}
	int mastersize = 0;
	int lower_num;
	int upper_den;
	if(k > n/2){											// take advantage of symmetrie
		k = n - k;
	}
	if((n-k) < k){
		upper_den = n-k;
		lower_num = k+1;
	}
	else{
		upper_den = k;
		lower_num = n-k+1;
	}
	int *denlist = new int[upper_den*max_factors];
	int densize = 0;
	int *numlist = new int[(n-lower_num+1)*max_factors];
	int numsize = 0;
	int i, j;
	int fac;
	int *list = NULL;
	for(i = 2; i <= upper_den; i++){							// factorize denominator
		list = NULL;
		list = prime_factors(i, fac);
		for(j = 0; j < fac; j++, densize++){
			denlist[densize] = list[j];
		}
		if(list != NULL){
			delete [] list;
		}
	}
	for(i = n; i >= lower_num; i--){							// factorize numerator
		list = NULL;
		list = prime_factors(i, fac);
		for(j = 0; j < fac; j++, numsize++){
			numlist[numsize] = list[j];
		}
		if(list != NULL){
			delete [] list;
		}		
	}
	Sort <int> sorter;
	sorter.init(densize, denlist);
	sorter.qsort();
	sorter.init(numsize, numlist);
	sorter.qsort();
	// at this point each factor in the denlist MUST be contained in the numlist. a matter of cancelling things out now
	j = 0;
	for(i = 0; i < densize; i++){
		while(denlist[i] != numlist[j]){
			j++;
		}
		numlist[j] = 1;
	}
	for(i = 0; i < numsize; i++){
		temp = result * numlist[i];
		if(temp < result){
			overflow++;
		}
		result = temp;
	}
	delete [] denlist;
	delete [] numlist;
	return(result);
}

//**************************************

long double n_choose_k(int n, int k)
// provides reasonable approximation in case n & k are large
{
	long double result = 1.0;
	if(k > n/2){											// take advantage of symmetrie
		k = n - k;
	}
	int i;
	double temp;
	for(i = 1; i <= k; i++){
		temp =  (double)(n-i+1) / (double)i;
		result *= temp;
	}
	return(result);
}

//**************************************
