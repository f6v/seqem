#ifndef __SORT_H__
#define __SORT_H__

//#include "main.h"
#include <stdio.h>
#include <stdlib.h>


//const float missing_covariate = -1024.0;

//struct Record{
//	int famid;                    // family id
//	int ordinal;				  // from 0..(families - 1)
//	float *cov;                   // one or many covariates...depending
//};

const int STACKSIZE = 10000;
const int SMALLSIZE = 20;

struct Stack{
	int a;						// lower border
	int b;						// upper border
};

template <class Tplt>
class Sort{
public:
	Sort();
	~Sort();
	void init(int size, Tplt *list);	// size = number of covariates, key = covariate to be sorted, list = pointer to list of struct Record. This has to exist we will not allocate memory here.
	void swap(int, int);							// switch Record "i" with "j" in "list"
	void push(int i, int j);							// push start/stop on stack
	void pop(int *i, int *j);                           // get stuff of the stack
	void qsort();				// sorts everthing in list*, on "key"
	void split(int first, int last, int *splitpt);
	void selection_sort(int first, int last);

	struct Stack st[STACKSIZE];

private:
	int top;
	Tplt *list;
	Tplt swapper;
	Tplt small;
	Tplt median;											// holds value all others compare against
	int size;									// number of total Records in list
};

#include "sort.cpp"											// yes, this absolutely needs to be there otherwise the Microsoft compiler will give "unresolved external symbol" error

#endif
