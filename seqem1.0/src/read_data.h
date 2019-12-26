#ifndef __READ_DATA__
#define __READ_DATA__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include "main.h"
#include "sort.h"
#include "ccds.h"
//#include "md5wrapper.h"

extern struct Individual *data;
extern struct Sync_sorter *sort_syncher;

class Read_data{
public:
	Read_data();
	~Read_data();
	void main_reader(int);				// will open each individual file and read the data
	unsigned char hash_atou8(char *input, int size);
	unsigned short hash_atou16(char *input, int size);
	unsigned int hash_atou32(char *input, int size);

};

//************************************************************
/*
const int STACKSIZE = 10000;
const int SMALLSIZE = 20;

struct Stack{
	int a;						// lower border
	int b;						// upper border
};

template <class T>
class Sort{
public:
	Sort();
	~Sort();
	void init(int size, T *list);	// size = number of covariates, key = covariate to be sorted, list = pointer to list of struct Record. This has to exist we will not allocate memory here.
	void swap(int, int);							// switch Record "i" with "j" in "list"
	void push(int i, int j);							// push start/stop on stack
	void pop(int *i, int *j);                           // get stuff of the stack
	void qsort();				// sorts everthing in list*, on "key"
	void split(int first, int last, int *splitpt);
	void selection_sort(int first, int last);

	struct Stack st[STACKSIZE];

private:
	int top;
	T *list;
	int size;											// number of total Records in list
};
*/

#endif
