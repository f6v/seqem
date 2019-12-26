#ifndef __BINS_H__
#define __BINS_H__


#include <stdio.h>
#include <stdlib.h>
//#include "main.h"


template <class Tplt>
class Binsearch{
public:
	Binsearch();
	~Binsearch();
	int search(Tplt *list, Tplt &item, int min, int max);								// searching for item in list, return -1 if not found, otherwise return position
																				// min & max are first and last position to search, usually 0..(array_size-1)
};

#include "binsearch.cpp"
#endif
