#ifndef __BINS_CPP__						// this is needed since it's a template
#define __BINS_CPP__

#include "binsearch.h"

//**************************************************
template <class Tplt>
Binsearch<Tplt>::Binsearch()
{
}

//**************************************************
template <class Tplt>
Binsearch<Tplt>::~Binsearch()
{

}

//**************************************************
template <class Tplt>
int Binsearch<Tplt>::search(Tplt *list, Tplt &item, int lower, int upper)
{
	int middle;
	int position = -1;
	while(upper >= lower){
		middle = (lower + upper) / 2;
		if(list[middle] == item){
			position = middle;
			break;
		}
		else if(list[middle] < item){
			lower = middle + 1;
		}
		else{
			upper = middle - 1;
		}
	}
//	printf("Binsearch: result = %i, string = %s\n",position, list[middle].ccds_name);
	return(position);
}

//**************************************************

#endif
