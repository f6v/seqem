#ifndef __SORT_CPP__						// this is needed since it's a template
#define __SORT_CPP__

#include "sort.h"

//**************************************************
template <class Tplt>
Sort<Tplt>::Sort()
{
	Sort::list = NULL;
	Sort::size = 0;
	Sort::top = 0;
}

//**************************************************
template <class Tplt>
Sort<Tplt>::~Sort()
{

}

//**************************************************
template <class Tplt>
void Sort<Tplt>::init(int s, Tplt *l)
// s = total number of families
// k = key. i.e. the covariate we sort on. 0..(cov. - 1)
// l = pointer to the list to be sorted
{
	Sort::size = s;					// we absolutely do no error checking. the data coming in is assumed to be correct
	Sort::list = l;
	Sort::top = -1;                       // for stack control. right now stack is empty. 1st element is "0"	
}

//**************************************************
template <class Tplt>
void Sort<Tplt>::split(int first, int last, int *splitpt)
{
	int i,j;
	int pivot = -1;
	int sm,lg,med;
	/// start find pivot element. median of first, last & middle
	med = (first + last) / 2;
	if(list[first] < list[med]){
		sm = first;
		lg = med;
	}
	else{
		sm = med;
		lg = first;
	}
	if(list[last] <= list[sm]){
		pivot = sm;
	}
	else if(list[last] <= list[lg]){
		pivot = last;
	}
	else{
		pivot = lg;
	}
	swap(pivot, first);										// swap pivot element into first place
	median = list[first];								// random array access are expensive, so we store it locally
	i = first + 1;
	j = last;
	while(i < j){
		while((list[j] >= median) && (j > (first + 1))){
			j--;
		}
		while((list[i] < median) && (i < last)){
			i++;
		}
		swap(i, j);											
	}
	swap(i, j);												// must undo last swap since we already had "i<j"
	if(median > list[j]){
		swap(first, j);										// "j" is the position of the last of the small elements, so we swap it with the pivot element which is at "first"
		*splitpt = j;
	}
	else if(median < list[j]){						// due to bad luck median was the smallest element
		*splitpt = first;
	}
	else{
		*splitpt = j;											// must know this to proceed further with this scheme
	}
	// start dbg
//	for(i = first; i <= last; i++){
//		local[i] = list[i].cov[key];
//	}
	// end dbg
}

//**************************************************
template <class Tplt>
void Sort<Tplt>::qsort()
{
	int first = 0;											// first position in array
	int last = size - 1;									// you guessed it! The last position in the array/list
	int splitpt;											// where we break the subarray apart. note that array[splitpt] will be in its correct place
	push(first, last);
	while(top != -1){										// as long as there's stuff on the stack...
		pop(&first, &last);									// get first and last from stack
		for(;;){
			if((last-first) > SMALLSIZE){					// proceed with splitting sub-array
				split(first, last, &splitpt);
				if(first < (splitpt - 1)){
					push(first, splitpt - 1);
				}
				first = splitpt + 1;
			}
			else{											// sort the remainder via selection sort. Ineffiecient but less overhead. works well for small arrays
				selection_sort(first, last);
				break;										// sub-array sorted, get next from stack
			}
		}
	}
}

//**************************************************
template <class Tplt>
void Sort<Tplt>::selection_sort(int first, int last)
{
    int i, j;
	int position;
       
    for(i = first; i < last; i++){							// does nothing if first == last
		small = list[i];
		position = i;
		for(j = i+1; j <= last; j++){
			if(small > list[j]){
				small = list[j];
				position = j;
			}
		}
		if(position != i){
			swap(i, position);
		}
    }

}

//**************************************************
template <class Tplt>
void Sort<Tplt>::swap(int a, int b)
{
	swapper = list[a];
	list[a] = list[b];
	list[b] = swapper;
}

//**************************************************
template <class Tplt>
void Sort<Tplt>::pop(int *a,int *b)                       // pop
{
    *a = st[top].a;								// need those puppies by refernce
    *b = st[top].b;
    top--;											// decrement top of stack
}


//**************************************************
template <class Tplt>
void Sort<Tplt>::push(int a,int b)                        // push
{
    top++;											// increment top of stack BEFORE pushing stuff on there
    st[top].a = a;
    st[top].b = b;
}


//**************************************************

#endif
