#ifndef __SORTER_T__
#define __SORTER_T__

#include "main.h"
#include "sort.h"
#include "ccds.h"
#include "sync_sorter.h"

extern struct Parameters *para;
extern struct Individual *data;
extern struct Sync_sorter *sort_syncher;

class Sorter_t{
public:
	Sorter_t();
	~Sorter_t();
	void run(int id);											// function that invokes the sorting
};

#endif
