#ifndef __SEARCHER_T__
#define __SEARCHER_T__

#include "main.h"
#include "binsearch.h"
#include "ccds.h"
#include "sync_searcher.h"

extern struct Individual *data;
extern struct Sync_searcher *search_syncher;

class Searcher_t{
public:
	Searcher_t();
	~Searcher_t();
	void run(int id);											// function that invokes the sorting
};

#endif
