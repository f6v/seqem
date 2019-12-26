#ifndef __MAIN_H__
#define __MAIN_H__

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
using namespace std;

#include "read_param.h"
#include "read_data.h"
#include "version.h"
#include "max.h"
#include "ccds.h"
#include "ccds_i.h"
#include "binsearch.h"
#include "combinatorics.h"
#include "invmat.h"
#include "analyze.h"
#include "sort.h"
#include "sync_em.h"
#include "sync_sorter.h"
#include "sync_searcher.h"
#include "sync_ana.h"
#include "sorter_t.h"
#include "searcher_t.h"

//#define PROFILE 1

//#ifdef PROFILE
//#include "profile.h"
//#endif

struct Parameters{									// Input parameter for program
	int num_header_files;
	char **header_file_names;
	int hwe;										// Hardy-Weinberg Eq.
	int debug;										// debug 0|1, prints info on each iteration in EM
	int cov_matrix;									// calculate covaraint amtrix? 0|1 (no | yes)
	int threads;									// how many threads to run
	int min_depth;									// minimum read depth of CCDS to be considered in calculations
	int stab;										// add dummy individuals in case of missing genotypes
	double naive_thres;								// threshold for naive sequencer if used (ie. stab > 0)
	double hwecutoff;								// cutoff for chisq test if data conforms to HWE assumption
	int min_indiv;									// minimum number of individuals for setting HWE=0
	int pickinit;									// initial value scheme to pick error and p(R)
	int fixerror;									// 0 or 1, fixes error estimate to whatever specified under "error"
	double prr;										// for HWE==0 only, initial genotype probability
	double prv;										// for HWE==0 only, initial genotype probability
	double pvv;										// for HWE==0 only, initial genotype probability
// if pickinit == 8, use naive sequencer for results
	double max_naive_error;							// max error to use naive sequencer
	int max_indiv;									// max number of individuals to use naive sequencer
	int max_rdepth;									// max read depth to use naive sequencer
	int and_or;										// to use logical AND or OR to combine above 3
//	int input;										// specifies input format
	int max_iterations;								// max number of iterations in EM algorithm
	double min_denominator;							// minimum allowable denominator for COV-matrix calculations (small denominator = trouble)
	double em_threshold;							// the difference between 2 consecutive iterations has to be less than that for EM to stop
	double forcepval;								// forces initial p-value on EM-algorithm instead of estimate suggested by data at hand
	double error;									// initial error estimate, can be user defined
	char outfile[80];								// file we will print resutls in
};

struct Individual{									// Info for ONE person/individual
	int pedigree_id;								// maybe needed later
	int individual_id;								// ID of individual within pedigree
	int father;										// usual pedigree info standard in all .ped files
	int mother;
	int sex;
	int aff;
	int total_reads;								// total number of reads
	char filename[max_name];						// name of file where info is stored
	CCDSI *ccds_data;								// we will have "total_reads" many of this
};

struct Node{										// for linked list
	CCDSI item;
	struct Node *next;
};

struct Main_list{
	CCDS *list;
	int size;
};

#endif
