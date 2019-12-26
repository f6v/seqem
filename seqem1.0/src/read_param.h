#ifndef __READPARAM_H__
#define __READPARAM_H__

#include <iostream>
#include <fstream>
#include <string.h>
#include <stddef.h>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>

#include "main.h"
#include "max.h"

#if defined (SOLARIS) || defined (LINUX) || (WIN32)
#include <malloc.h>	
#endif

#if defined (OSX)
#include <malloc/malloc.h>	
#endif

#if defined (SOLARIS) || defined (OSX) || defined (OSX) || defined (LINUX)
	#include <sstream>
	#include <backward/strstream>
	using std::ostrstream;
#endif

#if defined (WIN32)
	#include <strstream>
#endif

using namespace std;


const char header_files[max_name] = "Header_files:";
const char hwe[max_name] = "HWE:";
const char debug[max_name] = "Debug:";
const char maxiter[max_name] = "MAX_iterations:";
const char em_thres[max_name] = "EM_threshold:";
const char cov_matrix[max_name] = "Cov_matrix:";
const char outfilename[max_name] = "Outfile:";
const char errorsetting[max_name] = "Initial_error:";
const char fixerror[max_name] = "Fix_error:";
//const char infiletype[max_name] = "Input_format:";
const char min_den[max_name] = "Min_denominator:";
const char threads[max_name] = "Threads:";
const char force_p[max_name] = "Force_pval:";				// debug measure at Eden's request
const char min_depth[max_name] = "Min_read_depth:";			// EM will only consider individuals whose min read depth is at least this number. defaults to 0
const char initval[max_name] = "Pickinitval:";				// 0 or 1. This will enable the follwing options: Stabilizer:, HWEset2false:, Min_individuals4HWEzero: 
const char stabilizer[max_name] = "Stabilizer:";			// 2 naive sequencer parameters: 1. minimum read depth, defaults to 10 2. cutoff, defaults to 0.12 
const char hwechisq[max_name] = "HWEset2false:";			// p-value for chisq test. defaults to 0.005
const char minindiv[max_name] = "Min_individuals4HWEzero:";	// defaults to 15 unless specified otherwise
const char naiveparams[max_name] = "Naive_parameters:";		// cutoff values for using naive sequencer as results, must have Pickinitval: 8
const char prrs[max_name] = "P(rr):";
const char prvs[max_name] = "P(rv):";
class Readparam{

public:
	Readparam();
	~Readparam();
	void set_file_name(char[]);								 // sets nane of private member "filename" name
	int open_infile(ifstream&);								 // this one does NOT asks for filename, use set_file_name 1st
	struct Parameters* master_reader(char filename[]);		// master reader, handles all the parameter reading
	void read_hwe(ifstream&, int&);							// read if HEW in effect
	void read_debug(ifstream&, int&);						// print out a bunch of stuff for each iteration
	void read_iter(ifstream&, int&);
	void read_emthres(ifstream&, double&);
	void read_forcep(ifstream&, double&);
	void read_error(ifstream&, double&);
	void read_cov(ifstream&, int&);
	void read_threads(ifstream&, int&);
	void read_mindepth(ifstream&, int&);
	void read_pickinit(ifstream&, int&);
	void read_fixerror(ifstream&, int&);
	void read_stabilizer(ifstream&, int&, double&);
	void read_naive_params(ifstream&, double&, int&, int&, int&);
	void read_hweswitch(ifstream&, double&);
	void read_init_geno_pval(ifstream&, double&, double&);
	void read_minindiv(ifstream&, int&);
//	void read_inputformat(ifstream&, int&);
	void read_outfile(ifstream&, char[]);
	void read_minden(ifstream&, double&);					// only execute IF COV-matrix is requested, minimum denominator
	char** read_header_files(ifstream&, int&);				// reads the names of the header files

	long find_file_pos(ifstream&, const char[], long); 
	long find_file_pos(fstream&, const char[], long); 	
	void ptr_to_buffer(char *ptr, char buffer[], int size);

private:
	char filename[max_name];
};

#endif
