#ifndef __CCDS_H__
#define __CCDS_H__

#include <iostream>
#include <string>
#include <stdio.h>
#include "max.h"


class CCDS{										// one line from a header file. Each CCDS represents a single read from one person
public:
	CCDS();
	~CCDS();
	CCDS & operator=(const CCDS &rhs);
	int operator==(const CCDS&);
	int operator!=(const CCDS&);
	int operator>(const CCDS&);
	int operator>=(const CCDS&);
	int operator<(const CCDS&);
	int operator<=(const CCDS&);

	unsigned int ccds_name;							// location
	int start;										// where read starts
	int stop;										// where read stops...not important right now, but we might as well store the info
	int chr;										// chromosome. 1..22, X == 23, Y == 24
//	char *reference;								// reference or expected read
//	char *variant;									// what we read instead of the reference
	int num_reads;									// number of reads
	int num_var;									// number of variant reads
	double proportion;								// of the variant in relation to "num_reads"
	unsigned char *presentin;						// 0|1, only for master record, present in individual or not
	int individuals;								// only for master record, size of presentin
	int iterations;									// how many iterations until it converges, only in masterlist 
	unsigned char *genotype;						// 255 = undetermined, 0 = "ref-ref", 1 = "ref-var", 2 = "var-var", only for master record
	float **geno_prop;								// instead of genotype is lists the probabilities for the 3 genotypes
	float error;									// only keep track in master_list;
	float p;										// prob. of reference, only keep trach in master_list
	float prr;
	float prv;
	float pvv;
	float stat;										// HWE stat or other stuff for debugging and algorithm research purposes
};


#endif
