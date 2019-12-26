#ifndef __ANALYZE_H__
#define __ANALYZE_H__

#include "main.h"
#include "invmat.h"
#include "cdflib.h"
#include <string>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// global variables
// extern CCDS* masterlist;
// extern int mastersize;
extern struct Main_list master;
extern struct Individual *data;
extern struct Parameters *para;
extern struct Sync_ana *ana_syncher;

const int parameters = 4;											// namely: p(RR), p(RV), p(VV), error
const double infinity = 1.0e+200;
const double epsi = 1.0e-20;
const float minchisq = 10.0;

struct EM_info{
	int num_threads;
	int id;
	int start;
	int stop;
};

class Analyzer{
public:
	~Analyzer();
	Analyzer();
	void run(struct EM_info *arg);
	void print_result2();
	double absf(const double &val);
	double powexpo(double a, double b);								// returns b*ln(a)
	double diff_of_three(double largest, double a, double b, double c);					// returns how much larger the largest is in relation to 2nd largest number
	double largest_of_three(double a, double b, double c);				// which is the largest? a returns 1, b returns 2, c returns 3
	int largest_of_threei(double a, double b, double c);				// which is the largest? a returns 1, b returns 2, c returns 3

	Mat Cov;																										// complete cov. matrix
	Mat IDM;
	Mat Vm;																											// observed cov. matrix
	Mat Tm;																											// temp matrix
	Mat Jac;
	Mat Inv;
	Mat Dif;

private:
	void em();
	void cov_matrix_no_hwe(int valid_individuals, int iterations, double sum_h, double sum_i);
	void single_em_iteration_no_hwe(int valid_individuals, double input[], double output[]);
	FILE *fptr_seq;													// where we log the sequencing results
	FILE *fptr_matrix;												// where we write various covariance & other matrices
	double input[parameters];										// in case we want to calculate the covariance matirx
	double output[parameters];
	double *p_rr;													// genotype "r,r" = reference, refevernce
	double *p_rv;													// genotype "v,r" = variant, reference
	double *p_vv;													// genotype "v,v," = variant, variant
	double *expected_rr;											// expected number of RR reads
	double *expected_vv;											// expected number of VV reads
	double *expected_r_reads;										// expected r-reads from rr genotype
	double *expected_v_reads;										// expected v-reads from vv genotypes
	int *N;															// number of reads per individual      ie. N
	int *V;															// number of times variant is observed ie. N-X
	int *R;															// number of times reference showed up ie. X
	int *indices;													// index into data array for particular read, -1 == no read in this individual
	double **para_estimates;										// keep estimates for each iteration if HWE == 0
	double **indiv_geno;											// probalisitic genotypes for each individual. ie. indiv_geno[num_individuals][3]
	int max_iter;													// maximum number of iterations for EM-algorithm
	// the following only if calculating Jacobian & Covariance matrices
	int mastersize;
	CCDS *masterlist;
	double **prev;														// previous/initial matrix
	double **current;													// iteration matirx, to be compared to previous
	double **diff;														// difference between current and previous. Once the [i][j] position converged (ie. diff[i][j] < para->em_threshold) we no longer iterate that position
// testing algorithm
	int stab;															// add stabilizer individuals in case genotype(s) missing (like in rare genotypes or HWE==0)
	double naive_thres;													// threshold for naive sequencer IF used (ie. stab > 0)
	double hwecutoff;													// cutoff for chisq test in regards to HWE of SNP data
	int pickinit;														// way to pick initial value (error & p(R)) for EM
	int min_indiv;														// minimum unmber of individuals for setting HWE == 0
	int fixerror;														// 0|1 fixed error
	double init_prr, init_prv, init_pvv;								// initial, user determined, values IF hwe == 0
	// thread related
	int id;																// thread ID
	int num_threads;													// number of threads
	int start;															// first SNP to analyze
	int stop;															// last SNP to analyze
// debug related
	double force_initial_p;												// -1 normally. If in [0, 1] it forces a p-value instead of using num/den as reasonable estimate based on the data at hand
	double initial_error;												// initial error estimate
};

#endif

