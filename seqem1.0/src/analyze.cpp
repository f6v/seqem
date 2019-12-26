#include "analyze.h"

//***********************************

Analyzer::~Analyzer()
{
	int i;
	delete [] N;
	delete [] V;
	delete [] R;
	delete [] p_rr;
	delete [] p_rv;
	delete [] p_vv;
	delete [] expected_rr;
	delete [] expected_vv;
	delete [] expected_r_reads;
	delete [] expected_v_reads;
	delete [] indices;
	for(i = 0; i < para->num_header_files; i++){
		delete [] indiv_geno[i];
	} // for(i = 0; i < num_individuals; i++)
	delete [] indiv_geno;

	if(para->cov_matrix){
		for(i = 0; i < parameters; i++){
			delete [] prev[i];
			delete [] current[i];
			delete [] diff[i];
		}
		delete [] prev;
		delete [] current;
		delete [] diff;
		for(i = 0; i < max_iter; i++){
			delete [] para_estimates[i];
		}
		delete [] para_estimates;
	}
	fclose(fptr_seq);
	fclose(fptr_matrix);
}

//***********************************

Analyzer::Analyzer()
{
	double *p_rr = NULL;
	double *p_rv = NULL;
	double *p_vv = NULL;
	double *expected_rr = NULL;
	double *expected_vv = NULL;
	double *expected_r_reads = NULL;
	double *expected_v_reads = NULL;
	int *N = NULL;
	int *V = NULL;
	int *R = NULL;
	double **para_estimates = NULL;											// keep estimates for each iteration if HWE == 0
	int *indices = NULL;
	double **prev = NULL;
	double **current = NULL;
	double **diff = NULL;
	double **indiv_geno = NULL;
}

//***********************************

void Analyzer::run(struct EM_info *arg)
{
	int i;
	Analyzer::start = arg->start;
	Analyzer::stop = arg->stop;
	Analyzer::id = arg->id;
printf("I am EM-thread %i\n", id);
	Analyzer::num_threads = arg->num_threads;
	Analyzer::stab = para->stab;												// add stabilizer if genotype(s) are missing
	Analyzer::naive_thres = para->naive_thres;
	Analyzer::hwecutoff = para->hwecutoff;
	Analyzer::pickinit = para->pickinit;
	Analyzer::min_indiv = para->min_indiv;
	Analyzer::fixerror = para->fixerror;
	Analyzer::init_prr = para->prr;
	Analyzer::init_prv = para->prv;
	Analyzer::init_pvv = para->pvv;
	Analyzer::p_rr = new double[para->num_header_files+3];						// genotype "r,r" = reference, refevernce
	Analyzer::p_rv = new double[para->num_header_files+3];						// genotype "v,r" = variant, reference
	Analyzer::p_vv = new double[para->num_header_files+3];						// genotype "v,v," = variant, variant
	Analyzer::expected_rr = new double[para->num_header_files+3];				// expected number of RR reads
	Analyzer::expected_vv = new double[para->num_header_files+3];				// expected number of VV reads
	Analyzer::expected_r_reads = new double[para->num_header_files+3];			// expected r-reads from rr genotype
	Analyzer::expected_v_reads = new double[para->num_header_files+3];			// expected v-reads from vv genotypes
	Analyzer::N = new int[para->num_header_files+3];								// number of reads per individual      ie. N
	Analyzer::V = new int[para->num_header_files+3];								// number of times variant is observed ie. N-X
	Analyzer::R = new int[para->num_header_files+3];								// number of times reference showed up ie. X
	Analyzer::para_estimates = NULL;											// keep estimates for each iteration if HWE == 0
	Analyzer::indices = new int[para->num_header_files+3];							// index into data array for particular read, -1 == no read in this individual
	Analyzer::max_iter = para->max_iterations;									// max iterations for EM-algorithm
	Analyzer::force_initial_p = para->forcepval;								// forces initial p-value instead of the one derived from the data at hand
	if(para->cov_matrix){													// allocate memory to para_estimates, so we can remember those values instead of recalculating them
		para_estimates = new double*[max_iter];
		for(i = 0; i < max_iter; i++){
			para_estimates[i] = new double[4];								// p(RR), p(RV), p(VV), error-estimate
		}
//		fptr1 = fopen("cov_matrix.txt", "w");
		prev = new double*[parameters];
		current = new double*[parameters];
		diff = new double*[parameters];
		for(i = 0; i < parameters; i++){
			prev[i] = new double[parameters];
			current[i] = new double[parameters];
			diff[i] = new double[parameters];
		}
	}
	indiv_geno = new double*[para->num_header_files+3];						// The accumulator from equation (3)
	for(i = 0; i < (para->num_header_files+3); i++){
		indiv_geno[i] = new double[3];
	} // for(i = 0; i < num_individuals; i++)
	Analyzer::fptr_seq = NULL;
	Analyzer::fptr_matrix = NULL;
//	fptr_seq = fopen(para->outfile, "w");
	if(para->cov_matrix){
		fptr_matrix = fopen("matrices.txt", "w");
	}
	Analyzer::mastersize = master.size;
	Analyzer::masterlist = master.list;
	// done allocating memory
	em();
}
//***********************************

void Analyzer::em()
{
	// Global parameters:
	/*
CCDS* masterlist, 
int mastersize, 
Individual *data, 
Parameters *para
	*/
	// hwe = hardy-weinberg assumption, 0 = NO, 1 = YES
/*
General Parameterization	Assuming HWE
Parameters:
PRR = Prior probability of RR homozygote
PRV = Prior probability of RV homozygote 
PVV = Prior probability of VV homozygote	P = Prior probability of R allele
e = Probability that the wrong allele is read
Data:
There are m individuals indexed by i with data points (Xi, Ni), where:
Xi = Number of R reads for individual i
Ni = Total reads (both R and V) for individual i
EM Algorithm:
I. Initialize Parameter Estimates for Iteration 0:
1) Initialize genotype priors from the data:
 
where  
1) Initialize allele prior from the data:
 

2) Set   or another reasonable value.

II. E Step:
1) For each individual i in iteration r calculate the following joint probabilities for Xi and the genotype:
(A) = 
(B) = 
(C) = 
1) For each individual i in iteration r calculate the following joint probabilities for Xi and the genotype:
(A) = 
(B) = 
(C) = 

2) For each individual i, calculate the following:

(D) = (A) + (B) + (C) 	(Marginal probability of Xi)

(E) = (A)/(D) 		(Conditional probability of RR genotype given Xi)

(F) = (B)/(D) 	(Conditional probability of RV genotype given Xi)

(G) = (C)/(D)	(Conditional probability of VV genotype given Xi)

(H) = Ni*(E) 	(Expected number of total reads from RR genotype in i)	

(I) = Ni*(G)	(Expected number of total reads from VV genotype in i)

(J) = Xi*(E)	(Expected number of R reads from RR genotype in i)

(K) = Xi*(G)	(Expected number of R reads from VV genotype in i)

3) For quantities (E) through (K), take the sum over all individuals.  I will denote the resulting sums by SUM(E) through SUM(K).  These sums are the full sample expectations that are used in the M step.
III. M Step:
1) Update genotype prior estimates as follows:
 
1) Update the allele prior estimate as follows:
 

2) Update the error rate estimate as follows:
 

IV. If any parameter estimate changed by more than some tolerance level (e.g., 1e-8) from the previous iteration, increment the counter r and return to II with the updated parameter estimates.  Otherwise, stop.
 
*/
	int i,j,k;
	int temp;
	int length;
	int iter;														// iteration
	long long unsigned bin_coeff;									// binomial coefficient
	long double bc;													// binomila coeff. as double
	int overflow = 0;												// checks for overflow when calculating bin_coeff
	double nd, kd, td, dd;											// temp variables for calculating n choose k
	int ni, ki;														// temp variables for calculating n choose k
	double denom_sum;												// denominator sum in eq. (5)
	double num_sum;													// numerator sum
	double p, prev_p, p_obs;
	double prr,prv,pvv;												// genotype frequency estimates
	double prev_prr, prev_prv, prev_pvv;							// previous estimates
	double error, prev_error;
	double pinit, errorinit;
	double largest, diff;
	double sum_e, sum_g, sum_f, sum_k, sum_h, sum_j, sum_i;
	char *tmp;														// temp storage

	int maxN;														// max number of reads
	double *one_err_R;
	double *one_err_R_N;
	double *err_R;
	double *err_R_N;
	double *half;
	unsigned char *tally;

	double *one_err_Rs;
	double *one_err_R_Ns;
	double *err_Rs;
	double *err_R_Ns;
	double *halfs;
	unsigned char *tallys;

	double td1,td2,td3,td4,td5;
	double sumrr, sumrv, sumvv, exprr, exprv, expvv;
	double q;

	int continue_EM;
	double valid_individuals;											// number of individuals that actually have one or more reads
	Binsearch <CCDSI> searcher;
	CCDSI item;
//#ifdef PROFILE
//	PROFILE_START("EM_net");
//#endif
//	for(i = start; i <= stop; i++){									// start main loop
	i = start;
	int done = 0;
	int troubledexpo;
	int prelim_rr = 0;
	int prelim_rv = 0;
	int prelim_vv = 0;
	int observed_genotypes;
	double naive_geno_sum;
	double wrong_v = 0;
	double wrong_r = 0;
	double total_for_error_cal = 0;
	double expchi_rr, expchi_vv, expchi_rv;
	int hwe;
	double logsum;
	double upper_thres = 1 - naive_thres;
//	double estimated_maf;
//	double estimated_error;
	int number_of_individuals;
	// for cdfchi
	int code = 1;
	int status = 0;
	double chi_p, chi_q, df, bound, global;
	float convert;

	while(done == 0){
		valid_individuals = 0;
		num_sum = 0;
		denom_sum = 0;
//#ifdef PROFILE
//	PROFILE_START("EM_search");
//#endif
		one_err_R = NULL;
		one_err_R_N = NULL;
		err_R = NULL;
		err_R_N = NULL;
		half = NULL;
		tally = NULL;
		maxN = -1;
		// stabilizer related
		prelim_rr = 0;
		prelim_rv = 0;
		prelim_vv = 0;
		naive_geno_sum = 0;
		wrong_v = 0;
		wrong_r = 0;
		total_for_error_cal = 0;
		
		Analyzer::initial_error = para->error;										// initial error estimate for EM

		number_of_individuals = para->num_header_files;
		max_iter = para->max_iterations;
		hwe = para->hwe;
		masterlist[i].stat = -1.0;
		// emd stab
		for(j = 0; j < para->num_header_files; j++){
			item.location = (unsigned int)masterlist[i].ccds_name;
			temp = searcher.search(data[j].ccds_data, item, 0, data[j].total_reads-1);
			indices[j] = temp;
			if(indices[j] >= 0){
				N[j] = data[j].ccds_data[indices[j]].depth;
				if(N[j] <= 0){
					N[j] = -1;
					V[j] = -1;
					R[j] = -1;		
					indices[j] = -1;
				}
				else{
					V[j] = data[j].ccds_data[indices[j]].variants;
					R[j] = N[j] - V[j];
// naive sequencer related
					if(pickinit && (N[j] >= stab)){								// the minimum number of reads to make any kind of genotype call >= 1
						if((naive_thres*N[j]) > V[j]){
							prelim_rr++;
							wrong_v += V[j];
							total_for_error_cal += N[j];
						}
						else if(V[j] > (upper_thres*N[j])){
							prelim_vv++;
							wrong_r += R[j];
							total_for_error_cal += N[j];
						}
						else{
							prelim_rv++;
						}
					}
// end stab
//				data[j].ccds_data[indices[j]].num_var = V[j];
					valid_individuals += 1;
					denom_sum += N[j];
					num_sum += R[j];
					if(N[j] > maxN){
						maxN = N[j];
					}
				} // if(N[j] <= 0){
			} // if(indices[j] >= 0)
			else{
				N[j] = -1;
				V[j] = -1;
				R[j] = -1;
			}
		} // for(j = 0; j < para->num_header_files; j++)
//#ifdef PROFILE
//    PROFILE_STOP("EM_search");
//#endif
		iter = 0;

		prev_p = 1;
		if(denom_sum > 0){
			if((force_initial_p >= 0) && (force_initial_p <= 1.0)){
				p = force_initial_p;
			}
			else{
				p = num_sum / denom_sum;
			}
			if((hwe == 0) && (para->prr >= 0)){
				prr = para->prr;
				prv = para->prv;
				pvv = 1 - prr - prv;
			}
			else{
				prr = p*p;
				prv = 2*(1-p)*p;
				pvv = (1-p)*(1-p);
			}
		}
		else{
			fprintf(stderr,"main::analyze3. Denominator for p-estimation is zero.item_name= %u\n", item.location);
			p = -1.0;
		}
		error = initial_error;
		if(fixerror == 0){
			prev_error = 0.99;
		}
		else{
			prev_error = error;
		}
//		if(hwecutoff > 0){
		// initial values according to the naive sequencer
		if(pickinit){
			if(total_for_error_cal > 0){
				errorinit = (wrong_v + wrong_r) / total_for_error_cal;			// estimated error
			}
			else{
				errorinit = initial_error;										// in case we have zero denominator
			}
			naive_geno_sum = (double)(prelim_rr + prelim_rv + prelim_vv);
			if(naive_geno_sum > 0){
				if(denom_sum > 0){
					p_obs = num_sum / denom_sum;
					if(errorinit < 0.3){
						pinit = (p_obs - errorinit) / (1 - 2*errorinit);
					}
					else{
						pinit = p_obs;
					}
					if(p > 1){
						pinit = 0.999;
					}
				}
				else{
					pinit = (prelim_rv + 2*prelim_rr) / (2*naive_geno_sum);	// estimated "1 - MAF"
				}
			}
		
		// check for serious violation of HWE based on the naive numbers
			if((pickinit & 2) && (errorinit > 0)){
				error = errorinit;
			}
			if((pickinit & 4) && (pinit >0) && (pinit < 1)){
				p = pinit;
			}
			if((pickinit == 8) && (valid_individuals == (prelim_rr + prelim_rv + prelim_vv))){	
				td = denom_sum / naive_geno_sum;
				if(para->and_or == 1){ // AND
					if((para->max_naive_error > errorinit) && (para->max_indiv >= valid_individuals) && (para->max_rdepth > td)){
						hwe = 1;
						max_iter = 1;
						error = naive_thres;
						p = pinit;
					}
				}
				else if(para->and_or == 0){ // OR
					if((para->max_naive_error > errorinit) || (para->max_indiv >= valid_individuals) || (para->max_rdepth > td)){
						hwe = 1;
						max_iter = 1;

					}
				}
			}
			observed_genotypes = 0;
			if(prelim_rr > 0) observed_genotypes++;
			if(prelim_rv > 0) observed_genotypes++;
			if(prelim_vv > 0) observed_genotypes++;
			if((hwecutoff > 0) && (naive_geno_sum > (double)(min_indiv - 0.1)) && (pinit > 0) && (pinit < 1) && (pickinit & 1)){
				expchi_rr = naive_geno_sum*pinit*pinit;
				expchi_rv = naive_geno_sum*2*pinit*(1-pinit);
				expchi_vv = naive_geno_sum*(1-pinit)*(1-pinit);
				bound = 0;
				df = 1.0;
				global = (expchi_rr - prelim_rr) * (expchi_rr - prelim_rr) / expchi_rr;
				global += (expchi_rv - prelim_rv) * (expchi_rv - prelim_rv) / expchi_rv;
				global += (expchi_vv - prelim_vv) * (expchi_vv - prelim_vv) / expchi_vv;
				cdfchi(&code, &chi_p, &chi_q, &global, &df,  &status, &bound);
				masterlist[i].stat = (float)chi_q;
				if((status == 0) && (chi_q < hwecutoff)){
					hwe = 0;
				}
			}
			if((hwe == 0) && (observed_genotypes == 2)){
				fixerror = 1;
			}
			else if((hwe == 0) && (observed_genotypes == 1)){
				hwe = 1;
				fixerror = 1;
			}
			if((hwe == 0) && (naive_geno_sum > (double)(min_indiv - 0.1))){
				prr = (double)prelim_rr / naive_geno_sum;
				prv = (double)prelim_rv / naive_geno_sum;
				pvv = (double)prelim_vv / naive_geno_sum;	
			}
			else{
				prr = p*p;
				prv = 2*(1-p)*p;
				pvv = (1-p)*(1-p);			
			}

		} // if(pickinit)
		for(j = 0; j < para->num_header_files; j++){
			memset((void*)indiv_geno[j], 0, 3*sizeof(double));		// the 3 possible genotypes
		}
		maxN++;
		if(maxN < 101){
			maxN = 101;
		}
		one_err_R = new double[maxN];
		one_err_R_N = new double[maxN];
		err_R = new double[maxN];
		err_R_N = new double[maxN];
		half = new double[maxN];
		tally = new unsigned char[maxN];

		memset((void*)one_err_R, 0, maxN*sizeof(double));
		memset((void*)one_err_R_N, 0, maxN*sizeof(double));
		memset((void*)err_R, 0, maxN*sizeof(double));
		memset((void*)err_R_N, 0, maxN*sizeof(double));
		memset((void*)half, 0, maxN*sizeof(double));

		one_err_Rs = new double[maxN];
		one_err_R_Ns = new double[maxN];
		err_Rs = new double[maxN];
		err_R_Ns = new double[maxN];
		halfs = new double[maxN];
		tallys = new unsigned char[maxN];

		memset((void*)one_err_Rs, 0, maxN*sizeof(double));
		memset((void*)one_err_R_Ns, 0, maxN*sizeof(double));
		memset((void*)err_Rs, 0, maxN*sizeof(double));
		memset((void*)err_R_Ns, 0, maxN*sizeof(double));
		memset((void*)halfs, 0, maxN*sizeof(double));

		memset((void*)p_rr, 0, number_of_individuals*sizeof(double));
		memset((void*)p_rv, 0, number_of_individuals*sizeof(double));
		memset((void*)p_vv, 0, number_of_individuals*sizeof(double));
		memset((void*)expected_rr, 0, number_of_individuals*sizeof(double));
		memset((void*)expected_vv, 0, number_of_individuals*sizeof(double));
		memset((void*)expected_r_reads, 0, number_of_individuals*sizeof(double));
		memset((void*)expected_v_reads, 0, number_of_individuals*sizeof(double));


		if(para->debug & 2){
			printf("### START SNP = %i ###\n", masterlist[i].ccds_name);
		}
		if(p >= 0){
			continue_EM = 1;
		}
		else{
			continue_EM = 0;
		}
		if(para->debug & 2){
			printf("initial error = %7.5f initial p(r) = %7.5f HWE= %i\n",error,p,hwe);
			if(hwe == 0){
				printf("Initial genotype probabilities: P(rr)= %7.5f P(rv)= %7.5f P(vv)= %7.5f\n", prr,prv,pvv);
			}
		}
		while((continue_EM) && (iter < max_iter)){					// main EM iteration(s)
			sum_e = sum_g = sum_f = sum_k = sum_h = sum_j = sum_i = 0;
			memset((void*)tally, 0, maxN*sizeof(unsigned char));
			memset((void*)tallys, 0, maxN*sizeof(unsigned char));
			if(para->debug & 2) logsum = 0;
			for(j = 0; j < number_of_individuals; j++){
				if(indices[j] >= 0){								// only those individuals have valid data

	// the following section seems a bit mysterious. tally[] simply keeps track of various coefficients
	// so that they don't have to be recalculated for every individual. The pow() function is very costly
	// and should only be used sparingly. We store coeffiecients in one_err_R_N[], one_err_R[], err_R[], half[]
					if(N[j] > 60){									// large read depth need to use expensive log() function
						if(tally[R[j]] & 1){						// (1-err)^X_i
							td1 = one_err_R[R[j]];
						}
						else{
							td1 = powexpo(1-error,R[j]);
							one_err_R[R[j]] = td1;
							tally[R[j]] |= 1;
						}

						if(tally[N[j]-R[j]] & 2){					// (1-err)^(N_i - X_i)
							td2 = one_err_R_N[N[j]-R[j]];
						}
						else{
							td2 = powexpo(1-error, N[j]-R[j]);
							one_err_R_N[N[j]-R[j]] = td2;
							tally[N[j]-R[j]] |= 2;
						}

						if(tally[R[j]] & 4){						// err^X_i
							td3 = err_R[R[j]];
						}
						else{
							// td3 = pow(error, R[j]) * pvv;
							if(pvv < epsi){
								td3 = -infinity;
							}
							else{
								td3 = powexpo(error, R[j]) + log(pvv);
								err_R[R[j]] = td3;
								tally[R[j]] |= 4;
							}
						}

						if(tally[N[j]-R[j]] & 8){					// err^(N_i - X_i)
							td4 = err_R_N[N[j]-R[j]];
						}
						else{
							// td4 =  pow(error, N[j]-R[j]) * prr;
							if(prr < epsi){
								td4 = -infinity;
							}
							else{
								td4 =  powexpo(error, N[j]-R[j]) + log(prr);
								err_R_N[N[j]-R[j]] = td4;
								tally[N[j]-R[j]] |= 8;
							}
						}
// here the real calculations start
						if(hwe){
							if(tally[N[j]] & 16){					// 0.5^(N_i - 1)
								td5 = half[N[j]];
							}
							else{
								// td5 = pow(0.5, N[j]-1) * p * (1-p);
								if((p > 0) && (p < 1.0)){
									td5 = powexpo(0.5, N[j]-1) + log(p) + log((1-p));
									half[N[j]] = td5;
									tally[N[j]] |= 16;
								}
								else{
									td5 = -infinity;
								}
							}
							//p_rr[j] = bc * td1 * td4;
							//p_rv[j] = bc * td5;
							//p_vv[j] = bc * td2 * td3;

						}
						else{ // if(hwe)
							if(tally[N[j]] & 16){					// 0.5^N_i
								td5 = half[N[j]];
							}
							else{
								//td5 = pow(0.5, N[j]) * prv;
								if(prv > epsi){
									td5 = powexpo(0.5, N[j]) + log(prv);
									half[N[j]] = td5;
									tally[N[j]] |= 16;
								}
								else{
									td5 = -infinity;
								}
							}
						//p_rr[j] = bc * td1 * td4;
						//p_rv[j] = bc * td5;	
						//p_vv[j] = bc * td2 * td3;

						} // if(hwe){} else {}
						kd = td1 + td4;
						nd = td2 + td3;

						largest = largest_of_three(kd, nd, td5);
						kd -= largest;
						nd -= largest;
						td5 -= largest;
						p_rr[j] = 0;
						p_rv[j] = 0;
						p_vv[j] = 0;
						if(100.0 > absf(kd)){
							p_rr[j] = exp(kd);
						}
						if(100.0 > absf(nd)){
							p_vv[j] = exp(nd);
						}						
						if(100.0 > absf(td5)){
							p_rv[j] = exp(td5);
						}

	
					}
					else{ // if(N[j] > 60) ...so we're dealing with a "managable" read depth 
						if(tallys[R[j]] & 1){
							td1 = one_err_Rs[R[j]];
						}
						else{
							td1 = pow(1-error,R[j]);
							one_err_Rs[R[j]] = td1;
							tallys[R[j]] |= 1;
						}

						if(tallys[N[j]-R[j]] & 2){
							td2 = one_err_R_Ns[N[j]-R[j]];
						}
						else{
							td2 = pow(1-error, N[j]-R[j]);
							one_err_R_Ns[N[j]-R[j]] = td2;
							tallys[N[j]-R[j]] |= 2;
						}

						if(tallys[R[j]] & 4){
							td3 = err_Rs[R[j]];
						}
						else{
							td3 = pow(error, R[j]) * pvv;
							err_Rs[R[j]] = td3;
							tallys[R[j]] |= 4;
						}

						if(tallys[N[j]-R[j]] & 8){
							td4 = err_R_Ns[N[j]-R[j]];
						}
						else{
							td4 =  pow(error, N[j]-R[j]) * prr;
							err_R_Ns[N[j]-R[j]] = td4;
							tallys[N[j]-R[j]] |= 8;
						}
// here the real calculations start
						if(hwe){
							if(tallys[N[j]] & 16){
								td5 = halfs[N[j]];
							}
							else{
								td5 = pow(0.5, N[j]-1) * p * (1-p);
								halfs[N[j]] = td5;
								tallys[N[j]] |= 16;
							}
							p_rr[j] = td1 * td4;
							p_rv[j] = td5;
							p_vv[j] = td2 * td3;
						}
						else{ // if(hwe)
							if(tallys[N[j]] & 16){
								td5 = halfs[N[j]];
							}
							else{
								td5 = pow(0.5, N[j]) * prv;
								halfs[N[j]] = td5;
								tallys[N[j]] |= 16;
							}
							p_rr[j] = td1 * td4;
							p_rv[j] = td5;	
							p_vv[j] = td2 * td3;
						} // if(hwe){} else {}

					}
					denom_sum = p_rr[j] + p_rv[j] + p_vv[j];
					if(para->debug & 2){
						logsum += log(denom_sum);
					}
					p_rr[j] /= denom_sum;
					p_rv[j] /= denom_sum;
					p_vv[j] /= denom_sum;
					expected_rr[j] = N[j] * p_rr[j];
					expected_vv[j] = N[j] * p_vv[j];
					expected_r_reads[j] = R[j] * p_rr[j];
					expected_v_reads[j] = R[j] * p_vv[j];
					sum_e += p_rr[j];
					sum_f += p_rv[j];
					sum_g += p_vv[j];
					sum_h += expected_rr[j];
					sum_i += expected_vv[j];
					sum_j += expected_r_reads[j];
					sum_k += expected_v_reads[j];
				} // if(indices[j] >= 0)
			} // for(j = 0; j < para->num_header_files; j++)

			if(fixerror == 0){
				prev_error = error;
				if((sum_h + sum_i) > para->em_threshold){
					error = (sum_h + sum_k - sum_j) / (sum_h + sum_i);
				}
				else{
					error = 0;
				}
			}
			prev_p = p;
			if(hwe){
				p = (2*sum_e + sum_f) / (2*valid_individuals);
				if((absf(p - prev_p) < para->em_threshold) && (absf(error - prev_error) < para->em_threshold)){
					continue_EM = 0;
				}
			}
			else{ // if(hwe)
				prev_prr = prr;
				prev_prv = prv;
				prev_pvv = pvv;
				prr = sum_e / valid_individuals;
				prv = sum_f / valid_individuals;
				pvv = sum_g / valid_individuals;
				if((absf(prr - prev_prr) < para->em_threshold) && (absf(prv - prev_prv) < para->em_threshold) && (absf(pvv - prev_pvv) < para->em_threshold) && (absf(error - prev_error) < para->em_threshold)){
					continue_EM = 0;
				}
				if(para->cov_matrix && (iter < min_iter)){				// to estimate the covariate matrix we need a minimum number of iterations, otherwise division by zero -> crash
					continue_EM = 1;
				}
			} // if(hwe){} else {}
			if(para->debug & 2){
				printf("\n### DEBUG: iteration= %i, read in masterlist= %i, valid individuals= %5.0f error= %7.5f logsum= %8.6f\n", iter, i, valid_individuals, error, logsum);
				if(hwe){
					printf("DEBUG: p(r)= %7.5f\n", p);
				}
				else{
					printf("DEBUG: p(rr)= %7.5f p(rv)= %7.5f p(vv)= %7.5f\n", prr,prv,pvv);
				}
				printf("sum_e (rr)= %8.5f =", sum_e);
				for(j = 0; j < para->num_header_files; j++){
					if(indices[j] >= 0){
						printf(" + %8.5f", p_rr[j]);
					}
					else{
						printf(" 0");
					}
				}
				printf("\n");
				printf("sum_f (rv)= %8.5f =", sum_f);
				for(j = 0; j < para->num_header_files; j++){
					if(indices[j] >= 0){
						printf(" + %8.5f", p_rv[j]);
					}
					else{
						printf(" 0");
					}
				}
				printf("\n");

				printf("sum_g (vv)= %8.5f =", sum_g);
				for(j = 0; j < para->num_header_files; j++){
					if(indices[j] >= 0){
						printf(" + %8.5f", p_vv[j]);
					}
					else{
						printf(" 0");
					}
				}
				printf("\n");

				printf("Total reads by individual:\n");
				for(j = 0; j < para->num_header_files; j++){
					printf("%5i",N[j]);
				}
				printf("\n");
				printf("Variant reads by individual:\n");
				for(j = 0; j < para->num_header_files; j++){
					printf("%5i",V[j]);
				}
				printf("\n");

				printf("sum_h = %8.5f =", sum_h);
				for(j = 0; j < para->num_header_files; j++){
					if(indices[j] >= 0){
						printf(" + %8.5f", expected_rr[j]);
					}
					else{
						printf(" 0");
					}
				}
				printf("\n");

				printf("sum_i = %8.5f =", sum_i);
				for(j = 0; j < para->num_header_files; j++){
					if(indices[j] >= 0){
						printf(" + %8.5f", expected_vv[j]);
					}
					else{
						printf(" 0");
					}
				}
				printf("\n");

				printf("sum_j = %8.5f =", sum_j);
				for(j = 0; j < para->num_header_files; j++){
					if(indices[j] >= 0){
						printf(" + %8.5f", expected_r_reads[j]);
					}
					else{
						printf(" 0");
					}
				}
				printf("\n");

				printf("sum_k = %8.5f =", sum_k);
				for(j = 0; j < para->num_header_files; j++){
					if(indices[j] >= 0){
						printf(" + %8.5f", expected_v_reads[j]);
					}
					else{
						printf(" 0");
					}
				}
				printf("\n");

			} // if(para->debug)
			if(para->cov_matrix){
				if(0 == hwe){										// this should always be 1, otherwise calculating the matrix is not very interesting
					para_estimates[iter][0] = prr;
					para_estimates[iter][1] = prv;
					para_estimates[iter][2] = pvv;
					para_estimates[iter][3] = error;
				} // if(0 == hwe)
			} // if(para->cov_matrix)
			iter++;
			
		} // while((continue_EM) && (iter < max_iter))
		delete [] one_err_R;
		delete [] one_err_R_N;
		delete [] err_R;
		delete [] err_R_N;
		delete [] half;
		delete [] tally;

		delete [] one_err_Rs;
		delete [] one_err_R_Ns;
		delete [] err_Rs;
		delete [] err_R_Ns;
		delete [] halfs;
		delete [] tallys;

		if(para->cov_matrix){
			if(0 == para->hwe){										// this should always be 1, otherwise calculating the matrix is not very interesting
				fprintf(fptr_matrix,"Matrices for %s start= %i stop= %i\n", masterlist[i].ccds_name, masterlist[i].start, masterlist[i].stop);
				cov_matrix_no_hwe((int)valid_individuals, iter, sum_h, sum_i);
			}
		} // if(para->cov_matrix)

		if(p >= 0){
			masterlist[i].error = (float)error;
			if(para->hwe){
				masterlist[i].p = p;
			}
			else{
				masterlist[i].p = (float)(prr + 0.5*prv);
				masterlist[i].prr = (float)prr;
				masterlist[i].prv = (float)prv;
				masterlist[i].pvv = (float)pvv;
			}
			masterlist[i].iterations = iter;
			masterlist[i].genotype = NULL;
			masterlist[i].individuals = para->num_header_files;
			masterlist[i].geno_prop = new float*[3];
			masterlist[i].genotype = new unsigned char[para->num_header_files];
			for(j = 0; j < 3; j++){
				masterlist[i].geno_prop[j] = new float[para->num_header_files];
			}
//		memcpy((void*)masterlist[i].geno_prop[0], (void*)p_rr, para->num_header_files*sizeof(float));
//		memcpy((void*)masterlist[i].geno_prop[1], (void*)p_rv, para->num_header_files*sizeof(float));
//		memcpy((void*)masterlist[i].geno_prop[2], (void*)p_vv, para->num_header_files*sizeof(float));
			sumrr = 0;
			sumrv = 0;
			sumvv = 0;
			sum_e = 0;
			masterlist[i].stat = -1.0;
			for(j = 0; j < para->num_header_files; j++){
				convert = (float)p_rr[j];
				if(convert < 0) convert = 0;
				masterlist[i].geno_prop[0][j] = convert;
				exprr= convert;
				convert = (float)p_rv[j];
				if(convert < 0) convert = 0;
				masterlist[i].geno_prop[1][j] = convert;
				exprv = convert;
				convert = (float)p_vv[j];
				if(convert < 0) convert = 0;
				masterlist[i].geno_prop[2][j] = convert;
				expvv = convert;

				if(indices[j] >= 0){
					masterlist[i].genotype[j] = (unsigned char)largest_of_threei(p_rr[j], p_rv[j], p_vv[j]);			// 0 == rr, 1 = rv, 2 == vv
					sumrr += exprr;
					sumrv += exprv;
					sumvv += expvv;
					sum_e += 1.0;
				}
				else{
					masterlist[i].genotype[j] = 255;
				}
			} // for(j = 0; j < para->num_header_files; j++)
			p = masterlist[i].p;
			if((sum_e > minchisq) && (p > 0.001) && (p < (0.999))){
				exprr = p*p*sum_e;
				exprv = 2*p*(1-p)*sum_e;
				expvv = (1-p)*(1-p)*sum_e;
				global = (exprr - sumrr) * (exprr - sumrr) / exprr;
				global += (exprv - sumrv) * (exprv - sumrv) / exprv;
				global += (expvv - sumvv) * (expvv - sumvv) / expvv;
				bound = 0;
				df = 1.0;
				cdfchi(&code, &p, &q, &global, &df, &status, &bound);
				masterlist[i].stat = q;
			}
		} // if(p >= 0)
		else{
			masterlist[i].error = (float)error;
			if(para->hwe){
				masterlist[i].p = -1.0;
			}
			else{
				masterlist[i].p = -1.0;
				masterlist[i].prr = -1.0;
				masterlist[i].prv = -1.0;
				masterlist[i].pvv = -1.0;
			}
			masterlist[i].iterations = 0;
			masterlist[i].genotype = NULL;
			masterlist[i].individuals = para->num_header_files;
			masterlist[i].geno_prop = new float*[3];
			masterlist[i].genotype = new unsigned char[para->num_header_files];
			for(j = 0; j < 3; j++){
				masterlist[i].geno_prop[j] = new float[para->num_header_files];
			}
//		memcpy((void*)masterlist[i].geno_prop[0], (void*)p_rr, para->num_header_files*sizeof(float));
//		memcpy((void*)masterlist[i].geno_prop[1], (void*)p_rv, para->num_header_files*sizeof(float));
//		memcpy((void*)masterlist[i].geno_prop[2], (void*)p_vv, para->num_header_files*sizeof(float));

			for(j = 0; j < para->num_header_files; j++){
				masterlist[i].geno_prop[0][j] = -1.0;
				masterlist[i].geno_prop[1][j] = -1.0;
				masterlist[i].geno_prop[2][j] = -1.0;
				masterlist[i].genotype[j] = 255;
			}
		}
		i++;
		if(i > stop){
			pthread_mutex_lock(&ana_syncher->sync_ana_mutex);
			if(ana_syncher->current < ana_syncher->items){
				start = ana_syncher->current;
				temp = (int)((ana_syncher->items - start) /(3*ana_syncher->threads));
				if(temp < min_slice){
					temp = min_slice;
				}
				stop = start + temp;
				if(stop < ana_syncher->items){
					ana_syncher->current = stop + 1;
				}
				else{
					stop = ana_syncher->items - 1;
					ana_syncher->current = ana_syncher->items;													// the next thread that sees "this" will exit. We're finished
				}
			}
			else{
				done = 1;																						// nothing else left to do. Exit EM algorithm
			}
			pthread_mutex_unlock(&ana_syncher->sync_ana_mutex);		// unlock mutex
			i = start;
		}
	} // for(i = start; i <= stop; i++)
//#ifdef PROFILE
//    PROFILE_STOP("EM_net");
//#endif
}

//***********************************


void Analyzer::cov_matrix_no_hwe(int valid_individuals, int iterations, double sum_h, double sum_i)
{
	// Parameters:
	// The first parameters from: CCDS* masterlist, ..., int valid_individuals
	// are of no further interest except that we need them to perform a single iteration via single_em_iteration_no_hwe(...) 
	// in order to calculate the estimate of the Jacobian matrix.
	// double **para_estimates[iterations][4] has the p(RR), p(RV), p(VV) and error for each iteration of the EM. 
	// The "iterations-1" index is the last valid entry and is the state of affairs right before it converged, ie. Theta-hat in the documentation.
	// Our first task is to find a reasonable starting point for our calculations. We need to assure we have a denominator >= para->min_denominator
	// to avoid numerical instbility due to division by small numbers.

	double input[4];
	double output[4];
	double mle[4];														// maximum likelihood estimate
	int rows[4] = {1,1,1,1};														// 1 = row not converged, 0 = row has converged
	int i,j,k;
	int do_row;
	int iter;
	double differnce = 1.0;
	double threshold = sqrt(para->em_threshold);						// use relaxed convergence criterion for matrix
//	double threshold = para->em_threshold;
	double temp;
	int not_converged = parameters;										// number of rows in matrix that haven't converged yet

	for(i = 0; i < parameters; i++){
		memset((void*)diff[i], 1, parameters*sizeof(double));
	}
	// find last iteration such that min(iteration, final_MLE) => para->min_denominator
	iter = 0;
	memcpy((void*)mle, (void*)para_estimates[iterations-1], parameters*sizeof(double));	// the value in the last valid entry of **para_estimates
	while((iter < (iterations-2-max_iteration_offset_for_matrix)) && (differnce >= para->min_denominator)){
		j = 0;
		for(j = 0; j < parameters; j++){
			temp = mle[j] - para_estimates[iter][j];
			if(temp < 0){
				temp *= -1;
			}
			if(temp < differnce){
				differnce = temp;
			}
		}
		iter++;
	}
	iter--;
	if(iter < 0){
		fprintf(fptr_matrix,"main::cov_matrix_no_hwe. Converged on one step. Insufficient data to estimate covariance matrix.\n");
		return;
	}
	// At this point we have the closest "acceptable" approximation in "estimates[4]", p(RR), p(RV), p(VV), error
	// We now need to iterate the DM-matrix as described in the algorithm documentation. A good first step is to find
	// the initial DM matrix. This is what we will do using the values in estimates[] and mle[]
	for(j = 0; j < parameters; j++){
		memcpy((void*)input, (void*)mle, parameters*sizeof(double));												// same as MLE values
		input[j] = para_estimates[iter][j];																					// ...except in the j_th position
		single_em_iteration_no_hwe(valid_individuals, input, output);	// run one iteration and use output[]
		for(k = 0; k < parameters; k++){
			current[j][k] = (output[k] - mle[k]) / (para_estimates[iter][j] - mle[j]);
		}
	}
	i = 0;																											// iteration count
	iter++;																											// index into **para_estimates. We can live of those up to & including "iterations-2" 
	// iteration start
	while(not_converged && (iter < (iterations-2))){																// so we can have at bare minimum "max_iteration_offset_for_matrix-1" iterations before we get into zero denomiators
		for(j = 0; j < parameters; j++){
			if(rows[j]){
				do_row = parameters;
				for(k = 0; k < parameters; k++){
					if(diff[j][k] < threshold){																		// see if this row already converged
						do_row--;	
					} // if(diff[j][k] < para->em_threshold)
				} // for(k = 0; k < parameters; k++)
				if(do_row){																							// if not, deal with is here
					memcpy((void*)input, (void*)mle, parameters*sizeof(double));									// same as MLE values
					input[j] = para_estimates[iter][j];
					single_em_iteration_no_hwe(valid_individuals, input, output);
					for(k = 0; k < parameters; k++){
						if(diff[j][k] >= threshold){																// only worry about entries that haven't converged yet
							prev[j][k] = current[j][k];
							current[j][k] = (output[k] - mle[k]) / (para_estimates[iter][j] - mle[j]);
							temp = current[j][k] - prev[j][k];
							if(temp < 0) temp *= -1;
							diff[j][k] = temp;
						} // if(diff[j][k] < para->em_threshold)
					} // for(k = 0; k < parameters; k++)
				} // if(do_row)
				else{
					rows[j] = 0;																					// row has converged, don't mess with it anymore
					not_converged--;
				} // else{}
			} // if(rows[j])
		} // for(j = 0; j < parameters; j++)
		iter++;
		i++;
	} // while(not_converged && (iter < (iterations-1)))
	Jac.init(parameters, current);																					// approximation of Jacobian, called "DM" in the documentation
	// now calculate the covariance matix C. We will recycle "double **current" for that purpose
	for(i = 0; i < parameters; i++){
		memset((void*)current[i], 0, parameters*sizeof(double));													// set everything to zero for starters
	}
	current[parameters-1][parameters-1] = (mle[parameters-1] * (1 - mle[parameters-1])) / (sum_h + sum_i);			// ...see documentation of covaraince matrix for details
	current[0][0] = (mle[0]*(1-mle[0])) / (double)valid_individuals;
	current[1][1] = (mle[1]*(1-mle[1])) / (double)valid_individuals;
	current[2][2] = (mle[2]*(1-mle[2])) / (double)valid_individuals;
	current[1][0] = (-mle[0]*mle[1]) / (double)valid_individuals;
	current[2][0] = (-mle[0]*mle[2]) / (double)valid_individuals;
	current[2][1] = (-mle[1]*mle[2]) / (double)valid_individuals;
	current[0][1] = current[1][0];																					// symmetric matrix
	current[0][2] = current[2][0];
	current[1][2] = current[2][1];
	Cov.init(parameters, current);
	for(i = 0; i < parameters; i++){
		memset((void*)current[i], 0, parameters*sizeof(double));													// set everything to zero for starters
		current[i][i] = 1.0;																						// identity matrix
	}
	IDM.init(parameters, current);
	Vm.init(parameters, current);
	Tm.init(parameters, current);
	

	Dif.init(parameters, diff);

	// start debug
	fprintf(fptr_matrix, "Difference from convergence\n");
	Dif.print_matrix_to_file(fptr_matrix);
	fprintf(fptr_matrix, "Jacobian\n");
	Jac.print_matrix_to_file(fptr_matrix);
	fprintf(fptr_matrix,"Covariance\n");
	Cov.print_matrix_to_file(fptr_matrix);
	IDM -= Jac;
	double **inv = IDM.find_inverse();
	Inv.init(parameters, inv);
//	fprintf(fptr,"(I-Jac)^-1\n");
//	Inv.print_matrix_to_file(fptr);
	Vm = Cov;
	Vm *= Inv;
	fprintf(fptr_matrix, "Obs. Covariance\n");
	Vm.print_matrix_to_file(fptr_matrix);
	Vm = Cov;
	Vm *= Jac;
	Vm *= Inv;
	fprintf(fptr_matrix,"Inflation\n");
	Vm.print_matrix_to_file(fptr_matrix);
	// end debug
	// delete memory

}

//***********************************


void Analyzer::single_em_iteration_no_hwe(int valid_individuals, double input[], double output[] )
{
	double sum_e, sum_g, sum_f, sum_k, sum_h, sum_j, sum_i;
	int j;
	unsigned long long bin_coeff;
	int overflow;
	double prr = input[0];
	double prv = input[1];
	double pvv = input[2]; 
	double error = input[3];
	double denom_sum;
	double bc;
	double *p_rr = new double[para->num_header_files];					// genotype "r,r" = reference, refevernce
	double *p_rv = new double[para->num_header_files];						// genotype "v,r" = variant, reference
	double *p_vv = new double[para->num_header_files];						// genotype "v,v," = variant, variant
	double *expected_rr = new double[para->num_header_files];				// expected number of RR reads
	double *expected_vv = new double[para->num_header_files];				// expected number of VV reads
	double *expected_r_reads = new double[para->num_header_files];			// expected r-reads from rr genotype
	double *expected_v_reads = new double[para->num_header_files];			// expected v-reads from vv genotypes

	sum_e = sum_g = sum_f = sum_k = sum_h = sum_j = sum_i = 0;
	for(j = 0; j < para->num_header_files; j++){
		if(indices[j] >= 0){								// only those individuals have valid data
			bin_coeff = n_choose_k(N[j], R[j], overflow);
			bc = (double)bin_coeff;
			if(overflow){
				fprintf(stderr,"main::analysis3. Overflow in calculation of binomial coefficient n= %i, k = %i. Bye\n",N[j], R[j]);
				exit(0);
			} // if(overflow)
			p_rr[j] = bc * pow(1-error, R[j]) * pow(error, N[j]-R[j]) * prr;
			p_rv[j] = bc * pow(0.5, N[j]) * prv;	
			p_vv[j] = bc * pow(1-error, N[j]-R[j]) * pow(error, R[j]) * pvv;
			denom_sum = p_rr[j] + p_rv[j] + p_vv[j];
			p_rr[j] /= denom_sum;
			p_rv[j] /= denom_sum;
			p_vv[j] /= denom_sum;
			expected_rr[j] = N[j] * p_rr[j];
			expected_vv[j] = N[j] * p_vv[j];
			expected_r_reads[j] = R[j] * p_rr[j];
			expected_v_reads[j] = R[j] * p_vv[j];
			sum_e += p_rr[j];
			sum_f += p_rv[j];
			sum_g += p_vv[j];
			sum_h += expected_rr[j];
			sum_i += expected_vv[j];
			sum_j += expected_r_reads[j];
			sum_k += expected_v_reads[j];
		} // if(indices[j] >= 0)
	} // for(j = 0; j < para->num_header_files; j++)
	output[3] = (sum_h + sum_k - sum_j) / (sum_h + sum_i);						// error
	output[0] = sum_e / valid_individuals;										// prr
	output[1] = sum_f / valid_individuals;										// prv
	output[2] = sum_g / valid_individuals;										// pvv


}

//*************************************


double Analyzer::absf(const double &val)
{
	double result = val;
	if(result < 0) result *= -1;
	return(result);
}


//***********************************


double Analyzer::powexpo(double a, double b)
{
	double result = -infinity;
	if(a > 0){
		result = b*log(a);
	}
	return(result);
}

//***********************************

double Analyzer::largest_of_three(double a, double b, double c)
{
	double result;
	if( a > b){
		result = a;
	}
	else{
		result = b;
	}
	if(result == a){
		if(c > a){
			result = c;
		}
	}
	else{
		if(c > b){
			result = c;
		}
	}
	return(result);
}

//***********************************

double Analyzer::diff_of_three(double largest, double a, double b, double c)
{
	double result = largest - a;
	double temp = largest - b;
	if(temp > result){
		result = temp;
	}
	temp = largest - c;
	if(temp > result){
		result = temp;
	}
	return(result);
}

//***********************************

int Analyzer::largest_of_threei(double a, double b, double c)
{
	int result = 0;

	if(a > b){
		result = 0;
		if(c > a) result = 2;
	}
	else{
		result = 1;
		if(c > b) result = 2;
	}
	return(result);
}

//***********************************
